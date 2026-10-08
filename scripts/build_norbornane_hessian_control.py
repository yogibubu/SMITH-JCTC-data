#!/usr/bin/env python3
"""Compare four rank-complete Hessian representations for norbornane.

The stationary Cartesian geometry and analytic HF/STO-3G Hessian were
computed with PySCF 2.14.0. The fixed R/A/D primitive definitions are the
native redundant list for the same norbornane atom ordering. Set MATRIX_ROOT
to the local MATRIX checkout before running this script.
"""

from __future__ import annotations

import json
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from geometric.internal import Angle, Dihedral, Distance

from _matrix_checkout import add_matrix_packages_to_path
from build_application_figures import (
    block_coupling_metrics,
    coupling_matrix,
    coupling_metrics,
    internal_hessian,
    vibrational_cartesian_hessian,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "hessian_norbornane"
BOHR_TO_ANGSTROM = 0.529177210903


def main() -> None:
    add_matrix_packages_to_path()
    from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
    from matrix_smith import build_gic_b_matrix, write_gicforge_build_sections

    xyz = DATA / "norbornane_hfsto3g_stationary.xyz"
    lines = xyz.read_text().splitlines()
    n = int(lines[0])
    coords = np.array([[float(v) for v in line.split()[1:4]] for line in lines[2 : 2 + n]])
    raw_hessian = np.load(DATA / "norbornane_hfsto3g_hessian.npy")
    cart_hessian = raw_hessian.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    cart_hessian = (cart_hessian + cart_hessian.T) / 2

    primitives = json.loads((DATA / "native_primitives.json").read_text())
    coordinate_types = {"R": Distance, "A": Angle, "D": Dihedral}
    coords_bohr = coords / BOHR_TO_ANGSTROM
    b_primitive = np.asarray([
        coordinate_types[item["kind"]](*[atom - 1 for atom in item["atoms"]])
        .derivative(coords_bohr).ravel()
        for item in primitives
    ])
    eigenvalues, eigenvectors = np.linalg.eigh(b_primitive @ b_primitive.T)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    rank = int(np.count_nonzero(eigenvalues > 1e-10 * eigenvalues[0]))
    assert rank == 3 * n - 6 == 51
    b_delocalized = eigenvectors[:, :rank].T @ b_primitive
    assert np.allclose(b_delocalized @ b_delocalized.T,
                       np.diag(eigenvalues[:rank]), rtol=1e-9, atol=1e-9)

    # This stationary geometry is oriented so that the C2v operations flip x
    # and y. Verify the BB^T eigenvectors belong to pure symmetry subspaces.
    cartesian_operators = []
    max_operation_deviation = 0.0
    for signs in ((1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1)):
        rotation = np.diag(signs)
        distances = np.linalg.norm((coords @ rotation.T)[:, None, :] - coords[None, :, :], axis=2)
        permutation = distances.argmin(axis=1)
        assert len(set(permutation)) == n
        max_operation_deviation = max(max_operation_deviation,
                                      float(distances[np.arange(n), permutation].max()))
        operator = np.zeros((3 * n, 3 * n))
        for source, target in enumerate(permutation):
            operator[3 * target : 3 * target + 3, 3 * source : 3 * source + 3] = rotation
        cartesian_operators.append(operator)
    assert max_operation_deviation < 1e-8
    normalized = b_delocalized / np.linalg.norm(b_delocalized, axis=1)[:, None]
    characters = np.array([[1, a, b, a * b] for a in (1, -1) for b in (1, -1)])
    fractions = np.stack([
        np.einsum("ij,jk,ik->i", normalized,
                  sum(c * op for c, op in zip(character, cartesian_operators)) / 4,
                  normalized)
        for character in characters
    ], axis=1)
    minimum_irrep_purity = float(fractions.max(axis=1).min())
    assert minimum_irrep_purity > 0.9999
    delocalized_irreps = fractions.argmax(axis=1)
    delocalized_order = np.argsort(delocalized_irreps, kind="stable")
    b_delocalized = b_delocalized[delocalized_order]
    delocalized_counts = np.bincount(delocalized_irreps, minlength=4)

    b_sonic = {}
    labels = None
    with tempfile.TemporaryDirectory(prefix="norbornane-hessian-") as scratch:
        for name, symmetrize in (("local ESALC-I", False), ("symmetry-adapted ESALC-I", True)):
            xyzin = Path(scratch) / ("sym.xyzin" if symmetrize else "local.xyzin")
            imported = preprocess_to_enriched_xyz(xyz, xyzin)
            assert imported.proposed_point_group == "C2v"
            write_validation_section(xyzin)
            definition = write_gicforge_build_sections(xyzin, symmetrize=symmetrize)
            result = build_gic_b_matrix(definition, coordinates_angstrom=coords)
            b_sonic[name] = np.asarray(result.rows, dtype=float)
            assert b_sonic[name].shape == (rank, 3 * n)
            assert np.linalg.matrix_rank(b_sonic[name]) == rank
            if symmetrize:
                labels = tuple(result.irreps)
    irrep_order = {"A1": 0, "A2": 1, "B1": 2, "B2": 3}
    sonic_order = np.argsort([irrep_order[label] for label in labels], kind="stable")
    b_sonic["symmetry-adapted ESALC-I"] = b_sonic["symmetry-adapted ESALC-I"][sonic_order]
    labels = tuple(labels[index] for index in sonic_order)
    sonic_counts = [labels.count(label) for label in irrep_order]

    forces = {
        "Cartesian vibrational": vibrational_cartesian_hessian(cart_hessian, coords),
        "delocalized internals": internal_hessian(cart_hessian, b_delocalized),
        **{name: internal_hessian(cart_hessian, b) for name, b in b_sonic.items()},
    }
    maps = {name: coupling_matrix(force) for name, force in forces.items()}
    metrics = {name: coupling_metrics(c) for name, c in maps.items()}
    report = {
        "system": "norbornane",
        "method": "PySCF 2.14.0 RHF/STO-3G analytic Cartesian Hessian",
        "optimized_energy_hartree": -268.89382384573236,
        "final_cartesian_gradient_max_hartree_per_bohr": 6.1891e-6,
        "final_cartesian_gradient_rms_hartree_per_bohr": 2.5413e-6,
        "geometry_point_group": "C2v",
        "lowest_vibrational_cartesian_hessian_eigenvalue": float(
            np.linalg.eigvalsh(forces["Cartesian vibrational"])[0]
        ),
        "primitive_count": len(primitives),
        "primitive_counts": dict(Counter(item["kind"] for item in primitives)),
        "primitive_rank": rank,
        "delocalized_relative_eigenvalue_threshold": 1e-10,
        "delocalized_largest_discarded_eigenvalue": float(eigenvalues[rank]),
        "delocalized_minimum_irrep_purity": minimum_irrep_purity,
        "max_symmetry_operation_deviation_angstrom": max_operation_deviation,
        "sonic_rank": rank,
        "symmetry_irrep_counts": dict(Counter(labels)),
        "symmetry_block_metrics": block_coupling_metrics(maps["symmetry-adapted ESALC-I"], labels),
        "metrics": metrics,
    }
    (DATA / "analysis.json").write_text(json.dumps(report, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.7), constrained_layout=True)
    boundaries = {
        "delocalized internals": np.cumsum(delocalized_counts)[:-1],
        "symmetry-adapted ESALC-I": np.cumsum(sonic_counts)[:-1],
    }
    for letter, (ax, (name, matrix)) in zip("abcd", zip(axes.flat, maps.items())):
        im = ax.imshow(matrix, origin="lower", cmap="magma", vmin=0, vmax=0.7,
                       interpolation="nearest", rasterized=True)
        ax.set_title(f"({letter}) {name}\nmean = {metrics[name]['mean_off_diagonal']:.3f}",
                     fontsize=10, loc="left")
        ax.set_xlabel("Coordinate index")
        ax.set_ylabel("Coordinate index")
        for boundary in boundaries.get(name, ()):
            position = boundary - 0.5
            ax.axhline(position, color="white", linewidth=0.5, linestyle=":")
            ax.axvline(position, color="white", linewidth=0.5, linestyle=":")
    fig.colorbar(im, ax=axes, shrink=0.75, label=r"$C_{ij}$")
    fig.savefig(ROOT / "figures" / "norbornane_hessian_heatmap.pdf")
    fig.savefig(ROOT / "figures" / "norbornane_hessian_heatmap.png", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
