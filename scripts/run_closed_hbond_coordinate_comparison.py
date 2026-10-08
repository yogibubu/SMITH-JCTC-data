#!/usr/bin/env python3
"""Compare closed double-H-bond SONIC charts on one frozen geometry.

The comparison separates coordinate construction from optimizer behavior.  It
builds the current fragment-pose and augmented-graph pseudobond definitions,
then constructs two diagnostic hybrid bases from their analytic Wilson rows:

* one H-bond distance plus five residual fragment-pose directions;
* two H-bond distances plus four residual fragment-pose directions.

The hybrids are diagnostic linear tangent bases, not frozen SMITH coordinate
contracts.  They identify whether a contact-anchored fragment chart merits a
native nonlinear implementation before optimizer benchmarks are attempted.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np

from _matrix_checkout import add_matrix_packages_to_path, matrix_root


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GEOMETRY = ROOT / "data" / "formamide_water_double_hbond.xyz"
DEFAULT_OUTPUT = ROOT / "data" / "formamide_water_coordinate_comparison.json"
TOL = 1.0e-10


@dataclass(frozen=True)
class BasisMetrics:
    label: str
    status: str
    rows: int
    rank: int
    target_rank: int
    normalized_condition_number: float
    relative_block_condition_number: float
    b_nnz: int
    b_density: float
    mean_atoms_per_row: float
    maximum_atoms_per_row: int
    projector_distance_from_fragment_pose: float
    coordinate_families: dict[str, int]
    selected_pose_coordinates: tuple[str, ...] = ()
    selected_contact_coordinates: tuple[str, ...] = ()


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_head(path: Path) -> str:
    completed = subprocess.run(
        ("git", "rev-parse", "HEAD"),
        cwd=path,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def orthonormal_row_basis(rows: np.ndarray, tolerance: float = TOL) -> np.ndarray:
    _u, singular_values, vh = np.linalg.svd(np.asarray(rows, dtype=float), full_matrices=False)
    rank = int(np.sum(singular_values > tolerance))
    return vh[:rank]


def normalized_rows(rows: np.ndarray) -> np.ndarray:
    dense = np.asarray(rows, dtype=float)
    norms = np.linalg.norm(dense, axis=1)
    if np.any(norms <= TOL):
        raise RuntimeError("coordinate basis contains a numerically null Wilson row")
    return dense / norms[:, None]


def nonzero_condition(rows: np.ndarray, tolerance: float = TOL) -> float:
    singular_values = np.linalg.svd(np.asarray(rows, dtype=float), compute_uv=False)
    nonzero = singular_values[singular_values > tolerance]
    if not len(nonzero):
        return float("inf")
    return float(nonzero[0] / nonzero[-1])


def relative_condition(rows: np.ndarray, intramolecular_rows: np.ndarray) -> float:
    dense = normalized_rows(rows)
    intra_basis = orthonormal_row_basis(normalized_rows(intramolecular_rows))
    residual = dense - (dense @ intra_basis.T) @ intra_basis
    return nonzero_condition(residual)


def projector_distance(first: np.ndarray, second: np.ndarray) -> float:
    first_basis = orthonormal_row_basis(first)
    second_basis = orthonormal_row_basis(second)
    first_projector = first_basis.T @ first_basis
    second_projector = second_basis.T @ second_basis
    return float(np.linalg.norm(first_projector - second_projector, ord=2))


def atom_support_metrics(rows: np.ndarray) -> tuple[int, float, int]:
    dense = np.asarray(rows, dtype=float)
    cartesian_support = np.abs(dense) > 1.0e-12
    atom_support = cartesian_support.reshape(dense.shape[0], -1, 3).any(axis=2)
    atoms_per_row = atom_support.sum(axis=1)
    return (
        int(cartesian_support.sum()),
        float(np.mean(atoms_per_row)),
        int(np.max(atoms_per_row)),
    )


def basis_metrics(
    *,
    label: str,
    rows: np.ndarray,
    target_rank: int,
    intramolecular_rows: np.ndarray,
    reference_rows: np.ndarray,
    families: Counter[str],
    selected_pose: tuple[str, ...] = (),
    selected_contacts: tuple[str, ...] = (),
    status: str = "implemented",
) -> BasisMetrics:
    dense = np.asarray(rows, dtype=float)
    rank = int(np.linalg.matrix_rank(dense, tol=TOL))
    nnz, mean_atoms, maximum_atoms = atom_support_metrics(dense)
    return BasisMetrics(
        label=label,
        status=status,
        rows=int(dense.shape[0]),
        rank=rank,
        target_rank=target_rank,
        normalized_condition_number=nonzero_condition(normalized_rows(dense)),
        relative_block_condition_number=relative_condition(dense, intramolecular_rows),
        b_nnz=nnz,
        b_density=float(nnz / dense.size),
        mean_atoms_per_row=mean_atoms,
        maximum_atoms_per_row=maximum_atoms,
        projector_distance_from_fragment_pose=projector_distance(dense, reference_rows),
        coordinate_families=dict(sorted(families.items())),
        selected_pose_coordinates=selected_pose,
        selected_contact_coordinates=selected_contacts,
    )


def best_pose_completion(
    *,
    fixed_rows: np.ndarray,
    pose_rows: np.ndarray,
    pose_names: tuple[str, ...],
    number_to_select: int,
    target_rank: int,
) -> tuple[np.ndarray, tuple[str, ...]]:
    candidates: list[tuple[float, float, tuple[int, ...], np.ndarray]] = []
    for indices in combinations(range(len(pose_rows)), number_to_select):
        selected = pose_rows[np.asarray(indices, dtype=int)]
        rows = np.vstack((fixed_rows, selected))
        if np.linalg.matrix_rank(rows, tol=TOL) != target_rank:
            continue
        candidates.append(
            (
                nonzero_condition(normalized_rows(rows)),
                -float(np.linalg.det(normalized_rows(rows) @ normalized_rows(rows).T)),
                indices,
                rows,
            )
        )
    if not candidates:
        raise RuntimeError("fragment-pose candidates cannot complete the requested hybrid rank")
    _condition, _negative_volume, indices, rows = min(candidates, key=lambda item: item[:2])
    return rows, tuple(pose_names[index] for index in indices)


def distance(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.linalg.norm(first - second))


def angle(first: np.ndarray, center: np.ndarray, third: np.ndarray) -> float:
    left = first - center
    right = third - center
    cosine = float(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    matrix_checkout = add_matrix_packages_to_path()
    from matrix_chem import (
        preprocess_to_enriched_xyz,
        read_enriched_xyz,
        write_validation_section,
    )
    from matrix_fragments import write_fragment_build_section
    from matrix_smith import build_gic_b_matrix, write_gicforge_build_sections

    source = args.geometry.resolve()
    with TemporaryDirectory(prefix="smith-double-hbond-") as temporary:
        temporary_path = Path(temporary)
        definitions = {}
        matrices = {}
        for label, fragment_mode in (
            ("fragment_pose", "special-coordinates"),
            ("double_pseudobond", "pseudo-bonds"),
        ):
            xyzin = temporary_path / f"{label}.xyzin"
            preprocess_to_enriched_xyz(source, xyzin)
            write_validation_section(xyzin)
            write_fragment_build_section(xyzin)
            definition = write_gicforge_build_sections(
                xyzin,
                symmetrize=False,
                fragment_mode=fragment_mode,
            )
            geometry = read_enriched_xyz(xyzin)
            matrix = build_gic_b_matrix(
                definition,
                coordinates_angstrom=geometry.coordinates_angstrom,
            )
            definitions[label] = definition
            matrices[label] = matrix

    pose_definition = definitions["fragment_pose"]
    pseudo_definition = definitions["double_pseudobond"]
    pose_matrix = matrices["fragment_pose"]
    pseudo_matrix = matrices["double_pseudobond"]
    pose_rows_all = np.asarray(pose_matrix.rows, dtype=float)
    pseudo_rows_all = np.asarray(pseudo_matrix.rows, dtype=float)
    target_rank = pose_definition.target_rank
    if pseudo_definition.target_rank != target_rank:
        raise RuntimeError("fragment-pose and pseudobond definitions disagree on target rank")

    fragment_families = {"FRAG_TRANSLATION", "FRAG_ORIENTATION"}
    intra_indices = tuple(
        index
        for index, coordinate in enumerate(pose_definition.gics)
        if coordinate.family not in fragment_families
    )
    pose_indices = tuple(
        index
        for index, coordinate in enumerate(pose_definition.gics)
        if coordinate.family in fragment_families
    )
    contact_indices = tuple(
        index
        for index, coordinate in enumerate(pseudo_definition.gics)
        if coordinate.family == "HBOND_DISTANCE"
    )
    if len(intra_indices) != target_rank - 6 or len(pose_indices) != 6:
        raise RuntimeError("unexpected fragment-pose partition")
    if len(contact_indices) != 2:
        raise RuntimeError(
            f"expected two selected H-bond distances, found {len(contact_indices)}"
        )

    intra_rows = pose_rows_all[np.asarray(intra_indices, dtype=int)]
    pose_rows = pose_rows_all[np.asarray(pose_indices, dtype=int)]
    contact_rows = pseudo_rows_all[np.asarray(contact_indices, dtype=int)]
    pose_names = tuple(pose_definition.gics[index].name for index in pose_indices)
    contact_names = tuple(pseudo_definition.gics[index].name for index in contact_indices)

    one_contact_rows, one_contact_pose = best_pose_completion(
        fixed_rows=np.vstack((intra_rows, contact_rows[:1])),
        pose_rows=pose_rows,
        pose_names=pose_names,
        number_to_select=5,
        target_rank=target_rank,
    )
    two_contact_rows, two_contact_pose = best_pose_completion(
        fixed_rows=np.vstack((intra_rows, contact_rows)),
        pose_rows=pose_rows,
        pose_names=pose_names,
        number_to_select=4,
        target_rank=target_rank,
    )

    metrics = (
        basis_metrics(
            label="fragment_pose_exponential_mapping",
            rows=pose_rows_all,
            target_rank=target_rank,
            intramolecular_rows=intra_rows,
            reference_rows=pose_rows_all,
            families=Counter(coordinate.family for coordinate in pose_definition.gics),
        ),
        basis_metrics(
            label="one_contact_plus_residual_pose",
            rows=one_contact_rows,
            target_rank=target_rank,
            intramolecular_rows=intra_rows,
            reference_rows=pose_rows_all,
            families=Counter(
                {
                    **Counter(pose_definition.gics[index].family for index in intra_indices),
                    "HBOND_DISTANCE": 1,
                    "RESIDUAL_FRAGMENT_POSE": 5,
                }
            ),
            selected_pose=one_contact_pose,
            selected_contacts=contact_names[:1],
            status="diagnostic tangent basis; nonlinear chart not yet implemented",
        ),
        basis_metrics(
            label="two_contacts_plus_residual_pose",
            rows=two_contact_rows,
            target_rank=target_rank,
            intramolecular_rows=intra_rows,
            reference_rows=pose_rows_all,
            families=Counter(
                {
                    **Counter(pose_definition.gics[index].family for index in intra_indices),
                    "HBOND_DISTANCE": 2,
                    "RESIDUAL_FRAGMENT_POSE": 4,
                }
            ),
            selected_pose=two_contact_pose,
            selected_contacts=contact_names,
            status="diagnostic tangent basis; nonlinear chart not yet implemented",
        ),
        basis_metrics(
            label="current_double_pseudobond_augmented_graph",
            rows=pseudo_rows_all,
            target_rank=target_rank,
            intramolecular_rows=intra_rows,
            reference_rows=pose_rows_all,
            families=Counter(coordinate.family for coordinate in pseudo_definition.gics),
        ),
    )

    coords = np.asarray(pose_definition.reference_coordinates_angstrom, dtype=float)
    contacts = (
        {
            "donor": "N1-H3",
            "acceptor": "O7(water)",
            "h_acceptor_distance_angstrom": distance(coords[2], coords[6]),
            "donor_h_acceptor_angle_degrees": angle(coords[0], coords[2], coords[6]),
        },
        {
            "donor": "O7-H9(water)",
            "acceptor": "O6(formamide)",
            "h_acceptor_distance_angstrom": distance(coords[8], coords[5]),
            "donor_h_acceptor_angle_degrees": angle(coords[6], coords[8], coords[5]),
        },
    )
    payload = {
        "schema": "smith.closed_hbond_coordinate_comparison.v1",
        "system": "formamide-water cyclic double hydrogen bond",
        "geometry": str(source.relative_to(ROOT)),
        "geometry_sha256": file_sha256(source),
        "matrix_checkout": str(matrix_checkout),
        "matrix_git_head": git_head(matrix_root()),
        "symmetrized": False,
        "point_group": pose_definition.point_group,
        "target_rank": target_rank,
        "contacts": contacts,
        "detected_pseudo_bonds": [list(pair) for pair in pseudo_definition.pseudo_bonds],
        "detected_pseudo_bond_kinds": list(pseudo_definition.pseudo_bond_kinds),
        "pseudo_cycle_coordinate_families_present": any(
            coordinate.family.startswith("PSEUDO_CYCLE")
            for coordinate in pseudo_definition.gics
        ),
        "interpretation": (
            "All four bases close the same 21-dimensional vibrational tangent space. "
            "The two hybrid bases are linear diagnostics only.  Optimization evidence "
            "must be obtained after their nonlinear realization is implemented."
        ),
        "bases": [asdict(item) for item in metrics],
        "independent_optimizer_evidence": {
            "source_repository": "../LINK",
            "same_starting_geometry_within_each_level": True,
            "GFN2-xTB": {
                "pseudobond_steps": 5,
                "exponential_mapping_steps": 5,
                "source": "data/xtb_campaign/intermolecular_final_62c6941b_p210.csv",
            },
            "ORCA_HF_STO-3G": {
                "pseudobond_steps": 31,
                "pseudobond_energy_gradient_evaluations": 88,
                "exponential_mapping_steps": 26,
                "exponential_mapping_energy_gradient_evaluations": 59,
                "source": "manuscript/link_optimizer.tex, table orca-link-intermolecular",
            },
            "PySCF_PBE0_def2-SVP_D3BJ_ATM": {
                "pseudobond_steps": 14,
                "pseudobond_energy_gradient_evaluations": 15,
                "exponential_mapping_steps": 7,
                "exponential_mapping_energy_gradient_evaluations": 8,
                "source": "data/pyscf_pbe0_def2svp_d3bj_ladder2_20260818.json",
            },
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
