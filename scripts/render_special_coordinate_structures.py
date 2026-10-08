#!/usr/bin/env python3
"""Render real molecular geometries for the special-coordinate figure."""
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "special_coordinate_structures.png"
OUT_PDF = ROOT / "figures" / "special_coordinate_structures.pdf"

COLORS = {"H": "#f5f5f5", "C": "#596b78", "N": "#315a9b",
          "O": "#c74d4d", "Fe": "#ef9c27", "Pd": "#6da9a2",
          "B": "#d69b54", "Cl": "#58a65c"}
RADII = {"H": 0.23, "C": 0.38, "N": 0.36, "O": 0.34,
         "Fe": 0.52, "Pd": 0.56, "B": 0.34, "Cl": 0.43}

def parse_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    atoms, xyz = [], []
    for line in lines[2:2+n]:
        s = line.split()
        atoms.append(s[0])
        xyz.append([float(v) for v in s[1:4]])
    return atoms, np.asarray(xyz, float)

def parse_gjf(path):
    atoms, xyz = [], []
    for line in Path(path).read_text().splitlines():
        s = line.split()
        if len(s) >= 4 and re.fullmatch(r"[A-Z][a-z]?", s[0]):
            try:
                xyz.append([float(s[1]), float(s[2]), float(s[3])])
                atoms.append(s[0])
            except ValueError:
                pass
    return atoms, np.asarray(xyz, float)

def bonds(atoms, xyz):
    out = []
    for i in range(len(atoms)):
        for j in range(i):
            d = np.linalg.norm(xyz[i] - xyz[j])
            cutoff = 1.25 * (RADII.get(atoms[i], .35) + RADII.get(atoms[j], .35)) * 1.7
            if d > 0.55 and d < cutoff:
                out.append((j, i))
    return out

def draw_molecule(ax, atoms, xyz, special=(), arrows=(), plane=None):
    centered = xyz - xyz.mean(axis=0)
    if plane is None:
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        xy = centered @ vt[:2].T
    else:
        xy = centered[:, list(plane)]
    xy[:, 1] *= -1.0
    segs = [(xy[i], xy[j]) for i, j in bonds(atoms, xyz)]
    for p, q in segs:
        ax.plot([p[0], q[0]], [p[1], q[1]], color="#81919a",
                lw=2.0, alpha=.9, zorder=1)
    for i, j in special:
        ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]],
                color="#16a6aa", lw=2.0,
                ls=(0, (5, 4)), alpha=.95)
    for i, (a, p) in enumerate(zip(atoms, xy)):
        ax.scatter([p[0]], [p[1]], s=(RADII.get(a, .35)*46)**2,
                   color=COLORS.get(a, "#777"), edgecolor="#ffffff",
                   linewidth=.65, zorder=4)
    for start, vec, label in arrows:
        p = xy[start]
        v = np.asarray(vec)[:2]
        ax.arrow(p[0], p[1], v[0], v[1], color="#7046a8",
                 width=.015, head_width=.18, length_includes_head=True)
        q = p + v*1.12
        ax.text(q[0], q[1], label, color="#7046a8", fontsize=8, ha="center")
    span = np.ptp(xy, axis=0).max()
    center = xy.mean(axis=0)
    ax.set_xlim(center[0]-span*.62, center[0]+span*.62)
    ax.set_ylim(center[1]-span*.62, center[1]+span*.62)
    ax.set_aspect("equal")
    ax.set_axis_off()

def panel(fig, rect, letter, title, subtitle, atoms, xyz, special=(), arrows=(), plane=None):
    ax = fig.add_axes(rect)
    ax.patch.set_facecolor("#f5f7f8")
    for spine in ax.spines.values():
        spine.set_visible(False)
    draw_molecule(ax, atoms, xyz, special, arrows, plane)
    ax.text(.03, .94, letter, transform=ax.transAxes, fontsize=13,
              color="white", weight="bold", ha="center", va="center",
              bbox=dict(boxstyle="round,pad=.22", fc="#173b59", ec="none"))
    ax.text(.10, .945, title, transform=ax.transAxes, fontsize=12,
              color="#173b59", weight="bold", va="center")
    ax.text(.10, .885, subtitle, transform=ax.transAxes, fontsize=8.5,
              color="#5e7180", va="center")
    return ax

def main():
    f = ROOT / "data" / "additional_gaussian_inputs"
    d = ROOT / "data" / "noncovalent_gfn2xtb_20260824"
    atoms_fe, xyz_fe = parse_xyz(f / "ferrocene_D5h.xyzin")
    atoms_fw, xyz_fw = parse_xyz(d / "formic_acid_water.xyz")
    atoms_wd, xyz_wd = parse_xyz(d / "water_dimer.xyz")
    atoms_bw, xyz_bw = parse_gjf(d / "benzene_water" / "benzene_water_corrected_native.gjf")
    atoms_ap, xyz_ap = parse_gjf(ROOT / "si" / "eta3_allyl_palladium_g16.gjf")
    atoms_ff, xyz_ff = parse_xyz(ROOT / "data" / "formamide_water_double_hbond.xyz")

    fig = plt.figure(figsize=(12.0, 8.2), facecolor="white")
    fig.text(.5, .975, "Special coordinates remain attached to physical molecular structures",
             ha="center", va="top", fontsize=18, color="#173b59", weight="bold")
    rects = [(0.025,.525,.30,.37),(0.35,.525,.30,.37),(0.675,.525,.30,.37),
             (0.025,.075,.30,.37),(0.35,.075,.30,.37),(0.675,.075,.30,.37)]
    panel(fig, rects[0], "a", "Ferrocene", "ring centre / Fe distance",
          atoms_fe, xyz_fe, special=((1, 0),), plane=(0, 1))
    panel(fig, rects[1], "b", "Formic acid--water", "fragment pose / H-bond contacts",
          atoms_fw, xyz_fw, special=((2, 7), (4, 5)))
    panel(fig, rects[2], "c", "Water dimer", "relative translation and orientation",
          atoms_wd, xyz_wd, special=((2, 4),))
    panel(fig, rects[3], "d", "Benzene--water", "ring centre / fragment pose",
          atoms_bw, xyz_bw, special=((6, 12),), plane=(1, 2))
    panel(fig, rects[4], "e", r"$\eta^3$ allyl--PdCl", "ligand-centre coordinate",
          atoms_ap, xyz_ap, special=((0, 2),))
    panel(fig, rects[5], "f", "Formamide--water", "closed fragment graph / relative pose",
          atoms_ff, xyz_ff, special=((1, 6), (7, 0)))
    fig.savefig(OUT, dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")

if __name__ == "__main__":
    main()
