#!/usr/bin/env python3
"""Build a clean schematic atlas for fragment and structural-center coordinates."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)

BLUE = "#2e6f9e"
TEAL = "#159a9c"
PURPLE = "#7252a5"
INK = "#19324a"
GREY = "#647481"
PALE = "#f4f7f8"


def base(ax, letter, title, subtitle):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.015, 0.015), 0.97, 0.97,
                                boxstyle="round,pad=0.012,rounding_size=0.025",
                                linewidth=0.8, edgecolor="#cad5dc", facecolor=PALE))
    ax.text(0.055, 0.925, letter, fontsize=10, weight="bold", color="white",
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.15", fc=INK, ec="none"))
    ax.text(0.12, 0.94, title, fontsize=10.5, weight="bold", color=INK,
            ha="left", va="center")
    ax.text(0.12, 0.855, subtitle, fontsize=7.6, color=GREY,
            ha="left", va="center")


def nodes(ax, xy, bonds=True):
    xy = np.asarray(xy)
    if bonds:
        for i in range(len(xy)):
            for j in range(i + 1, len(xy)):
                if np.linalg.norm(xy[i] - xy[j]) < 0.48:
                    ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]],
                            color="#82919b", lw=1.5, zorder=1)
    ax.scatter(xy[:, 0], xy[:, 1], s=44, facecolor="#536875",
               edgecolor="white", linewidth=0.7, zorder=3)


def panel_ferrocene(ax):
    base(ax, "a", "Ferrocene", "ring center / Fe distance")
    t = np.linspace(0, 2 * np.pi, 7)[:-1]
    for y in (0.31, 0.61):
        xy = np.c_[0.38 + 0.18 * np.cos(t), y + 0.07 * np.sin(t)]
        ax.plot(*xy[[*range(6), 0]].T, color="#8999a4", lw=1.2)
        ax.scatter(xy[:, 0], xy[:, 1], s=30, color="#536875", zorder=3)
    ax.scatter([0.38], [0.46], s=130, color="#e39a2d", edgecolor="white", zorder=4)
    ax.plot([0.38, 0.38], [0.38, 0.54], ls="--", lw=2, color=TEAL)
    ax.scatter([0.38], [0.46], marker="D", s=26, color=TEAL, zorder=5)
    ax.text(0.68, 0.47, "center\ncoordinate", color=TEAL, fontsize=8,
            ha="center", va="center")


def panel_formic(ax):
    base(ax, "b", "Formic acid--water", "fragment pose / pseudo-contacts")
    xy = np.array([[0.22, 0.47], [0.33, 0.66], [0.37, 0.28],
                   [0.70, 0.63], [0.80, 0.39], [0.63, 0.22]])
    nodes(ax, xy)
    for i, j in ((1, 3), (2, 5), (3, 4)):
        ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]],
                ls="--", lw=1.8, color=TEAL)
    ax.scatter([0.61], [0.46], marker="D", s=38, color=TEAL, zorder=5)
    ax.text(0.73, 0.18, "typed contact rows", color=TEAL, fontsize=7.7,
            ha="center")


def panel_water(ax):
    base(ax, "c", "Water dimer", "six relative rigid-body rows")
    left = np.array([[0.25, 0.46], [0.17, 0.60], [0.33, 0.60]])
    right = left + np.array([0.43, -0.02])
    nodes(ax, left)
    nodes(ax, right)
    ax.plot([0.33, 0.76], [0.60, 0.58], ls="--", lw=1.8, color=TEAL)
    for y, c in ((0.32, BLUE), (0.69, PURPLE)):
        ax.annotate("", xy=(0.80, y), xytext=(0.39, y),
                    arrowprops=dict(arrowstyle="-|>", lw=1.8, color=c))
    ax.text(0.59, 0.18, "translation + orientation", color=BLUE,
            fontsize=7.7, ha="center")


def panel_benzene(ax):
    base(ax, "d", "Benzene--water", "ring center / fragment pose")
    t = np.linspace(0, 2 * np.pi, 7)[:-1]
    ring = np.c_[0.35 + 0.17 * np.cos(t), 0.48 + 0.13 * np.sin(t)]
    ax.plot(*ring[[*range(6), 0]].T, color="#82919b", lw=1.5)
    ax.scatter(ring[:, 0], ring[:, 1], s=38, color="#536875", zorder=3)
    ax.scatter([0.35], [0.48], marker="D", s=36, color=TEAL, zorder=4)
    water = np.array([[0.72, 0.62], [0.65, 0.72], [0.80, 0.72]])
    nodes(ax, water)
    ax.annotate("", xy=(0.64, 0.61), xytext=(0.42, 0.52),
                arrowprops=dict(arrowstyle="-|>", lw=1.8, color=BLUE))
    ax.text(0.56, 0.20, "center + pose", color=TEAL, fontsize=8, ha="center")


def panel_allyl(ax):
    base(ax, "e", r"$\eta^3$ allyl--PdCl", "supplied ligand center")
    metal = np.array([0.50, 0.38])
    allyl = np.array([[0.35, 0.61], [0.50, 0.70], [0.65, 0.61]])
    nodes(ax, allyl)
    ax.scatter(*metal, s=105, color="#9b63b5", edgecolor="white", zorder=4)
    ax.plot([metal[0], 0.50], [metal[1], 0.61], ls="--", lw=1.8, color=TEAL)
    ax.scatter([0.50], [0.61], marker="D", s=34, color=TEAL, zorder=5)
    ax.scatter([0.50], [0.82], s=42, color="#6ca874", edgecolor="white", zorder=4)
    ax.plot([0.50, 0.50], [0.70, 0.82], color="#82919b", lw=1.5)
    ax.text(0.77, 0.42, "$\\eta^3$\ncenter", color=TEAL, fontsize=8,
            ha="center", va="center")


def panel_pose(ax):
    base(ax, "f", "Closed fragment graph", "quaternion / exponential-map pose")
    for x in (0.27, 0.70):
        xy = np.array([[x - 0.10, 0.43], [x, 0.57], [x + 0.10, 0.43],
                       [x, 0.30]])
        nodes(ax, xy)
    ax.plot([0.38, 0.59], [0.48, 0.48], ls="--", lw=1.8, color=TEAL)
    ax.annotate("", xy=(0.61, 0.73), xytext=(0.39, 0.73),
                arrowprops=dict(arrowstyle="-|>", lw=1.8, color=PURPLE))
    ax.text(0.50, 0.18, "relative pose", color=PURPLE, fontsize=8, ha="center")


def main():
    fig, axes = plt.subplots(2, 3, figsize=(8.2, 5.3), dpi=220)
    panels = (panel_ferrocene, panel_formic, panel_water,
              panel_benzene, panel_allyl, panel_pose)
    for ax, fun in zip(axes.flat, panels):
        fun(ax)
    fig.suptitle("Special coordinates remain attached to physical objects",
                 fontsize=14, weight="bold", color=INK, y=0.995)
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.045, top=0.91,
                        wspace=0.055, hspace=0.10)
    fig.savefig(OUT / "special_coordinate_structures.pdf", bbox_inches="tight")
    fig.savefig(OUT / "special_coordinate_structures.png", dpi=300,
                bbox_inches="tight", transparent=False)


if __name__ == "__main__":
    main()
