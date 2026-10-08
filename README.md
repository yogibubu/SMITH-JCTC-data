# Data and scripts for ct-2026-01770x.R1

This directory contains machine-readable numerical tables, starting structures,
representative input and output files, and analysis scripts for **Local,
Non-Redundant Molecular Coordinates from Disjoint Primitive Families:
Symmetry-Adapted Construction and Applications**. It is separate from the
Supporting Information PDF.

## Software

- Coordinate construction: [SMITH-standalone, fixed public revision
  358d265](https://github.com/yogibubu/SMITH-standalone/tree/358d26537535c6f2510231b51beccb595d3dc829).
- The external optimization examples also use Gaussian, geomeTRIC, PySCF,
  GFN2-xTB, and LINK, as specified for each result in the article and Supporting Information.
- The analysis scripts in `scripts/` are supplied for inspection and reuse.
  Scripts importing the MATRIX packages require `MATRIX_ROOT` to point to a
  compatible source checkout; Gaussian and other external programs require
  their own licenses or installations.

## Data

| Location | Contents |
| --- | --- |
| `data/main_tables/` | Seven main-text tables as CSV, including optimization steps and transition-state outcomes. |
| `data/si_tables/` | All 19 SI tables as CSV, including ranks, derivatives, Hessian coupling, construction times, and optimization summaries. |
| `data/baker_ts_starting_geometries/` | The 25 Cartesian Baker starting structures in XYZ format, with charge and multiplicity in `index.csv`. |
| `data/hessian_cyclobutane/` | Stationary HF/STO-3G formatted checkpoint and derived delocalized-coordinate data. |
| `data/hessian_norbornane/` | Stationary geometry, Cartesian Hessian array, primitive definitions, and analysis results. |
| `data/gaussian_matched_coordinate_inputs/` | Representative paired Gaussian inputs for covalent systems. |
| `data/noncovalent_gfn2xtb_20260824/` | Starting structures and representative GFN2-xTB preprocessing outputs. |
| `data/updated_benchmarks/goppeg_all_sycart_32_steps/` | GOPPEG starting structure, complete ESALC-C input, trajectory, and converged numerical summary. |
| `data/application_results.json` | Cyclobutane Hessian and puckering data used in the application figures. |

`source/` contains the manuscript and SI TeX from which the CSV tables and
Baker XYZ files were exported. Run `python3 scripts/export_si_data.py` to
regenerate them and `python3 scripts/check_release.py` for basic integrity
checks. The CSV exports preserve the SI table values; symbolic quantities
remain text fields. No LINK transition-state trial data are included.

The Gaussian paired inputs are provided as representative files. The
published SI tables and the machine-readable CSV exports record the
numerical results used here.
