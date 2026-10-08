#!/usr/bin/env python3
"""Check essential data files and a few reported numerical values."""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def rows(path: Path) -> list[list[str]]:
    with path.open(newline="") as handle:
        return list(csv.reader(handle))


def main() -> None:
    baker = rows(ROOT / "data/baker_ts_starting_geometries/index.csv")
    assert len(baker) == 26
    for number, name, filename, charge, multiplicity, atoms in baker[1:]:
        xyz = (ROOT / "data/baker_ts_starting_geometries" / filename).read_text().splitlines()
        assert int(xyz[0]) == int(atoms) == len(xyz) - 2
        assert f"charge={charge}" in xyz[1]
        assert f"multiplicity={multiplicity}" in xyz[1]
    assert len(list((ROOT / "data/main_tables").glob("*.csv"))) == 7
    assert len(list((ROOT / "data/si_tables").glob("*.csv"))) == 19
    for table in (*((ROOT / "data/main_tables").glob("*.csv")), *((ROOT / "data/si_tables").glob("*.csv"))):
        data = rows(table)
        assert len(data) >= 2, table
        assert len({len(row) for row in data}) == 1, table
    norbornane = json.loads((ROOT / "data/hessian_norbornane/analysis.json").read_text())
    assert norbornane["sonic_rank"] == 51
    goppeg = json.loads((ROOT / "data/updated_benchmarks/goppeg_all_sycart_32_steps/result.json").read_text())
    assert goppeg["converged"] is True
    assert goppeg["target_rank"] == 237
    assert goppeg["optimization_steps"] == 32
    assert goppeg["energy_gradient_evaluations"] == 47
    assert abs(goppeg["maximum_cartesian_gradient_hartree_per_bohr"] - 1.8498e-4) < 1e-8
    print("25 Baker structures, 26 tables, norbornane rank, and GOPPEG result verified")


if __name__ == "__main__":
    main()
