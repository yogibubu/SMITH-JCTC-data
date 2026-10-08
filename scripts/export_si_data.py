#!/usr/bin/env python3
"""Export SI tables and Baker starting structures into machine-readable files."""

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SI_SOURCE = ROOT / "source/Supporting_Information_clean.tex"
MAIN_SOURCE = ROOT / "source/Manuscript_clean.tex"


def plain(cell: str) -> str:
    cell = re.sub(r"\\(?:toprule|midrule|bottomrule|endfirsthead|endhead|endfoot|endlastfoot)", "", cell)
    cell = re.sub(r"\\(?:texttt|mathrm|mathbf|emph|textbf)\{([^{}]*)\}", r"\1", cell)
    cell = cell.replace(r"\Bmat", "B").replace(r"\ESALCI", "ESALC-I")
    cell = cell.replace(r"\rightarrow", " -> ").replace(r"\#", "#")
    cell = cell.replace(r"\times", "x").replace(r"\pm", "+/-")
    cell = cell.replace(r"\leq", "<=").replace(r"\geq", ">=")
    cell = cell.replace(r"\%", "%").replace("~", " ")
    cell = cell.replace(r"\(", "").replace(r"\)", "")
    cell = re.sub(r"\\(?:mathsf|mathrm|text)\{([^{}]*)\}", r"\1", cell)
    cell = re.sub(r"\\([A-Za-z]+)", r"\1", cell)
    cell = cell.replace(r"\ ", " ")
    cell = cell.replace("{", "").replace("}", "")
    cell = cell.replace("$", "")
    cell = re.sub(r"\s+", " ", cell)
    return cell.strip()


def export_tables(source: str, folder: str) -> None:
    output = ROOT / "data" / folder
    output.mkdir(parents=True, exist_ok=True)
    total = 0
    blocks = [m.group() for m in re.finditer(r"\\begin\{table\}.*?\\end\{table\}", source, re.S)]
    blocks += [m.group() for m in re.finditer(r"\\begin\{longtable\}.*?\\end\{longtable\}", source, re.S)]
    for table in blocks:
        label_match = re.search(r"\\label\{([^}]+)\}", table)
        tabular_match = re.search(r"\\begin\{(tabularx|tabular|longtable)\}(.*?)\\end\{\1\}", table, re.S)
        if not label_match or not tabular_match:
            continue
        label = label_match.group(1).replace("tab:", "")
        content = tabular_match.group(2)
        if r"\toprule" in content:
            content = content.split(r"\toprule", 1)[1]
        rows = []
        for part in re.split(r"(?<!\\)\\\\", content):
            if "&" not in part:
                continue
            cells = [plain(c) for c in re.split(r"(?<!\\)&", part)]
            if any(cells):
                rows.append(cells)
        if tabular_match.group(1) == "longtable":
            header = rows[0]
            rows = [header] + [row for row in rows[1:] if row != header]
        if not rows:
            continue
        path = output / f"{label}.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerows(rows)
        total += 1
    print(f"Exported {total} tables to {folder}")


def export_baker(source: str) -> None:
    output = ROOT / "data/baker_ts_starting_geometries"
    output.mkdir(parents=True, exist_ok=True)
    pattern = re.compile(
        r"\\subsection\*\{Baker TS\s+(\d+):\s+\\texttt\{([^}]+)\}\}"
        r"\s*Charge/multiplicity:\s*\\texttt\{([^}]+)\}\."
        r"\s*\\begin\{verbatim\}(.*?)\\end\{verbatim\}",
        re.S,
    )
    index = []
    for match in pattern.finditer(source):
        number, name, charge_mult, body = match.groups()
        atoms = [line.strip() for line in body.strip().splitlines() if line.strip()]
        charge, multiplicity = charge_mult.split()
        filename = f"{int(number):02d}_{name}.xyz"
        comment = f"Baker TS {int(number):02d}; charge={charge}; multiplicity={multiplicity}; angstrom"
        (output / filename).write_text(f"{len(atoms)}\n{comment}\n" + "\n".join(atoms) + "\n")
        index.append((int(number), name, filename, charge, multiplicity, len(atoms)))
    assert len(index) == 25, f"Expected 25 Baker structures, got {len(index)}"
    with (output / "index.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("case", "name", "xyz_file", "charge", "multiplicity", "atoms"))
        writer.writerows(index)
    print("Exported 25 Baker starting structures")


if __name__ == "__main__":
    si = SI_SOURCE.read_text()
    export_tables(si, "si_tables")
    export_tables(MAIN_SOURCE.read_text(), "main_tables")
    export_baker(si)
