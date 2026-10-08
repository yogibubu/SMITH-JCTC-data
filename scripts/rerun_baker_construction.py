"""Rebuild the Baker TS ORACLE/SONIC coordinate sets with the installed code.

This is a construction audit, not a geometry-optimization benchmark. It writes
all outputs outside the nano-MATRIX source checkout and records failures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tomllib


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], log_path: Path) -> int:
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    process = subprocess.run(command, text=True, capture_output=True, env=env)
    log_path.write_text(
        "$ " + " ".join(command) + "\n" + process.stdout + process.stderr,
        encoding="utf-8",
    )
    return process.returncode


def field(text: str, name: str) -> str | None:
    match = re.search(rf"^{re.escape(name)} ([^\n]+)$", text, re.MULTILINE)
    return match.group(1).strip() if match else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    manifest_path = repo / "benchmarks/smith/baker_ts/manifest.toml"
    manifest = tomllib.loads(manifest_path.read_text(encoding="utf-8"))
    commit = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    tracked_changes = subprocess.check_output(
        ["git", "-C", str(repo), "status", "--short", "--untracked-files=no"],
        text=True,
    )
    summary: dict = {
        "kind": "ORACLE/SMITH construction audit; no QM optimization",
        "source_commit": commit,
        "source_tracked_changes": tracked_changes.splitlines(),
        "manifest_sha256": sha256(manifest_path),
        "cases": [],
    }
    for case in manifest["cases"]:
        number = case["id"]
        source = manifest_path.parent / case["geometry"]
        case_dir = out / f"{number:02d}_{case['name']}"
        case_dir.mkdir(exist_ok=True)
        xyzin = case_dir / "geometry.xyzin"
        report = case_dir / "oracle.json"
        record = {
            "id": number,
            "name": case["name"],
            "source_geometry": str(source),
            "source_sha256": sha256(source),
        }
        command = [
            "matrix", "oracle", "analyze", str(source), "-o", str(xyzin),
            "--report", str(report), "--freeze-onic-contract",
            "--stationary-point", "transition-state", "--json",
        ]
        record["oracle_exit"] = run(command, case_dir / "oracle.log")
        if record["oracle_exit"] == 0:
            command = [
                "matrix", "smith", "build", str(xyzin), "--oracle-contract",
                "--symmetrize", "--stationary-point", "transition-state",
                "--no-diagnostics",
            ]
            record["smith_exit"] = run(command, case_dir / "smith.log")
            if record["smith_exit"] == 0:
                state = xyzin.read_text(encoding="utf-8")
                for name in (
                    "POINT_GROUP", "TARGET_RANK", "RANK", "GIC_COUNT",
                    "OUT_OF_PLANE_MODE", "SOOP_FINAL_RANK",
                    "SOOP_FINAL_NORMALIZED_CONDITION",
                ):
                    record[name.lower()] = field(state, name)
                record["xyzin_sha256"] = sha256(xyzin)
        summary["cases"].append(record)
        print(f"{number:02d} {case['name']}: ORACLE={record['oracle_exit']} "
              f"SMITH={record.get('smith_exit', 'skipped')} "
              f"rank={record.get('rank', 'n/a')}/{record.get('target_rank', 'n/a')}",
              flush=True)
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
