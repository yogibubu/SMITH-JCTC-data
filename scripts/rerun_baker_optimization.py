"""Run current ORACLE/SONIC Baker TS coordinates with LINK and GDV.

The campaign uses a user-authorized 2-thread, 4-GB envelope per job. It does
not use the historical frozen-coordinate benchmark manifest. Each case has
its own exact input-bound calculation-authorization bundle and full log.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tomllib

from matrix_core.calculation_launch import (
    CalculationResources,
    authorize_calculation_launch,
    build_calculation_launch_plan,
    write_calculation_authorization_bundle,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_case(case: dict, *, input_root: Path, out_root: Path,
             executable: Path, repository_commit: str,
             coordinate_kind: str, trust_radius: float,
             hessian_update: str, step_model: str,
             frozen_chart_manifest: Path | None) -> dict:
    number = int(case["id"])
    name = str(case["name"])
    source = input_root / f"{number:02d}_{name}" / "geometry.xyzin"
    if not source.is_file():
        raise FileNotFoundError(source)
    directory = out_root / f"{number:02d}_{name}"
    directory.mkdir(parents=True, exist_ok=False)
    command = [
        "matrix", "link", "optimize", str(source),
        "--run-dir", str(directory),
        "--backend", "gdv",
        "--route", "#p PBE1PBE/def2SVP Force",
        "--charge", str(case["charge"]),
        "--multiplicity", str(case["multiplicity"]),
        "--coordinate-kind", coordinate_kind,
        "--stationary-point", "transition-state",
        "--max-steps", "50",
        "--trust-radius", str(trust_radius),
        "--max-trust-radius", str(trust_radius),
        "--hessian-update", hessian_update,
        "--step-model", step_model,
        "--allow-approximate-final-hessian", "--skip-irc",
        "--initial-hessian-seed-backend", "gdv",
        "--initial-hessian-seed-route", "#p PBE1PBE/def2SVP Freq",
        "--executable", str(executable),
        "--backend-workers", "2", "--backend-memory-gb", "4",
        "--timeout", "3600",
    ]
    if frozen_chart_manifest is not None:
        command.extend((
            "--frozen-chart-manifest", str(frozen_chart_manifest),
            "--frozen-chart-case", str(number),
        ))
    plan = build_calculation_launch_plan(
        backend="gdv", input_path=source, command=command,
        resources=CalculationResources(1, 2, 4.0, 1), workdir=directory,
    )
    approval = authorize_calculation_launch(
        plan, approved_plan_sha256=plan.sha256,
        authorized_by="Vincenzo Barone, approval in ct-2026-01770x revision chat",
    )
    bundle = write_calculation_authorization_bundle(
        directory / "launch_authorization.json", plan, approval
    )
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    environment = dict(os.environ)
    environment["MATRIX_CALCULATION_AUTHORIZATION_BUNDLE"] = str(bundle)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    with (directory / "launch.log").open("w") as stream:
        result = subprocess.run(
            command, stdout=stream, stderr=subprocess.STDOUT,
            env=environment, timeout=3700, check=False,
        )
    summary_path = directory / "optimizer_summary.json"
    record = {
        "case_id": number, "name": name, "coordinate_kind": coordinate_kind,
        "trust_radius_angstrom": trust_radius,
        "hessian_update": hessian_update,
        "step_model": step_model,
        "source_commit": repository_commit,
        "source_xyzin_sha256": sha256(source),
        "launch_plan_sha256": plan.sha256,
        "returncode": result.returncode,
        "launch_log_sha256": sha256(directory / "launch.log"),
    }
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text())
        record.update({
            "summary_sha256": sha256(summary_path),
            "status": summary.get("status"),
            "converged": summary.get("converged"),
            "optimization_steps": summary.get("optimization_steps"),
            "qm_evaluations": summary.get("qm_evaluations"),
            "final_energy_hartree": summary.get("final_energy_hartree"),
            "approximate_hessian_index": summary.get("final_hessian_index"),
        })
    (directory / "campaign_record.json").write_text(
        json.dumps(record, indent=2) + "\n"
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cases", default="all")
    parser.add_argument("--coordinate-kind", choices=("sonic", "sycart", "cartesian", "simple-redundant"),
                        default="sonic")
    parser.add_argument("--trust-radius", type=float, default=0.3)
    parser.add_argument("--hessian-update", choices=("auto", "bofill", "sella_ts_bfgs", "sr1"), default="auto")
    parser.add_argument("--step-model", choices=("dense-rfo", "quadraticized-rfo"), default="dense-rfo")
    parser.add_argument("--frozen-reference", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        raise ValueError("approved campaign permits at most four concurrent jobs")
    if not 0.0 < args.trust_radius <= 0.3:
        raise ValueError("trial trust radius must be between zero and the approved 0.3 A ceiling")
    ids = set(range(1, 26)) if args.cases == "all" else {
        int(item) for item in args.cases.split(",")
    }
    if not ids or not ids.issubset(set(range(1, 26))):
        raise ValueError("case IDs must be in 1..25")
    repo = args.repo.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = repo / "benchmarks/smith/baker_ts/manifest.toml"
    cases = [case for case in tomllib.loads(manifest_path.read_text())["cases"]
             if int(case["id"]) in ids]
    commit = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    results: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                run_case, case, input_root=args.input_root.resolve(),
                out_root=out, executable=args.executable.resolve(),
                repository_commit=commit, coordinate_kind=args.coordinate_kind,
                trust_radius=args.trust_radius,
                hessian_update=args.hessian_update,
                step_model=args.step_model,
                frozen_chart_manifest=(
                    repo / "benchmarks/smith/baker_ts/gdv_sonic_reference/manifest.tsv"
                    if args.frozen_reference else None
                ),
            ): case for case in cases
        }
        for future in as_completed(futures):
            case = futures[future]
            try:
                record = future.result()
            except Exception as error:
                record = {"case_id": case["id"], "name": case["name"],
                          "returncode": -1, "error": str(error)}
            results.append(record)
            print(f"{record['case_id']:02d} {record['name']}: "
                  f"{record.get('status', record.get('error', 'failed'))}, "
                  f"steps={record.get('optimization_steps', 'n/a')}", flush=True)
            report = {
                "schema": "ct-2026-01770x.baker_updated.v1",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "repository_commit": commit,
                "coordinate_kind": args.coordinate_kind,
                "trust_radius_angstrom": args.trust_radius,
                "hessian_update": args.hessian_update,
                "step_model": args.step_model,
                "frozen_reference": args.frozen_reference,
                "manifest_sha256": sha256(manifest_path),
                "workers": args.workers,
                "cases": sorted(results, key=lambda item: item["case_id"]),
            }
            (out / "campaign_summary.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )


if __name__ == "__main__":
    main()
