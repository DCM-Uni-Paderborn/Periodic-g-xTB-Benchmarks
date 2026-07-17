#!/usr/bin/env python3
"""Wait for selected EOS points, fit the local well, and run own-minimum SPs.

The helper is intentionally campaign-bound and writes a resumable JSON status.
It is used only for local refinement grids whose points were launched by a
separate, stamped production runner.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


BUILD = Path("/home/kuehne88/work/codex-gxtb-post5582-clean-20260714")
MANIFEST = BUILD / "campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--solid", action="append", required=True)
    parser.add_argument("--scale", action="append", type=float, required=True)
    parser.add_argument("--status-name", required=True)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--timeout-hours", type=float, default=36.0)
    parser.add_argument("--wait-for-existing-equilibrium", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    scripts = root / "Goldzak12" / "scripts"
    sys.path.insert(0, str(scripts))

    import run_goldzak12_benchmark as base
    import run_goldzak12_eos_benchmark as eos
    import run_goldzak12_k_convergence as kconv

    status_path = root / "Goldzak12" / "data" / args.status_name
    solids = tuple(dict.fromkeys(args.solid))
    scales = tuple(sorted(set(args.scale)))
    refs = {ref.solid: ref for ref in base.LC10_PAPER_REFERENCES}
    unknown = sorted(set(solids) - set(refs))
    if unknown:
        raise ValueError(f"unknown LC10 solids: {unknown}")

    campaign, paths = base.validated_gxtb_campaign_from_manifest(
        MANIFEST,
        BUILD / "sources/cp2k",
        BUILD / "sources/save_tblite",
    )

    def write_status(phase: str, **extra: object) -> None:
        status_path.parent.mkdir(parents=True, exist_ok=True)
        status_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "phase": phase,
                    "mesh": args.mesh,
                    "solids": list(solids),
                    "fit_scales": list(scales),
                    "campaign_fingerprint_sha256": campaign["fingerprint_sha256"],
                    "updated_at_utc": utc_now(),
                    **extra,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

    def eos_output(solid: str, scale: float) -> Path:
        project = eos.eos_project(solid, "GXTB", args.mesh, scale)
        return (
            base.ROOT
            / "runs"
            / "eos"
            / "GXTB"
            / solid
            / args.mesh
            / eos.scale_tag(scale, "GXTB")
            / f"{project}.out"
        )

    def active_equilibrium_processes() -> list[str]:
        if not args.wait_for_existing_equilibrium:
            return []
        active: list[str] = []
        for solid in solids:
            pattern = f"{solid}_GXTB_independent_eos_minimum_{args.mesh}.inp"
            result = subprocess.run(
                ["pgrep", "-af", pattern],
                check=False,
                capture_output=True,
                text=True,
            )
            active.extend(line for line in result.stdout.splitlines() if line.strip())
        return active

    deadline = time.monotonic() + args.timeout_hours * 3600.0
    while True:
        pending: list[str] = []
        issues: list[str] = []
        for solid in solids:
            for scale in scales:
                output = eos_output(solid, scale)
                label = f"{solid}/{args.mesh}/{scale:.5f}"
                if not base.output_ok(output):
                    pending.append(label)
                    continue
                issue = base.completed_stamp_campaign_issue(
                    output,
                    campaign,
                    executable_role="cp2k",
                    require_completed=True,
                )
                if issue:
                    pending.append(label)
                    issues.append(f"{label}: {issue}")
        active = active_equilibrium_processes()
        if not pending and not active:
            break
        write_status(
            "waiting_for_eos",
            pending_points=pending,
            campaign_issues=issues,
            active_existing_equilibrium_processes=active,
        )
        if time.monotonic() >= deadline:
            write_status(
                "timeout",
                pending_points=pending,
                campaign_issues=issues,
                active_existing_equilibrium_processes=active,
            )
            raise TimeoutError("timed out waiting for selected EOS points")
        time.sleep(max(5, args.poll_seconds))

    fits: list[dict[str, object]] = []
    for solid in solids:
        raw = eos.load_eos_points(refs[solid], "GXTB", args.mesh, scales, campaign)
        points = [
            (a, scale, float(energy), ok)
            for a, scale, energy, ok in raw
            if energy is not None and ok
        ]
        fit = eos.fit_gxtb_eos(points)
        if fit.get("fit_status") != "quadratic" or not fit.get("a_eos_A"):
            write_status("fit_failed", solid=solid, fit=fit, points=points)
            raise RuntimeError(f"non-reportable local EOS fit for {solid}/{args.mesh}: {fit}")
        fits.append(
            {
                "solid": solid,
                "method": "GXTB",
                "eos_mesh": args.mesh,
                **fit,
                "n_requested": len(scales),
                "n_completed": len(points),
                "n_converged_raw": len(points),
                "n_explicit_excluded": 0,
                "n_unresolved_branch_candidates": 0,
                "fit_scope": "local_refinement_bracket",
                "fit_scales": [f"{value:.5f}" for value in scales],
            }
        )

    kconv.PAPER_SYSTEMS = solids
    specs = kconv.equilibrium_specs(base.ROOT, fits, ("GXTB",))
    write_status(
        "running_equilibrium",
        fits=fits,
        outputs=[str(spec[2]) for spec in specs],
    )
    eos.run_jobs(
        specs,
        paths["cp2k"],
        jobs=len(specs),
        threads=1,
        force=True,
        retry_scf=False,
        campaign_fingerprint=campaign,
        execution_pool=None,
        campaign_bind_all_methods=True,
    )
    kconv.validate_campaign_outputs(specs, campaign)
    records = kconv.collect_independent_values(base.ROOT, fits, campaign, ("GXTB",))
    write_status(
        "complete",
        fits=fits,
        records=records,
        outputs=[str(spec[2]) for spec in specs],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
