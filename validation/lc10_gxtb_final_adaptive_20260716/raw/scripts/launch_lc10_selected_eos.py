#!/usr/bin/env python3
"""Launch a campaign-bound selected LC10 EOS mesh on Terok."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


BUILD = Path("/home/kuehne88/work/codex-gxtb-post5582-clean-20260714")
MANIFEST = BUILD / "campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--solid", action="append", required=True)
    parser.add_argument("--scale", action="append", type=float, required=True)
    parser.add_argument("--status-name", required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    sys.path.insert(0, str(root / "Goldzak12" / "scripts"))

    import run_goldzak12_benchmark as base
    import run_goldzak12_eos_benchmark as eos
    import run_goldzak12_k_convergence as kconv

    solids = tuple(dict.fromkeys(args.solid))
    scales = tuple(sorted(set(args.scale)))
    status = root / "Goldzak12" / "data" / args.status_name
    campaign, paths = base.validated_gxtb_campaign_from_manifest(
        MANIFEST,
        BUILD / "sources/cp2k",
        BUILD / "sources/save_tblite",
    )
    specs = eos.eos_job_specs(args.mesh, scales, ("GXTB",), solids)

    def write_status(phase: str) -> None:
        status.parent.mkdir(parents=True, exist_ok=True)
        status.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "phase": phase,
                    "mesh": args.mesh,
                    "solids": list(solids),
                    "scales": list(scales),
                    "outputs": [str(spec[2]) for spec in specs],
                    "campaign_fingerprint_sha256": campaign["fingerprint_sha256"],
                    "updated_at_utc": datetime.now(timezone.utc).isoformat(),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

    write_status("running")
    eos.run_jobs(
        specs,
        paths["cp2k"],
        jobs=len(specs),
        threads=1,
        force=False,
        retry_scf=False,
        campaign_fingerprint=campaign,
        execution_pool=None,
        campaign_bind_all_methods=True,
    )
    kconv.validate_campaign_outputs(specs, campaign)
    write_status("complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
