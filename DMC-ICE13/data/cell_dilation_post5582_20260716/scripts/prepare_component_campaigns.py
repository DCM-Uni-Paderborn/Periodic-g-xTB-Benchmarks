#!/usr/bin/env python3
"""Clone the controlled DMC dilation inputs for component diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


VARIANTS = {
    "no_exchange": "gxtb_no_exchange.toml",
    "frozen_qvszp": "gxtb_frozen_qvszp.toml",
    "no_anisotropic_multipole": "gxtb_no_anisotropic_multipole.toml",
    "no_acp": "gxtb_no_acp.toml",
    "no_exchange_no_acp": "gxtb_no_exchange_no_acp.toml",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_inputs", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("remote_parameter_dir")
    args = parser.parse_args()

    records = []
    for variant, parameter_name in VARIANTS.items():
        target_dir = args.output_root / variant / "inputs"
        target_dir.mkdir(parents=True, exist_ok=True)
        remote_param = f"{args.remote_parameter_dir}/{parameter_name}"
        for source in sorted(args.source_inputs.glob("*.inp")):
            text = source.read_text()
            project_match = re.search(r"(?m)^\s*PROJECT\s+(\S+)\s*$", text)
            if not project_match:
                raise RuntimeError(f"PROJECT missing in {source}")
            project = project_match.group(1)
            text, project_count = re.subn(
                r"(?m)^(\s*PROJECT\s+)\S+(\s*)$",
                rf"\g<1>{project}_{variant}\g<2>",
                text,
                count=1,
            )
            text, param_count = re.subn(
                r"(?m)^(\s*METHOD\s+GXTB\s*)$",
                rf"\g<1>\n          PARAM {remote_param}",
                text,
                count=1,
            )
            if project_count != 1 or param_count != 1:
                raise RuntimeError(
                    f"controlled edit failed for {source}: project={project_count}, param={param_count}"
                )
            target = target_dir / source.name
            target.write_text(text)
            records.append(
                {
                    "variant": variant,
                    "case": source.stem,
                    "source": str(source),
                    "source_sha256": sha256(source),
                    "generated": str(target),
                    "generated_sha256": sha256(target),
                    "remote_parameter": remote_param,
                }
            )

    provenance = {
        "description": "DMC-ICE13 cell-dilation component diagnostics cloned from the post-#5582 full-model controlled inputs",
        "controlled_change": "only PARAM is added and PROJECT is suffixed; geometry, cell, PBC mask, k mesh, SCC and all numerical settings are byte-for-byte inherited otherwise",
        "variants": VARIANTS,
        "cases": records,
    }
    (args.output_root / "campaign_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
