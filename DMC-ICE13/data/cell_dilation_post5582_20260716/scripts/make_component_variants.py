#!/usr/bin/env python3
"""Create schema-supported g-xTB component-deletion diagnostics.

The variants retain the complete exported g-xTB parameter file except for the
single named change.  They are diagnostic models, not reparameterizations.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def remove_top_level_block(text: str, start: str, stop: str) -> str:
    pattern = re.compile(
        rf"(?ms)^\[{re.escape(start)}\]\n.*?(?=^\[{re.escape(stop)}\]\n)"
    )
    result, count = pattern.subn("", text)
    if count != 1:
        raise RuntimeError(f"expected one [{start}]...before [{stop}] block, got {count}")
    return result


def zero_environment_coefficients(text: str) -> tuple[str, int, int]:
    pattern = re.compile(r"(?m)^coeffs_env\s*=\s*\[([^\n]*)\]$")
    total_values = 0

    def replace(match: re.Match[str]) -> str:
        nonlocal total_values
        values = [item.strip() for item in match.group(1).split(",") if item.strip()]
        if not values:
            raise RuntimeError("empty coeffs_env array")
        total_values += len(values)
        zeros = ", ".join("0.0000000000000000" for _ in values)
        return f"coeffs_env = [ {zeros} ]"

    result, arrays = pattern.subn(replace, text)
    if arrays == 0:
        raise RuntimeError("no coeffs_env arrays found")
    return result, arrays, total_values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    source = args.source.read_text()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    no_exchange = remove_top_level_block(source, "exchange", "spin")
    no_multipole = remove_top_level_block(source, "multipole", "exchange")
    no_acp = remove_top_level_block(source, "acp", "element")
    no_exchange_no_acp = remove_top_level_block(no_exchange, "acp", "element")
    frozen_qvszp, arrays, values = zero_environment_coefficients(source)

    variants = {
        "gxtb_full.toml": source,
        "gxtb_no_exchange.toml": no_exchange,
        "gxtb_frozen_qvszp.toml": frozen_qvszp,
        "gxtb_no_anisotropic_multipole.toml": no_multipole,
        "gxtb_no_acp.toml": no_acp,
        "gxtb_no_exchange_no_acp.toml": no_exchange_no_acp,
    }
    for name, content in variants.items():
        if not content.endswith("\n"):
            content += "\n"
        (args.output_dir / name).write_text(content)

    provenance = {
        "source": str(args.source),
        "variants": {
            "gxtb_no_exchange.toml": {
                "change": "remove complete top-level [exchange] parameter block",
                "interpretation": "exchange-free schema-supported calculator diagnostic",
            },
            "gxtb_frozen_qvszp.toml": {
                "change": "set every coeffs_env array entry to zero without changing array lengths",
                "arrays_changed": arrays,
                "values_changed": values,
                "interpretation": "freeze q-vSZP coefficient response to effective charge/environment",
            },
            "gxtb_no_anisotropic_multipole.toml": {
                "change": "remove complete top-level [multipole] parameter block",
                "interpretation": "anisotropic multipole-free diagnostic; NOT a switch for all periodic electrostatics or images",
            },
            "gxtb_no_acp.toml": {
                "change": "remove top-level [acp] model marker while retaining the otherwise unchanged parameter export",
                "interpretation": "ACP-free schema-supported calculator diagnostic",
            },
            "gxtb_no_exchange_no_acp.toml": {
                "change": "remove both top-level [exchange] and [acp] model blocks",
                "interpretation": "two-component interaction diagnostic; not a separately parameterized physical model",
            },
        },
    }
    (args.output_dir / "variant_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
