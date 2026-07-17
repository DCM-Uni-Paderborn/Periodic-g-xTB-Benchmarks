#!/usr/bin/env python3
"""Exercise the verifier's raw-output hash gate without modifying evidence."""

from __future__ import annotations

import verify_matrix as verify


def main() -> None:
    case = next(case for case in verify.expanded_cases() if case["name"] == "ch4_full" and case["ranks"] == 1)
    run_dir = verify.ROOT / "runs" / "ch4_full_p1_partial_replicated"
    verify.checked_metadata(run_dir, "PARTIAL_REPLICATED", case)

    real_sha256 = verify.sha256

    def forged_sha256(path):
        if path.name == "cp2k.out":
            return "0" * 64
        return real_sha256(path)

    verify.sha256 = forged_sha256
    try:
        verify.checked_metadata(run_dir, "PARTIAL_REPLICATED", case)
    except RuntimeError as error:
        if "raw-output hash mismatch" not in str(error):
            raise
        print("PASS: forged cp2k.out digest rejected fail-closed")
    else:
        raise RuntimeError("forged cp2k.out digest was incorrectly accepted")


if __name__ == "__main__":
    main()
