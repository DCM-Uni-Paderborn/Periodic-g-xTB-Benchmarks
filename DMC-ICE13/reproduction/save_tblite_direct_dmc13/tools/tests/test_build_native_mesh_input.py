#!/usr/bin/env python3
"""Regression tests for the canonical native BvK mesh rewriter."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOL = Path(__file__).resolve().parents[1] / "build_native_mesh_input.py"


def input_text(mesh: int, shift: str) -> str:
    return (
        "&FORCE_EVAL\n"
        "  &DFT\n"
        "    &KPOINTS\n"
        f"      SCHEME MACDONALD {mesh} {mesh} {mesh} {shift} {shift} {shift}\n"
        "    &END KPOINTS\n"
        "  &END DFT\n"
        "&END FORCE_EVAL\n"
    )


class NativeMeshInputTests(unittest.TestCase):
    def rewrite(self, source_mesh: int, source_shift: str, target_mesh: int) -> tuple[str, dict[str, object]]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.inp"
            target = root / "target.inp"
            provenance = root / "rewrite.json"
            source.write_text(input_text(source_mesh, source_shift), encoding="utf-8")
            subprocess.run(
                [
                    sys.executable,
                    str(TOOL),
                    str(source),
                    str(target),
                    str(target_mesh),
                    "--provenance",
                    str(provenance),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            return target.read_text(encoding="utf-8"), json.loads(provenance.read_text(encoding="utf-8"))

    def test_odd_to_even_recomputes_shift(self) -> None:
        target, provenance = self.rewrite(3, "0.0", 4)
        self.assertIn("SCHEME MACDONALD 4 4 4 0.375 0.375 0.375", target)
        self.assertEqual(provenance["source_shift"], [0.0, 0.0, 0.0])
        self.assertEqual(provenance["target_shift"], [0.375, 0.375, 0.375])

    def test_even_to_odd_clears_shift(self) -> None:
        target, provenance = self.rewrite(8, "0.4375", 9)
        self.assertIn("SCHEME MACDONALD 9 9 9 0.0 0.0 0.0", target)
        self.assertEqual(provenance["target_shift"], [0.0, 0.0, 0.0])

    def test_six_mesh_uses_exact_canonical_float(self) -> None:
        target, _ = self.rewrite(5, "0.0", 6)
        self.assertIn(
            "SCHEME MACDONALD 6 6 6 0.4166666666666667 0.4166666666666667 0.4166666666666667",
            target,
        )

    def test_only_scheme_line_changes(self) -> None:
        source = input_text(4, "0.375")
        target, provenance = self.rewrite(4, "0.375", 5)
        changed = [
            index
            for index, (left, right) in enumerate(zip(source.splitlines(), target.splitlines()), start=1)
            if left != right
        ]
        self.assertEqual(changed, [4])
        self.assertEqual(provenance["changed_line"], 4)


if __name__ == "__main__":
    unittest.main()
