from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "DMC-ICE13" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import make_dmc_ice13_plots as plots  # noqa: E402


def profile(offset: float = 0.0) -> dict[str, float]:
    return {phase: offset + index / 10 for index, phase in enumerate(plots.PHASES)}


def complete_result(offset: float = 0.0) -> dict[str, object]:
    return {"complete": True, "relative_kjmol": profile(offset)}


class DMCGXTBPaperPlotTests(unittest.TestCase):
    def test_k333_loader_includes_only_complete_gxtb_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            old_data = plots.DATA
            plots.DATA = Path(tmp)
            try:
                (plots.DATA / "kpoint_results.json").write_text(
                    json.dumps(
                        {
                            "results": {
                                "k333": {
                                    "GFN1": complete_result(1.0),
                                    "GFN2": complete_result(2.0),
                                    "GXTB": complete_result(3.0),
                                },
                                "gamma": {
                                    "GFN1": complete_result(1.0),
                                    "GFN2": complete_result(2.0),
                                    "GXTB": {
                                        "complete": False,
                                        "relative_kjmol": profile(3.0),
                                    },
                                },
                            }
                        }
                    )
                )
                primary = plots.load_primary_gfn_relative_energies({})
                gamma = plots.load_gamma_gfn_relative_energies({})
            finally:
                plots.DATA = old_data

        self.assertEqual(list(primary), ["GFN1-xTB", "GFN2-xTB", "g-xTB"])
        self.assertEqual(list(gamma), ["GFN1-xTB (Gamma)", "GFN2-xTB (Gamma)"])

    def test_missing_phase_rejects_nominally_complete_gxtb(self) -> None:
        incomplete = complete_result()
        del incomplete["relative_kjmol"]["XVII"]  # type: ignore[index]
        profiles = plots.load_complete_gfn_profiles({"GXTB": incomplete})
        self.assertEqual(profiles, {})

    def test_summary_rows_report_twelve_nonreference_points(self) -> None:
        dmc = profile()
        rows = plots.build_summary_rows({"g-xTB": profile(1.0)}, dmc)
        self.assertEqual(rows[0]["method"], "g-xTB")
        self.assertEqual(rows[0]["N"], 12)
        self.assertEqual(rows[0]["MAE"], "1.0000")

    def test_profile_inputs_and_plot_clause_add_gxtb_dynamically(self) -> None:
        primary = {
            "GFN1-xTB": profile(1.0),
            "GFN2-xTB": profile(2.0),
            "g-xTB": profile(3.0),
        }
        gamma = {
            "GFN1-xTB (Gamma)": profile(1.1),
            "GFN2-xTB (Gamma)": profile(2.1),
            "g-xTB (Gamma)": profile(3.1),
        }
        dft = {method: profile(4.0 + index) for index, method in enumerate(plots.PARENT_DFT_METHODS)}
        with tempfile.TemporaryDirectory() as tmp:
            rel_path = Path(tmp) / "relative.dat"
            err_path = Path(tmp) / "errors.dat"
            rel_columns, err_columns = plots.write_profile_plot_data(
                rel_path,
                err_path,
                profile(),
                primary,
                gamma,
                {**primary, **dft},
            )
            rel_header = rel_path.read_text().splitlines()[0]
            err_header = err_path.read_text().splitlines()[0]

        self.assertIn("g_xTB", rel_header)
        self.assertIn("g_xTB_Gamma_error", err_header)
        self.assertEqual(rel_columns["g-xTB"], 7)
        clause = plots.gnuplot_profile_clause(
            Path("relative.dat"),
            rel_columns,
            ["g-xTB", "g-xTB (Gamma)", "PBE-D4"],
            include_dmc=True,
        )
        self.assertIn("title 'g-xTB'", clause)
        self.assertIn("title 'g-xTB (Γ-point)'", clause)
        self.assertIn(f"using 1:{err_columns['g-xTB']}", plots.gnuplot_profile_clause(
            Path("errors.dat"), err_columns, ["g-xTB"], include_dmc=False
        ))


if __name__ == "__main__":
    unittest.main()
