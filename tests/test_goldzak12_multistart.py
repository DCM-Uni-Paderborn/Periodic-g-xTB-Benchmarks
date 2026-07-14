from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "Goldzak12" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import classify_gxtb_multistart_branches as classify  # noqa: E402
import run_goldzak12_benchmark as base  # noqa: E402
import run_gxtb_multistart_branches as runner  # noqa: E402


def fake_campaign(cp2k_sha256: str = "cp2k-fixture") -> dict[str, object]:
    dependencies = {"tblite": "d" * 40}
    return base.make_campaign_identity(
        campaign_id="multistart-fixture",
        cp2k_executable_sha256=cp2k_sha256,
        cp2k_loaded_library_sha256="libcp2k-fixture",
        cp2k_cmake_cache_sha256="cp2k-cache-fixture",
        cp2k_embedded_source_revision="c92cc08b4",
        cp2k_source_revision="c" * 40,
        save_tblite_executable_sha256="save-cli-fixture",
        save_tblite_source_revision="s" * 40,
        save_tblite_library_sha256="save-library-fixture",
        save_tblite_cmake_cache_sha256="save-cache-fixture",
        dependency_lock_sha256=base.campaign_fingerprint_sha256(dependencies),
    )


def synthetic_output(solid: str, energy: float, *, continuation: bool) -> str:
    if solid == "LiH":
        electron_count = 8
        atoms = [("Li", 0.8, 0.2)] * 4 + [("H", 1.2, -0.2)] * 4
    else:
        electron_count = 32
        atoms = [("Mg", 1.0, 1.0)] * 4 + [("O", 7.0, -1.0)] * 4
    lines = []
    if continuation:
        lines.append(" WFN_RESTART| Reading restart file")
    lines += [
        f" Number of electrons: {electron_count}",
        f" Number of occupied orbitals: {electron_count // 2}",
        f" Number of molecular orbitals: {electron_count}",
        " Fermi energy: -0.2000000000",
        " MO| 1 -0.500000 0.000000 2.000000",
        f" MO| Total occupied: {electron_count // 2}",
        " MULLIKEN POPULATION ANALYSIS",
    ]
    for index, (element, population, charge) in enumerate(atoms, start=1):
        lines.append(f" {index:3d} {element:2s} 1 {population: .10f} {charge: .10f}")
    lines += [
        "     8 GXTB-FDIIS  0.10E+01    0.1     5.0000E-10     "
        f"{energy: .12f}  0.00E+00",
        " *** SCF run converged in     8 steps ***",
        f" Total energy: {energy: .12f}",
        " Electronic entropic energy: 0.000000000000",
        f" Total energy (extrapolated to T->0): {energy: .12f}",
        f" ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] {energy: .12f}",
        " PROGRAM ENDED",
    ]
    return "\n".join(lines) + "\n"


def fixture_candidate(
    root: Path,
    *,
    ref: base.Reference,
    scale: float,
    mode: str,
    campaign: dict[str, object],
    plan: dict[str, object],
    plan_sha256: str,
    parent_manifest: Path | None,
    parent_restart: Path | None,
) -> tuple[Path, Path]:
    inp, out, restart, manifest = runner.candidate_paths(root, ref.solid, mode, scale)
    text = runner.multistart_input(ref, scale, inp.stem, parent_restart)
    base.write_file(inp, text)
    mode_offset = {"cold": 0.0, "ascending": 1.0e-9, "descending": 2.0e-9}[mode]
    energy = (-32.8 if ref.solid == "LiH" else -1102.0) + 4.0 * (scale - 1.0) ** 2
    energy += mode_offset
    base.write_file(out, synthetic_output(ref.solid, energy, continuation=mode != "cold"))
    base.write_file(restart, f"fixture restart {ref.solid} {mode} {scale:.5f}\n")
    parent_manifest_record = runner.artifact(parent_manifest) if parent_manifest else None
    parent_restart_record = runner.artifact(parent_restart) if parent_restart else None
    signature = base.job_signature(
        Path("/fixture/cp2k"),
        inp,
        executable_identity={
            "path": "/fixture/cp2k",
            "sha256": str(campaign["cp2k_executable_sha256"]),
        },
        command_contract={
            "driver": "cp2k",
            "diagnostic": "lc12_gxtb_multistart",
            "policy_id": plan["policy_id"],
            "plan_sha256": plan_sha256,
            "solid": ref.solid,
            "mesh": "k444",
            "scale": scale,
            "mode": mode,
            "parent_restart_sha256": (
                parent_restart_record["sha256"] if parent_restart_record else None
            ),
            "parent_manifest_sha256": (
                parent_manifest_record["sha256"] if parent_manifest_record else None
            ),
            "omp_threads": 1,
            "production_eligible": False,
        },
        campaign_fingerprint=campaign,
    )
    base.write_job_stamp(out, signature, completed=True, return_code=0)
    payload = {
        "schema_version": 1,
        "diagnostic": "lc12_gxtb_multistart",
        "production_eligible": False,
        "completed": True,
        "return_code": 0,
        "solid": ref.solid,
        "mesh": "k444",
        "scale": scale,
        "lattice_a_A": ref.a_exp * scale,
        "mode": mode,
        "plan_sha256": plan_sha256,
        "policy_id": plan["policy_id"],
        "campaign_state_at_execution": "qualification_pending",
        "campaign_identity": campaign,
        "job_signature": signature,
        "parent_candidate_manifest": parent_manifest_record,
        "parent_wfn_restart": parent_restart_record,
        "input": runner.artifact(inp),
        "output": runner.artifact(out),
        "wfn_restart": runner.artifact(restart),
        "selection_status": "unclassified; never an EOS point",
    }
    base.write_file(manifest, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return manifest, restart


def build_solid_fixture(
    root: Path,
    solid: str,
    campaign: dict[str, object],
    plan: dict[str, object],
    plan_sha256: str,
) -> None:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    record = plan["solids"][solid]  # type: ignore[index]
    scales = tuple(float(value) for value in record["scales"])
    cold: dict[float, tuple[Path, Path]] = {}
    for scale in scales:
        cold[scale] = fixture_candidate(
            root,
            ref=refs[solid],
            scale=scale,
            mode="cold",
            campaign=campaign,
            plan=plan,
            plan_sha256=plan_sha256,
            parent_manifest=None,
            parent_restart=None,
        )
    parent_manifest, parent_restart = cold[scales[0]]
    for scale in scales[1:]:
        parent_manifest, parent_restart = fixture_candidate(
            root,
            ref=refs[solid],
            scale=scale,
            mode="ascending",
            campaign=campaign,
            plan=plan,
            plan_sha256=plan_sha256,
            parent_manifest=parent_manifest,
            parent_restart=parent_restart,
        )
    parent_manifest, parent_restart = cold[scales[-1]]
    for scale in reversed(scales[:-1]):
        parent_manifest, parent_restart = fixture_candidate(
            root,
            ref=refs[solid],
            scale=scale,
            mode="descending",
            campaign=campaign,
            plan=plan,
            plan_sha256=plan_sha256,
            parent_manifest=parent_manifest,
            parent_restart=parent_restart,
        )


class Goldzak12GXTBMultiStartTests(unittest.TestCase):
    def test_versioned_plan_has_exact_candidate_counts(self) -> None:
        plan, digest = runner.load_plan(runner.DEFAULT_PLAN)
        self.assertEqual(len(digest), 64)
        self.assertEqual(plan["required_cp2k_ancestor"], "c92cc08b45378b85150447011b5a4bb552f5b797")
        lih = plan["solids"]["LiH"]["scales"]  # type: ignore[index]
        mgo = plan["solids"]["MgO"]["scales"]  # type: ignore[index]
        self.assertEqual((len(lih), 3 * len(lih) - 2), (18, 52))
        self.assertEqual((len(mgo), 3 * len(mgo) - 2), (20, 58))

    def test_cp2k_ancestor_gate_rejects_a_pre_required_branch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)
            subprocess.run(["git", "init", "-q", str(source)], check=True)
            subprocess.run(
                ["git", "-C", str(source), "config", "user.email", "fixture@example.test"],
                check=True,
            )
            subprocess.run(
                ["git", "-C", str(source), "config", "user.name", "Fixture"],
                check=True,
            )
            marker = source / "marker"
            marker.write_text("before\n")
            subprocess.run(["git", "-C", str(source), "add", "marker"], check=True)
            subprocess.run(
                ["git", "-C", str(source), "commit", "-q", "-m", "before"],
                check=True,
            )
            before = subprocess.check_output(
                ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
            ).strip()
            marker.write_text("required\n")
            subprocess.run(["git", "-C", str(source), "add", "marker"], check=True)
            subprocess.run(
                ["git", "-C", str(source), "commit", "-q", "-m", "required"],
                check=True,
            )
            required = subprocess.check_output(
                ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
            ).strip()
            runner.require_cp2k_ancestor(source, required)
            subprocess.run(
                ["git", "-C", str(source), "checkout", "-q", "-b", "old", before],
                check=True,
            )
            with self.assertRaisesRegex(ValueError, "not descended"):
                runner.require_cp2k_ancestor(source, required)

    def test_input_contract_rejects_full_grid_and_pins_restart_logging(self) -> None:
        ref = next(ref for ref in base.REFERENCES if ref.solid == "LiH")
        text = runner.multistart_input(ref, 1.0, "fixture", Path("parent-RESTART.kp"))
        self.assertIn("FULL_GRID F", text)
        self.assertIn("LOG_PRINT_KEY T", text)
        self.assertIn("SCF_GUESS RESTART", text)
        with self.assertRaisesRegex(ValueError, "SPGLIB-reduced mesh contract"):
            runner.validate_multistart_input(text.replace("FULL_GRID F", "FULL_GRID T"), True)

    def test_native_mixer_residual_parser_uses_last_fdiis_iteration(self) -> None:
        text = synthetic_output("LiH", -32.8, continuation=False)
        residual, label = runner.final_native_mixer_residual(text)
        self.assertEqual(label, "FDIIS")
        self.assertAlmostEqual(float(residual), 5.0e-10)

    def test_physical_gate_rejects_charge_inverted_lih(self) -> None:
        plan, _ = runner.load_plan(runner.DEFAULT_PLAN)
        policy = plan["classification_policy"]
        solid_policy = plan["solids"]["LiH"]  # type: ignore[index]
        atoms = [
            {"element": "Li", "atomic_population": 1.2, "net_charge": -0.2}
            for _ in range(4)
        ] + [
            {"element": "H", "atomic_population": 0.8, "net_charge": 0.2}
            for _ in range(4)
        ]
        gates, _ = classify.physical_gates(
            "LiH",
            {
                "mulliken_atoms": atoms,
                "electron_count": 8,
                "fermi_energy_hartree": -0.2,
                "mo_occupations_printed": True,
            },
            policy,
            solid_policy,
        )
        self.assertFalse(gates["electronegativity_polarity"])
        self.assertTrue(gates["nonnegative_atomic_populations"])

    def test_full_lih_fixture_yields_hash_pinned_quadratic_path(self) -> None:
        plan, plan_sha256 = runner.load_plan(runner.DEFAULT_PLAN)
        campaign = fake_campaign()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            build_solid_fixture(root, "LiH", campaign, plan, plan_sha256)
            result = classify.classify_solid(
                root,
                "LiH",
                plan,
                plan_sha256,
                campaign,
                "qualification_pending",
            )
            self.assertEqual(result["expected_candidate_count"], 52)
            self.assertEqual(result["numerically_valid_candidate_count"], 52)
            self.assertEqual(result["physically_valid_candidate_count"], 52)
            self.assertTrue(result["full_continuous_physical_path"])
            self.assertEqual(result["fit"]["fit_status"], "quadratic")  # type: ignore[index]
            self.assertEqual(len(result["selected_path"]), 18)  # type: ignore[arg-type]

            *_, target_manifest = runner.candidate_paths(
                root, "LiH", "cold", 1.0
            )
            _, foreign_output, _, _ = runner.candidate_paths(
                root, "LiH", "cold", 1.02
            )
            tampered = json.loads(target_manifest.read_text())
            tampered["output"] = runner.artifact(foreign_output)
            target_manifest.write_text(json.dumps(tampered, indent=2, sort_keys=True) + "\n")
            rejected = classify.classify_solid(
                root,
                "LiH",
                plan,
                plan_sha256,
                campaign,
                "qualification_pending",
            )
            self.assertEqual(rejected["numerically_valid_candidate_count"], 51)
            cold = next(
                candidate
                for candidate in rejected["candidates_by_scale"]["1.00000"]  # type: ignore[index]
                if candidate["mode"] == "cold"
            )
            self.assertIn("not bound to the candidate path", " ".join(cold["issues"]))

    def test_missing_candidate_fails_closed(self) -> None:
        plan, plan_sha256 = runner.load_plan(runner.DEFAULT_PLAN)
        campaign = fake_campaign()
        with tempfile.TemporaryDirectory() as tmp:
            result = classify.classify_solid(
                Path(tmp),
                "MgO",
                plan,
                plan_sha256,
                campaign,
                "qualification_pending",
            )
            self.assertFalse(result["full_continuous_physical_path"])
            self.assertIn("no physically valid candidate", str(result["failure"]))

    def test_resume_revalidates_restart_hash_before_skip(self) -> None:
        plan, plan_sha256 = runner.load_plan(runner.DEFAULT_PLAN)
        ref = next(ref for ref in base.REFERENCES if ref.solid == "LiH")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k"
            cp2k.write_bytes(b"fixture executable")
            campaign = fake_campaign(base.sha256(cp2k))

            def successful_run(
                executable: Path, inp: Path, out: Path, threads: int
            ) -> int:
                del executable, threads
                base.write_file(out, synthetic_output("LiH", -32.8, continuation=False))
                restart = out.with_name(out.stem + "-RESTART.kp")
                base.write_file(restart, f"restart for {inp.name}\n")
                return 0

            arguments = {
                "root": root,
                "ref": ref,
                "scale": 1.0,
                "mode": "cold",
                "restart": None,
                "parent_manifest": None,
                "cp2k": cp2k,
                "campaign_identity": campaign,
                "campaign_state": "qualification_pending",
                "plan": plan,
                "plan_sha256": plan_sha256,
                "threads": 1,
                "retry_failed": False,
            }
            with patch.object(base, "run_cp2k", side_effect=successful_run) as mocked:
                first = runner.run_candidate(**arguments)
                self.assertTrue(first["completed"])
                self.assertEqual(mocked.call_count, 1)
                _, _, restart, _ = runner.candidate_paths(root, "LiH", "cold", 1.0)
                restart.write_text("tampered restart\n")
                with self.assertRaisesRegex(RuntimeError, "WFN restart artifact hash mismatch"):
                    runner.run_candidate(**arguments)
                self.assertEqual(mocked.call_count, 1)

    def test_explicit_failed_retry_archives_before_replacement(self) -> None:
        plan, plan_sha256 = runner.load_plan(runner.DEFAULT_PLAN)
        ref = next(ref for ref in base.REFERENCES if ref.solid == "LiH")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k"
            cp2k.write_bytes(b"fixture executable")
            campaign = fake_campaign(base.sha256(cp2k))
            calls = 0

            def staged_run(
                executable: Path, inp: Path, out: Path, threads: int
            ) -> int:
                nonlocal calls
                del executable, threads
                calls += 1
                if calls == 1:
                    base.write_file(out, " SCF run NOT converged\n")
                    return 1
                base.write_file(out, synthetic_output("LiH", -32.8, continuation=False))
                restart = out.with_name(out.stem + "-RESTART.kp")
                base.write_file(restart, f"restart for {inp.name}\n")
                return 0

            arguments = {
                "root": root,
                "ref": ref,
                "scale": 1.0,
                "mode": "cold",
                "restart": None,
                "parent_manifest": None,
                "cp2k": cp2k,
                "campaign_identity": campaign,
                "campaign_state": "qualification_pending",
                "plan": plan,
                "plan_sha256": plan_sha256,
                "threads": 1,
            }
            with patch.object(base, "run_cp2k", side_effect=staged_run):
                failed = runner.run_candidate(**arguments, retry_failed=False)
                self.assertFalse(failed["completed"])
                preserved = runner.run_candidate(**arguments, retry_failed=False)
                self.assertFalse(preserved["completed"])
                self.assertEqual(calls, 1)
                completed = runner.run_candidate(**arguments, retry_failed=True)
                self.assertTrue(completed["completed"])
                self.assertEqual(calls, 2)
            archive = (
                runner.candidate_paths(root, "LiH", "cold", 1.0)[3].parent
                / "attempt_archive"
                / "attempt_001"
            )
            self.assertTrue((archive / "archive_manifest.json").is_file())
            self.assertTrue((archive / "candidate_manifest.json").is_file())

    def test_campaign_lock_is_nonblocking(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with runner.campaign_lock(root):
                with self.assertRaisesRegex(RuntimeError, "already locked"):
                    with runner.campaign_lock(root):
                        self.fail("nested campaign lock unexpectedly succeeded")

    def test_classifier_uses_relocatable_root_snapshots(self) -> None:
        plan, plan_sha256 = runner.load_plan(runner.DEFAULT_PLAN)
        campaign = fake_campaign()
        with tempfile.TemporaryDirectory() as tmp:
            top = Path(tmp)
            root = top / "campaign"
            root.mkdir()
            build_solid_fixture(root, "LiH", campaign, plan, plan_sha256)
            build_solid_fixture(root, "MgO", campaign, plan, plan_sha256)
            plan_snapshot = root / "plan_snapshot.json"
            runner.atomic_write_text(plan_snapshot, runner.DEFAULT_PLAN.read_text())
            relocated_plan = top / "relocated-plan.json"
            runner.atomic_write_text(relocated_plan, runner.DEFAULT_PLAN.read_text())
            build_snapshot = root / "build_manifest_snapshot.json"
            build = {
                "campaign_id": campaign["campaign_id"],
                "campaign_state": "qualification_pending",
                "cp2k": {
                    "binary_sha256": campaign["cp2k_executable_sha256"],
                    "loaded_library_sha256": campaign[
                        "cp2k_loaded_library_sha256"
                    ],
                    "cmake_cache_sha256": campaign["cp2k_cmake_cache_sha256"],
                    "reported_revision": campaign[
                        "cp2k_embedded_source_revision"
                    ],
                    "revision": campaign["cp2k_source_revision"],
                },
                "save_tblite": {
                    "cli_sha256": campaign["save_tblite_executable_sha256"],
                    "revision": campaign["save_tblite_source_revision"],
                    "static_library_sha256": campaign[
                        "save_tblite_library_sha256"
                    ],
                    "cmake_cache_sha256": campaign[
                        "save_tblite_cmake_cache_sha256"
                    ],
                },
                "fetched_dependencies": {"tblite": "d" * 40},
            }
            runner.atomic_write_text(
                build_snapshot,
                json.dumps(build, indent=2, sort_keys=True) + "\n",
            )
            campaign_manifest = {
                "schema_version": 1,
                "diagnostic": "lc12_gxtb_multistart_campaign",
                "production_eligible": False,
                "completed": True,
                "campaign_state_at_execution": "qualification_pending",
                "campaign_identity": campaign,
                "plan": runner.artifact(plan_snapshot),
                "plan_sha256": plan_sha256,
                "build_manifest": runner.artifact(build_snapshot),
                "required_cp2k_ancestor": plan["required_cp2k_ancestor"],
                "systems": [],
            }
            runner.atomic_write_text(
                root / "campaign_manifest.json",
                json.dumps(campaign_manifest, indent=2, sort_keys=True) + "\n",
            )
            with patch.object(
                sys,
                "argv",
                [
                    "classify_gxtb_multistart_branches.py",
                    "--campaign-root",
                    str(root),
                    "--plan",
                    str(relocated_plan),
                ],
            ):
                self.assertEqual(classify.main(), 0)
            selection = json.loads((root / "branch_selection.json").read_text())
            self.assertTrue(selection["completed"])
            self.assertEqual(selection["plan"]["path"], str(plan_snapshot.resolve()))


if __name__ == "__main__":
    unittest.main()
