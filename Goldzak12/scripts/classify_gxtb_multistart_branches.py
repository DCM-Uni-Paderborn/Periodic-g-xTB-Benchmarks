#!/usr/bin/env python3
"""Classify the hash-pinned LiH/MgO g-XTB multi-start branch map.

This script never promotes calculations into the LC12 EOS tree.  It verifies
every expected cold/ascending/descending candidate, applies numerical and
physical gates, and finds the lowest-energy path that is continuous across the
entire versioned scale grid.  Missing evidence fails closed.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping

import diagnose_gxtb_wfn_hysteresis as wfn
import run_goldzak12_benchmark as base
import run_goldzak12_eos_benchmark as eos
import run_gxtb_multistart_branches as runner


SCHEMA_VERSION = 1
MODE_RANK = {"cold": 0, "ascending": 1, "descending": 2}


def verify_artifact(record: object, label: str) -> tuple[Path | None, str | None]:
    if not isinstance(record, Mapping):
        return None, f"missing {label} artifact record"
    path = Path(str(record.get("path", "")))
    expected = str(record.get("sha256", ""))
    if not path.is_file():
        return None, f"missing {label} artifact {path}"
    observed = base.sha256(path)
    if not expected or observed != expected:
        return None, f"{label} artifact hash mismatch for {path}"
    return path, None


def mean_element_charges(atoms: list[dict[str, object]]) -> dict[str, float]:
    grouped: dict[str, list[float]] = {}
    for atom in atoms:
        grouped.setdefault(str(atom["element"]), []).append(float(atom["net_charge"]))
    return {
        element: sum(values) / len(values)
        for element, values in sorted(grouped.items())
    }


def physical_gates(
    solid: str,
    observables: Mapping[str, object],
    policy: Mapping[str, object],
    solid_policy: Mapping[str, object],
) -> tuple[dict[str, bool], dict[str, object]]:
    atoms_raw = observables.get("mulliken_atoms", [])
    atoms = [dict(atom) for atom in atoms_raw] if isinstance(atoms_raw, list) else []
    populations = [
        float(atom["atomic_population"])
        for atom in atoms
        if "atomic_population" in atom
    ]
    charges = [float(atom["net_charge"]) for atom in atoms if "net_charge" in atom]
    means = mean_element_charges(atoms) if len(atoms) == 8 else {}
    cation = str(solid_policy["cation"])
    anion = str(solid_policy["anion"])
    element_counts = {
        element: sum(atom.get("element") == element for atom in atoms)
        for element in {cation, anion}
    }
    polarity_tolerance = float(policy["polarity_charge_tolerance_e"])
    equivalent_spreads: dict[str, float] = {}
    for element in sorted(means):
        values = [float(atom["net_charge"]) for atom in atoms if atom["element"] == element]
        equivalent_spreads[element] = max(values) - min(values)
    fermi = observables.get("fermi_energy_hartree")
    expected_electrons = int(policy["required_electron_count"][solid])  # type: ignore[index]
    max_charge = float(policy["maximum_absolute_charge_e"][solid])  # type: ignore[index]
    gates = {
        "eight_mulliken_atoms": len(atoms) == 8,
        "expected_four_plus_four_stoichiometry": (
            set(means) == {cation, anion}
            and element_counts == {cation: 4, anion: 4}
        ),
        "nonnegative_atomic_populations": (
            len(populations) == 8
            and all(
                math.isfinite(value)
                and value >= float(policy["minimum_atomic_population_e"])
                for value in populations
            )
        ),
        "neutral_cell": (
            len(charges) == 8
            and abs(sum(charges)) <= float(policy["total_charge_e_max"])
        ),
        "bounded_atomic_charges": (
            len(charges) == 8
            and all(math.isfinite(value) and abs(value) <= max_charge for value in charges)
        ),
        "equivalent_atoms_agree": (
            bool(equivalent_spreads)
            and all(
                spread <= float(policy["equivalent_atom_charge_spread_e_max"])
                for spread in equivalent_spreads.values()
            )
        ),
        "electronegativity_polarity": (
            cation in means
            and anion in means
            and means[cation] >= -polarity_tolerance
            and means[anion] <= polarity_tolerance
        ),
        "expected_electron_count": observables.get("electron_count") == expected_electrons,
        "finite_fermi_energy": isinstance(fermi, float) and math.isfinite(fermi),
        "mo_occupations_printed": observables.get("mo_occupations_printed") is True,
    }
    descriptors = {
        "mean_element_charge_e": means,
        "equivalent_atom_charge_spread_e": equivalent_spreads,
        "total_mulliken_charge_e": sum(charges) if len(charges) == 8 else None,
        "minimum_atomic_population_e": min(populations) if len(populations) == 8 else None,
        "maximum_absolute_atomic_charge_e": (
            max(abs(value) for value in charges) if len(charges) == 8 else None
        ),
        "fermi_energy_hartree": fermi,
        "electron_count": observables.get("electron_count"),
    }
    return gates, descriptors


def candidate_record(
    manifest_path: Path,
    *,
    expected_campaign: Mapping[str, object],
    expected_campaign_state: str,
    expected_plan_sha256: str,
    plan: Mapping[str, object],
    solid: str,
    scale: float,
    mode: str,
    expected_input: Path,
    expected_output: Path,
    expected_restart: Path,
    expected_parent_manifest: Path | None,
    expected_parent_restart: Path | None,
    expected_execution_contract: Mapping[str, object] | None = None,
) -> dict[str, object]:
    record: dict[str, object] = {
        "solid": solid,
        "scale": scale,
        "mode": mode,
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": None,
        "numerically_valid": False,
        "physically_valid": False,
        "issues": [],
    }
    issues: list[str] = []
    if not manifest_path.is_file():
        issues.append("missing candidate manifest")
        record["issues"] = issues
        return record
    record["manifest_sha256"] = base.sha256(manifest_path)
    try:
        manifest = json.loads(manifest_path.read_text())
    except json.JSONDecodeError as exc:
        issues.append(f"invalid candidate manifest: {exc}")
        record["issues"] = issues
        return record
    declarations = {
        "schema_version": manifest.get("schema_version"),
        "diagnostic": manifest.get("diagnostic"),
        "production_eligible": manifest.get("production_eligible"),
        "solid": manifest.get("solid"),
        "mesh": manifest.get("mesh"),
        "scale": manifest.get("scale"),
        "mode": manifest.get("mode"),
        "plan_sha256": manifest.get("plan_sha256"),
        "campaign_state_at_execution": manifest.get("campaign_state_at_execution"),
        "campaign_identity": manifest.get("campaign_identity"),
    }
    expected_declarations = {
        "schema_version": SCHEMA_VERSION,
        "diagnostic": "lc12_gxtb_multistart",
        "production_eligible": False,
        "solid": solid,
        "mesh": "k444",
        "scale": scale,
        "mode": mode,
        "plan_sha256": expected_plan_sha256,
        "campaign_state_at_execution": expected_campaign_state,
        "campaign_identity": expected_campaign,
    }
    if declarations != expected_declarations:
        issues.append("candidate declarations differ from campaign/plan target")
    inp, input_issue = verify_artifact(manifest.get("input"), "input")
    out, output_issue = verify_artifact(manifest.get("output"), "output")
    restart, restart_issue = verify_artifact(manifest.get("wfn_restart"), "WFN restart")
    issues.extend(issue for issue in (input_issue, output_issue, restart_issue) if issue)
    for observed, expected, label in (
        (inp, expected_input, "input"),
        (out, expected_output, "output"),
        (restart, expected_restart, "WFN restart"),
    ):
        if observed is not None and observed.resolve() != expected.resolve():
            issues.append(
                f"{label} artifact is not bound to the candidate path: "
                f"{observed} != {expected}"
            )
    parent_manifest_record = manifest.get("parent_candidate_manifest")
    parent_restart_record = manifest.get("parent_wfn_restart")
    if mode == "cold":
        if parent_manifest_record is not None or parent_restart_record is not None:
            issues.append("cold candidate unexpectedly declares a parent")
    else:
        parent_manifest, parent_manifest_issue = verify_artifact(
            parent_manifest_record, "parent candidate manifest"
        )
        parent_restart, parent_restart_issue = verify_artifact(
            parent_restart_record, "parent WFN restart"
        )
        issues.extend(
            issue for issue in (parent_manifest_issue, parent_restart_issue) if issue
        )
        if parent_manifest is not None and (
            expected_parent_manifest is None
            or parent_manifest.resolve() != expected_parent_manifest.resolve()
        ):
            issues.append("continuation parent candidate is not the preceding chain node")
        if parent_restart is not None and (
            expected_parent_restart is None
            or parent_restart.resolve() != expected_parent_restart.resolve()
        ):
            issues.append("continuation WFN is not the preceding chain node restart")
    expected_parent_manifest_sha = (
        base.sha256(expected_parent_manifest)
        if expected_parent_manifest is not None and expected_parent_manifest.is_file()
        else None
    )
    expected_parent_restart_sha = (
        base.sha256(expected_parent_restart)
        if expected_parent_restart is not None and expected_parent_restart.is_file()
        else None
    )
    signature = manifest.get("job_signature")
    execution = plan["execution"]
    assert isinstance(execution, Mapping)
    expected_command_contract = {
        "driver": "cp2k",
        "diagnostic": "lc12_gxtb_multistart",
        "policy_id": plan["policy_id"],
        "plan_sha256": expected_plan_sha256,
        "solid": solid,
        "mesh": "k444",
        "scale": scale,
        "mode": mode,
        "parent_restart_sha256": expected_parent_restart_sha,
        "parent_manifest_sha256": expected_parent_manifest_sha,
        "omp_threads": int(execution["omp_threads_per_job"]),
        "production_eligible": False,
    }
    if isinstance(signature, Mapping):
        executable = str(signature.get("executable", ""))
        expected_signature = {
            "schema_version": 1,
            "executable": executable,
            "executable_sha256": expected_campaign["cp2k_executable_sha256"],
            "input": str(expected_input.resolve()),
            "input_sha256": base.sha256(expected_input) if expected_input.is_file() else "",
            "command_contract": expected_command_contract,
            "campaign_identity": expected_campaign,
        }
        if not executable or signature != expected_signature:
            issues.append("candidate job signature differs from the complete expected contract")
    else:
        expected_signature = None
        issues.append("candidate lacks a hash-pinned job signature")
    if inp is not None:
        try:
            runner.validate_multistart_input(inp.read_text(), continuation=mode != "cold")
        except ValueError as exc:
            issues.append(str(exc))
    if out is not None:
        stamp_issue = base.completed_stamp_campaign_issue(
            out,
            dict(expected_campaign),
            executable_role="cp2k",
            require_completed=True,
        )
        if stamp_issue:
            issues.append(stamp_issue)
        if expected_signature is None or not base.job_stamp_matches(
            out, expected_signature
        ):
            issues.append("output stamp differs from the complete candidate job signature")
        execution_record = manifest.get("execution_provenance")
        if expected_execution_contract is None:
            if execution_record is not None:
                issues.append("direct candidate unexpectedly records MPI/affinity provenance")
        else:
            _, execution_issue = runner.classify_execution_artifact(
                execution_record,
                out,
                expected_execution_contract,
            )
            if execution_issue:
                issues.append(execution_issue)
    if issues or out is None or restart is None:
        record["issues"] = issues
        return record
    numerical, diagnostics = runner.numerical_gates(
        out,
        restart,
        continuation=mode != "cold",
    )
    if manifest.get("completed") is not True:
        issues.append("candidate manifest is not completed")
    if not all(value is True for value in numerical.values()):
        issues.append("one or more recomputed numerical gates failed")
    record["numerical_gates"] = numerical
    record["numerical_diagnostics"] = diagnostics
    record["numerically_valid"] = not issues
    if issues:
        record["issues"] = issues
        return record
    observables = diagnostics["observables"]
    assert isinstance(observables, Mapping)
    policy = plan["classification_policy"]
    solids = plan["solids"]
    assert isinstance(policy, Mapping) and isinstance(solids, Mapping)
    solid_policy = solids[solid]
    assert isinstance(solid_policy, Mapping)
    gates, descriptors = physical_gates(solid, observables, policy, solid_policy)
    record["physical_gates"] = gates
    record["descriptors"] = descriptors
    record["energy_hartree"] = observables.get("total_energy_extrapolated_t0_hartree")
    record["physically_valid"] = all(gates.values())
    if not record["physically_valid"]:
        issues.append("one or more physical gates failed")
    record["issues"] = issues
    return record


def charge_distance(left: Mapping[str, object], right: Mapping[str, object]) -> float:
    left_desc = left["descriptors"]
    right_desc = right["descriptors"]
    assert isinstance(left_desc, Mapping) and isinstance(right_desc, Mapping)
    left_charge = left_desc["mean_element_charge_e"]
    right_charge = right_desc["mean_element_charge_e"]
    assert isinstance(left_charge, Mapping) and isinstance(right_charge, Mapping)
    if set(left_charge) != set(right_charge):
        return math.inf
    return max(abs(float(left_charge[key]) - float(right_charge[key])) for key in left_charge)


def assign_clusters(
    candidates: list[dict[str, object]], policy: Mapping[str, object]
) -> list[dict[str, object]]:
    energy_tolerance = float(policy["cluster_energy_hartree_max"])
    charge_tolerance = float(policy["cluster_mean_element_charge_e_max"])
    clusters: list[list[dict[str, object]]] = []
    valid = sorted(
        (candidate for candidate in candidates if candidate["numerically_valid"]),
        key=lambda item: (
            float(item.get("energy_hartree", math.inf)),
            MODE_RANK[str(item["mode"])],
        ),
    )
    for candidate in valid:
        destination = None
        for index, cluster in enumerate(clusters):
            representative = cluster[0]
            energy_delta = abs(
                float(candidate["energy_hartree"]) - float(representative["energy_hartree"])
            )
            if (
                energy_delta <= energy_tolerance
                and charge_distance(candidate, representative) <= charge_tolerance
            ):
                destination = index
                break
        if destination is None:
            clusters.append([candidate])
            destination = len(clusters) - 1
        else:
            clusters[destination].append(candidate)
        candidate["cluster_id"] = destination + 1
    return [
        {
            "cluster_id": index + 1,
            "members": [str(item["mode"]) for item in cluster],
            "minimum_energy_hartree": min(float(item["energy_hartree"]) for item in cluster),
            "all_members_physically_valid": all(item["physically_valid"] for item in cluster),
        }
        for index, cluster in enumerate(clusters)
    ]


def adjacent(
    left: Mapping[str, object],
    right: Mapping[str, object],
    policy: Mapping[str, object],
) -> bool:
    if charge_distance(left, right) > float(policy["adjacent_mean_element_charge_change_e_max"]):
        return False
    left_desc = left["descriptors"]
    right_desc = right["descriptors"]
    assert isinstance(left_desc, Mapping) and isinstance(right_desc, Mapping)
    return abs(
        float(left_desc["fermi_energy_hartree"])
        - float(right_desc["fermi_energy_hartree"])
    ) <= float(policy["adjacent_fermi_change_hartree_max"])


def lowest_continuous_path(
    by_scale: Mapping[float, list[dict[str, object]]],
    scales: tuple[float, ...],
    policy: Mapping[str, object],
) -> tuple[list[dict[str, object]] | None, str | None]:
    choices = {
        scale: [candidate for candidate in by_scale[scale] if candidate["physically_valid"]]
        for scale in scales
    }
    missing = [scale for scale in scales if not choices[scale]]
    if missing:
        return None, "no physically valid candidate at scales " + ", ".join(
            f"{scale:.5f}" for scale in missing
        )
    costs: list[dict[int, tuple[float, int | None]]] = []
    for index, scale in enumerate(scales):
        scale_choices = choices[scale]
        minimum = min(float(item["energy_hartree"]) for item in scale_choices)
        current: dict[int, tuple[float, int | None]] = {}
        if index == 0:
            for node, candidate in enumerate(scale_choices):
                current[node] = (float(candidate["energy_hartree"]) - minimum, None)
        else:
            previous_choices = choices[scales[index - 1]]
            previous_costs = costs[-1]
            for node, candidate in enumerate(scale_choices):
                options = [
                    (previous_costs[parent][0], parent)
                    for parent, previous in enumerate(previous_choices)
                    if parent in previous_costs and adjacent(previous, candidate, policy)
                ]
                if options:
                    prior_cost, parent = min(options, key=lambda value: (value[0], value[1]))
                    current[node] = (
                        prior_cost + float(candidate["energy_hartree"]) - minimum,
                        parent,
                    )
        if not current:
            return None, f"no continuous physical edge into scale {scale:.5f}"
        costs.append(current)
    last_choices = choices[scales[-1]]
    node = min(
        costs[-1],
        key=lambda item: (
            costs[-1][item][0],
            MODE_RANK[str(last_choices[item]["mode"])],
            str(last_choices[item]["manifest_path"]),
        ),
    )
    selected: list[dict[str, object]] = []
    for index in range(len(scales) - 1, -1, -1):
        selected.append(choices[scales[index]][node])
        parent = costs[index][node][1]
        if parent is not None:
            node = parent
    selected.reverse()
    return selected, None


def classify_solid(
    root: Path,
    solid: str,
    plan: Mapping[str, object],
    plan_sha256: str,
    campaign_identity: Mapping[str, object],
    campaign_state: str,
    execution_contract: Mapping[str, object] | None = None,
) -> dict[str, object]:
    solids = plan["solids"]
    policy = plan["classification_policy"]
    assert isinstance(solids, Mapping) and isinstance(policy, Mapping)
    solid_policy = solids[solid]
    assert isinstance(solid_policy, Mapping)
    scales = tuple(float(value) for value in solid_policy["scales"])
    by_scale: dict[float, list[dict[str, object]]] = {scale: [] for scale in scales}
    expected: list[tuple[str, float]] = [("cold", scale) for scale in scales]
    expected += [("ascending", scale) for scale in scales[1:]]
    expected += [("descending", scale) for scale in scales[:-1]]
    for mode, scale in expected:
        inp, out, restart, manifest = runner.candidate_paths(root, solid, mode, scale)
        index = scales.index(scale)
        expected_parent_manifest = None
        expected_parent_restart = None
        if mode == "ascending":
            parent_scale = scales[index - 1]
            parent_mode = "cold" if index == 1 else "ascending"
            _, _, expected_parent_restart, expected_parent_manifest = runner.candidate_paths(
                root, solid, parent_mode, parent_scale
            )
        elif mode == "descending":
            parent_scale = scales[index + 1]
            parent_mode = "cold" if index == len(scales) - 2 else "descending"
            _, _, expected_parent_restart, expected_parent_manifest = runner.candidate_paths(
                root, solid, parent_mode, parent_scale
            )
        by_scale[scale].append(
            candidate_record(
                manifest,
                expected_campaign=campaign_identity,
                expected_campaign_state=campaign_state,
                expected_plan_sha256=plan_sha256,
                plan=plan,
                solid=solid,
                scale=scale,
                mode=mode,
                expected_input=inp,
                expected_output=out,
                expected_restart=restart,
                expected_parent_manifest=expected_parent_manifest,
                expected_parent_restart=expected_parent_restart,
                expected_execution_contract=execution_contract,
            )
        )
    clusters = {
        f"{scale:.5f}": assign_clusters(by_scale[scale], policy) for scale in scales
    }
    selected, failure = lowest_continuous_path(by_scale, scales, policy)
    refs = {ref.solid: ref for ref in base.REFERENCES}
    fit: dict[str, object] | None = None
    if selected is not None:
        points = [
            (
                refs[solid].a_exp * float(candidate["scale"]),
                float(candidate["scale"]),
                float(candidate["energy_hartree"]),
                True,
            )
            for candidate in selected
        ]
        fit = eos.fit_gxtb_eos(points)
        if fit.get("fit_status") != "quadratic":
            failure = f"selected full path failed EOS gate: {fit.get('fit_status')}"
    return {
        "solid": solid,
        "mesh": "k444",
        "requested_scales": list(scales),
        "expected_candidate_count": len(expected),
        "numerically_valid_candidate_count": sum(
            candidate["numerically_valid"]
            for candidates in by_scale.values()
            for candidate in candidates
        ),
        "physically_valid_candidate_count": sum(
            candidate["physically_valid"]
            for candidates in by_scale.values()
            for candidate in candidates
        ),
        "clusters_by_scale": clusters,
        "candidates_by_scale": {
            f"{scale:.5f}": by_scale[scale] for scale in scales
        },
        "full_continuous_physical_path": failure is None,
        "failure": failure,
        "selected_path": (
            [
                {
                    "scale": candidate["scale"],
                    "mode": candidate["mode"],
                    "cluster_id": candidate.get("cluster_id"),
                    "energy_hartree": candidate["energy_hartree"],
                    "descriptors": candidate["descriptors"],
                    "candidate_manifest": {
                        "path": candidate["manifest_path"],
                        "sha256": candidate["manifest_sha256"],
                    },
                }
                for candidate in selected
            ]
            if selected is not None
            else None
        ),
        "fit": fit,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--plan", type=Path, default=runner.DEFAULT_PLAN)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        plan, plan_sha256 = runner.load_plan(args.plan)
        root = args.campaign_root.resolve(strict=True)
        with runner.campaign_lock(root):
            campaign_path = root / "campaign_manifest.json"
            campaign = json.loads(campaign_path.read_text())
            declarations = {
                "schema_version": campaign.get("schema_version"),
                "diagnostic": campaign.get("diagnostic"),
                "production_eligible": campaign.get("production_eligible"),
                "plan_sha256": campaign.get("plan_sha256"),
                "required_cp2k_ancestor": campaign.get("required_cp2k_ancestor"),
                "execution_contract": campaign.get("execution_contract"),
            }
            expected_declarations = {
                "schema_version": SCHEMA_VERSION,
                "diagnostic": "lc12_gxtb_multistart_campaign",
                "production_eligible": False,
                "plan_sha256": plan_sha256,
                "required_cp2k_ancestor": plan["required_cp2k_ancestor"],
                "execution_contract": campaign.get("execution_contract"),
            }
            if declarations != expected_declarations:
                raise ValueError("campaign declarations differ from the classification plan")
            if not isinstance(campaign.get("completed"), bool):
                raise ValueError("campaign completion flag is missing or invalid")
            campaign_completed = campaign["completed"] is True
            campaign_identity = campaign["campaign_identity"]
            base.validate_campaign_identity(campaign_identity)
            campaign_state = str(campaign.get("campaign_state_at_execution", ""))
            if campaign_state not in runner.CAMPAIGN_STATES:
                raise ValueError("campaign has an invalid execution state")
            execution_contract_raw = campaign.get("execution_contract")
            if execution_contract_raw is not None and not isinstance(
                execution_contract_raw, Mapping
            ):
                raise ValueError("campaign has an invalid execution contract")
            execution_contract = execution_contract_raw
            expected_execution_sha = (
                runner.execution.canonical_sha256(dict(execution_contract))
                if execution_contract is not None
                else None
            )
            if campaign.get("execution_contract_sha256") != expected_execution_sha:
                raise ValueError("campaign execution-contract hash is missing or invalid")

            campaign_plan, campaign_plan_issue = verify_artifact(
                campaign.get("plan"), "campaign plan snapshot"
            )
            if campaign_plan_issue:
                raise ValueError(campaign_plan_issue)
            expected_plan_snapshot = root / "plan_snapshot.json"
            if (
                campaign_plan is None
                or campaign_plan.resolve() != expected_plan_snapshot.resolve()
            ):
                raise ValueError("campaign plan is not bound to its root snapshot")
            snapshot_plan, snapshot_plan_sha256 = runner.load_plan(campaign_plan)
            if snapshot_plan_sha256 != plan_sha256 or snapshot_plan != plan:
                raise ValueError("campaign plan snapshot and requested plan differ")

            build_snapshot, build_issue = verify_artifact(
                campaign.get("build_manifest"), "build-manifest snapshot"
            )
            if build_issue:
                raise ValueError(build_issue)
            expected_build_snapshot = root / "build_manifest_snapshot.json"
            if (
                build_snapshot is None
                or build_snapshot.resolve() != expected_build_snapshot.resolve()
            ):
                raise ValueError("campaign build manifest is not bound to its root snapshot")
            build_declarations = json.loads(build_snapshot.read_text())
            snapshot_identity = base.campaign_identity_from_manifest(
                build_declarations,
                build_snapshot,
                allowed_campaign_states=(campaign_state,),
            )
            if snapshot_identity != campaign_identity:
                raise ValueError("build snapshot identity differs from the campaign identity")

            systems = [
                classify_solid(
                    root,
                    solid,
                    plan,
                    plan_sha256,
                    campaign_identity,
                    campaign_state,
                    execution_contract,
                )
                for solid in ("LiH", "MgO")
            ]
            completed = campaign_completed and all(
                system["full_continuous_physical_path"] for system in systems
            )
            payload = {
                "schema_version": SCHEMA_VERSION,
                "classification": "lc12_gxtb_multistart_branch_selection",
                "production_eligible": False,
                "automatic_promotion": False,
                "completed": completed,
                "campaign_execution_completed": campaign_completed,
                "campaign_identity": campaign_identity,
                "campaign_state_at_execution": campaign_state,
                "execution_contract": execution_contract,
                "campaign_manifest": {
                    "path": str(campaign_path.resolve()),
                    "sha256": base.sha256(campaign_path),
                },
                "plan": {
                    "path": str(campaign_plan.resolve()),
                    "sha256": base.sha256(campaign_plan),
                    "canonical_sha256": plan_sha256,
                    "policy_id": plan["policy_id"],
                },
                "selection_policy": (
                    "lowest total relative energy among paths that pass every "
                    "numerical/physical gate and the adjacent charge/Fermi "
                    "continuity thresholds"
                ),
                "systems": systems,
                "manual_production_review_prerequisites": {
                    "both_full_paths_and_quadratic_fits": completed,
                    "current_build_manifest_must_be_production_ready_and_identity_identical": True,
                    "selection_and_all_candidate_hashes_must_be_revalidated": True,
                    "k555_final_single_points_are_not_started_by_this_script": True,
                },
            }
            output = args.output or (root / "branch_selection.json")
            runner.atomic_write_text(
                output,
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
            )
    except (
        OSError,
        KeyError,
        TypeError,
        ValueError,
        RuntimeError,
        json.JSONDecodeError,
    ) as exc:
        parser.error(str(exc))
    print(f"completed={completed} selection={output}", flush=True)
    return 0 if completed else 2


if __name__ == "__main__":
    raise SystemExit(main())
