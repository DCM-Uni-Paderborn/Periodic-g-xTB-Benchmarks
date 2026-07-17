#!/usr/bin/env python3
"""Fail-closed verifier for mixer symmetry-star storage qualification."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MATRIX = json.loads((ROOT / "test_matrix.json").read_text())
INPUT_ROOT = ROOT / "test_inputs"
RUN_ROOT = ROOT / os.environ.get("RUN_ROOT", "runs_v2_exact_binding")
SUMMARY_PATH = ROOT / os.environ.get(
    "SUMMARY_FILE", f"{RUN_ROOT.name}_summary.tsv"
)
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
QUALIFY_RE = re.compile(
    rf"GXTB-QUALIFICATION_ONLY MIXER-STAR iter=(\d+)"
    rf"\s+denseCov=\s*({FLOAT})\s+streamCov=\s*({FLOAT})"
    rf"\s+streamRoundtrip=\s*({FLOAT})\s+covDelta=\s*({FLOAT})"
    rf"\s+denseFullComplex=(\d+)\s+streamedPeakComplex=(\d+)"
)
STREAM_RE = re.compile(
    rf"GXTB-MIXER-STAR-STREAMED denseFullComplexAvoided=(\d+),"
    rf" peakComplex=(\d+), covariance=\s*({FLOAT}), roundtrip=\s*({FLOAT})"
)
MODE_RE = re.compile(
    r"GXTB-KGROUP-PARTIAL-ROOT groups=(\d+), nred=(\d+), nfull=(\d+), batch=(\d+);"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_ordered_pe_list(value: str) -> tuple[int, ...]:
    fields = value.split(",")
    if not fields or any(not field.strip().isdigit() for field in fields):
        raise RuntimeError("invalid literal ordered PE list")
    cpus = tuple(int(field.strip()) for field in fields)
    if len(set(cpus)) != len(cpus):
        raise RuntimeError("ordered PE list contains duplicate CPUs")
    return cpus


def reported_binding_rank_ids(text: str) -> list[int]:
    return sorted(
        {
            int(value)
            for value in re.findall(
                r"\b(?:MCW\s+)?rank\s+(\d+)\s+bound\b",
                text,
                flags=re.IGNORECASE,
            )
        }
    )


def expanded_cases() -> list[dict]:
    result = []
    for case in MATRIX["cases"]:
        for ranks in case["ranks"]:
            result.append({**case, "ranks": int(ranks)})
    return result


def validate_rank_process_provenance(
    item: dict, rank: int, expected_executable: str
) -> None:
    """Require one immutable, terminally resolved Linux task per rank."""
    required_fields = {
        "pid",
        "rank",
        "raw_rank",
        "rank_identity_source",
        "rank_observation_status",
        "is_cp2k_rank",
        "executable",
        "arguments",
        "sample_count",
        "process_starttime",
        "observed_process_starttimes",
        "process_identity_status",
        "snapshot_consistency_status",
        "stat_state",
        "state",
        "observed_process_states",
        "observed_rank_observation_statuses",
        "process_terminally_confirmed",
        "process_terminal_confirmation",
        "current_sample_matches_assigned_singleton",
        "affinity_violation_ever",
        "rank_identity_changed_ever",
        "process_starttime_changed_ever",
        "process_reappeared_after_terminal_ever",
        "executable_changed_ever",
        "cpu_mask_changed_during_sample_ever",
        "snapshot_unavailable_ever",
    }
    if not required_fields.issubset(item):
        raise RuntimeError("missing CP2K rank process provenance")
    pid = item.get("pid")
    raw_rank = item.get("raw_rank")
    starttime = item.get("process_starttime")
    sample_count = item.get("sample_count")
    statuses = item.get("observed_rank_observation_statuses")
    states = item.get("observed_process_states")
    confirmation = item.get("process_terminal_confirmation")
    arguments = item.get("arguments")
    argument_zero_matches = False
    if (
        isinstance(arguments, list)
        and arguments
        and all(isinstance(argument, str) for argument in arguments)
    ):
        try:
            argument_zero_matches = (
                str(Path(arguments[0]).resolve(strict=True))
                == expected_executable
            )
        except (FileNotFoundError, OSError):
            argument_zero_matches = False
    terminal_confirmation = confirmation == "process_disappeared" or (
        isinstance(confirmation, str)
        and confirmation in {"terminal_state_Z", "terminal_state_X"}
    )
    unavailable = item.get("rank_environment_unavailable_ever") is True
    sticky_failure_fields = (
        "affinity_violation_ever",
        "rank_identity_changed_ever",
        "process_starttime_changed_ever",
        "process_reappeared_after_terminal_ever",
        "executable_changed_ever",
        "cpu_mask_changed_during_sample_ever",
        "snapshot_unavailable_ever",
    )
    if (
        item.get("is_cp2k_rank") is not True
        or not isinstance(pid, int)
        or isinstance(pid, bool)
        or pid <= 0
        or item.get("rank") != rank
        or item.get("executable") != expected_executable
        or not argument_zero_matches
        or not isinstance(sample_count, int)
        or isinstance(sample_count, bool)
        or sample_count < 1
        or not isinstance(starttime, int)
        or isinstance(starttime, bool)
        or starttime < 0
        or item.get("observed_process_starttimes") != [starttime]
        or any(
            not isinstance(value, int) or isinstance(value, bool)
            for value in item.get("observed_process_starttimes", [])
        )
        or any(item.get(field) is not False for field in sticky_failure_fields)
        or item.get("process_terminally_confirmed") is not True
        or not terminal_confirmation
        or item.get("process_identity_status")
        not in {"stable", "terminal_state", "disappeared_after_sample"}
        or item.get("snapshot_consistency_status")
        not in {"consistent", "process_disappeared"}
        or not isinstance(item.get("stat_state"), str)
        or not re.fullmatch(r"[RSDZTtXIWP]", str(item.get("stat_state")))
        or not isinstance(item.get("state"), str)
        or not re.fullmatch(
            r"[RSDZTtXIWP](?:\s+\([^\r\n]*\))?", str(item.get("state"))
        )
        or not isinstance(states, list)
        or not states
        or any(
            not isinstance(state, str)
            or not re.fullmatch(r"[RSDZTtXIWP](?:\s+\([^\r\n]*\))?", state)
            for state in states
        )
        or len(states) != len(set(states))
        or len(states) > sample_count
        or item.get("state") not in states
        or not isinstance(statuses, list)
        or not statuses
        or any(not isinstance(status, str) for status in statuses)
        or len(statuses) != len(set(statuses))
        or len(statuses) > sample_count
        or statuses[0] != "explicit"
        or item.get("rank_observation_status") != statuses[-1]
    ):
        raise RuntimeError("invalid CP2K rank process provenance")

    if unavailable:
        if (
            raw_rank is not None
            or item.get("rank_identity_source")
            != "pending_terminal_environment_loss"
        ):
            raise RuntimeError("invalid CP2K rank process provenance")
    elif (
        raw_rank != rank
        or isinstance(raw_rank, bool)
        or statuses != ["explicit"]
        or item.get("rank_observation_status") != "explicit"
        or item.get("rank_identity_source") != "explicit_environment"
    ):
        raise RuntimeError("invalid CP2K rank process provenance")

    identity_status = item.get("process_identity_status")
    consistency_status = item.get("snapshot_consistency_status")
    stat_state = item.get("stat_state")
    if (
        identity_status == "disappeared_after_sample"
        and (
            confirmation != "process_disappeared"
            or consistency_status != "process_disappeared"
        )
    ) or (
        identity_status == "terminal_state"
        and (
            confirmation != f"terminal_state_{stat_state}"
            or consistency_status != "consistent"
        )
    ) or (
        identity_status == "stable" and consistency_status != "consistent"
    ) or (
        identity_status == "stable" and stat_state in {"Z", "X"}
    ):
        raise RuntimeError("inconsistent CP2K rank terminal provenance")


def validate_rank_environment_evidence(item: dict, expected_mask: str) -> None:
    required_fields = {
        "rank_environment_unavailable_ever",
        "rank_environment_unavailable_pending",
        "rank_environment_terminally_confirmed",
        "rank_environment_terminal_confirmation",
        "rank_environment_unavailable_sample_count",
        "rank_environment_events",
    }
    if not required_fields.issubset(item):
        raise RuntimeError("missing terminal rank-environment evidence")
    unavailable_sample_count = item.get(
        "rank_environment_unavailable_sample_count"
    )
    if not isinstance(unavailable_sample_count, int) or isinstance(
        unavailable_sample_count, bool
    ):
        raise RuntimeError("invalid terminal rank-environment sample count")
    unavailable = item.get("rank_environment_unavailable_ever") is True
    if not unavailable:
        if (
            item.get("rank_environment_unavailable_ever") is not False
            or item.get("rank_environment_unavailable_pending") is not False
            or item.get("rank_environment_terminally_confirmed") is not False
            or item.get("rank_environment_terminal_confirmation") is not None
            or item.get("rank_environment_unavailable_sample_count") != 0
            or item.get("rank_environment_events") != []
        ):
            raise RuntimeError("inconsistent terminal rank-environment evidence")
        return

    pid = item.get("pid")
    sample_count = item.get("sample_count")
    starttime = item.get("process_starttime")
    statuses = item.get("observed_rank_observation_statuses")
    events = item.get("rank_environment_events")
    confirmation = item.get("rank_environment_terminal_confirmation")
    terminal_confirmation = confirmation == "process_disappeared" or (
        isinstance(confirmation, str)
        and confirmation in {"terminal_state_Z", "terminal_state_X"}
    )
    if (
        not isinstance(pid, int)
        or isinstance(pid, bool)
        or not isinstance(sample_count, int)
        or isinstance(sample_count, bool)
        or not isinstance(starttime, int)
        or isinstance(starttime, bool)
        or item.get("observed_process_starttimes") != [starttime]
        or item.get("process_starttime_changed_ever") is not False
        or item.get("rank_environment_unavailable_pending") is not False
        or item.get("rank_environment_terminally_confirmed") is not True
        or not terminal_confirmation
        or confirmation != item.get("process_terminal_confirmation")
        or item.get("rank_identity_source")
        != "pending_terminal_environment_loss"
        or item.get("raw_rank") is not None
        or not isinstance(statuses, list)
        or not statuses
        or statuses[0] != "explicit"
        or not isinstance(events, list)
        or not events
        or item.get("rank_environment_unavailable_sample_count") != len(events)
    ):
        raise RuntimeError("invalid terminal rank-environment evidence")

    expected_statuses = ["explicit"]
    previous_sample_index = 0
    for event in events:
        if not isinstance(event, dict):
            raise RuntimeError("invalid terminal rank-environment evidence")
        event_status = event.get("environment_status")
        event_sample_index = event.get("sample_index")
        if (
            not isinstance(event_sample_index, int)
            or isinstance(event_sample_index, bool)
            or event_sample_index <= previous_sample_index
            or event_sample_index > sample_count
            or event_sample_index < 2
            or not isinstance(event.get("pid"), int)
            or isinstance(event.get("pid"), bool)
            or event.get("pid") != pid
            or not isinstance(event.get("process_starttime"), int)
            or isinstance(event.get("process_starttime"), bool)
            or event.get("process_starttime") != starttime
            or event.get("cpus_allowed_list") != expected_mask
            or event_status not in {"environment_empty", "environment_unreadable"}
            or event.get("terminal_resolution") != confirmation
            or event.get("state") not in item.get("observed_process_states", [])
        ):
            raise RuntimeError("invalid terminal rank-environment event sequence")
        previous_sample_index = event_sample_index
        if event_status not in expected_statuses:
            expected_statuses.append(str(event_status))
    if statuses != expected_statuses:
        raise RuntimeError("invalid terminal rank-environment status history")
    if (
        events[-1].get("sample_index") != sample_count
        or events[-1].get("state") != item.get("state")
    ):
        raise RuntimeError("invalid terminal rank-environment final sample")


def revalidated_rank_evidence(
    metadata: dict, ranks: int, expected_cpus: tuple[int, ...]
) -> list[dict]:
    children = metadata.get("observed_child_processes")
    if not isinstance(children, list) or not all(
        isinstance(item, dict) for item in children
    ):
        raise RuntimeError("invalid child-process affinity evidence")
    by_pid: dict[int, dict] = {}
    for item in children:
        pid = item.get("pid")
        if not isinstance(pid, int) or isinstance(pid, bool) or pid in by_pid:
            raise RuntimeError("invalid or duplicate observed child PID")
        by_pid[pid] = item

    samples = metadata.get("concurrent_duplicate_rank_samples")
    if not isinstance(samples, list):
        raise RuntimeError("missing concurrent-rank sample evidence")
    derived_duplicate_ids: set[int] = set()
    previous_sample_index = 0
    for sample in samples:
        if not isinstance(sample, dict):
            raise RuntimeError("invalid concurrent-rank sample evidence")
        sample_index = sample.get("sample_index")
        rank_pid_groups = sample.get("rank_pid_groups")
        if (
            not isinstance(sample_index, int)
            or isinstance(sample_index, bool)
            or sample_index <= previous_sample_index
            or not isinstance(rank_pid_groups, list)
            or not rank_pid_groups
        ):
            raise RuntimeError("invalid concurrent-rank sample evidence")
        previous_sample_index = sample_index
        ranks_in_sample: list[int] = []
        for group in rank_pid_groups:
            if not isinstance(group, dict):
                raise RuntimeError("invalid concurrent-rank sample evidence")
            rank = group.get("rank")
            pids = group.get("pids")
            if (
                not isinstance(rank, int)
                or isinstance(rank, bool)
                or not 0 <= rank < ranks
                or not isinstance(pids, list)
                or len(pids) < 2
                or any(
                    not isinstance(pid, int) or isinstance(pid, bool) for pid in pids
                )
                or pids != sorted(set(pids))
                or any(
                    pid not in by_pid or by_pid[pid].get("rank") != rank
                    for pid in pids
                )
            ):
                raise RuntimeError("invalid concurrent-rank sample evidence")
            ranks_in_sample.append(rank)
            derived_duplicate_ids.add(rank)
        if ranks_in_sample != sorted(set(ranks_in_sample)):
            raise RuntimeError("invalid concurrent-rank sample evidence")

    duplicate_ids = metadata.get("concurrent_duplicate_rank_ids_ever")
    if (
        not isinstance(duplicate_ids, list)
        or any(
            not isinstance(rank, int)
            or isinstance(rank, bool)
            or not 0 <= rank < ranks
            for rank in duplicate_ids
        )
        or duplicate_ids != sorted(derived_duplicate_ids)
        or metadata.get("concurrent_duplicate_rank_processes_ever")
        is not bool(derived_duplicate_ids)
    ):
        raise RuntimeError("inconsistent concurrent-rank summary")

    groups: dict[int, list[dict]] = {}
    expected_executable = str(Path(str(metadata.get("cp2k", ""))).resolve())
    for item in children:
        rank = item.get("rank")
        if not isinstance(rank, int) or isinstance(rank, bool) or not 0 <= rank < ranks:
            raise RuntimeError("invalid observed MPI rank identity")
        expected_mask = str(expected_cpus[rank])
        validate_rank_process_provenance(item, rank, expected_executable)
        validate_rank_environment_evidence(item, expected_mask)
        if (
            item.get("observed_rank_ids") != [rank]
            or any(
                not isinstance(value, int) or isinstance(value, bool)
                for value in item.get("observed_rank_ids", [])
            )
            or item.get("observed_cpu_masks") != [expected_mask]
            or item.get("cpus_allowed_list") != expected_mask
            or not isinstance(item.get("sample_count"), int)
            or isinstance(item.get("sample_count"), bool)
            or item["sample_count"] < 1
            or item.get("current_sample_matches_assigned_singleton") is not True
        ):
            raise RuntimeError("invalid singleton rank-affinity child history")
        groups.setdefault(rank, []).append(item)
    if sorted(groups) != list(range(ranks)):
        raise RuntimeError("logical MPI rank set is incomplete")
    if any(len(generations) != 1 for generations in groups.values()):
        raise RuntimeError("multiple MPI rank PID generations are not scaling-eligible")

    recomputed: list[dict] = []
    for rank in range(ranks):
        generations = sorted(groups[rank], key=lambda item: int(item["pid"]))
        mask_history = sorted(
            {
                str(mask)
                for item in generations
                for mask in item.get("observed_cpu_masks", [])
            }
        )
        exact = (
            mask_history == [str(expected_cpus[rank])]
            and rank not in derived_duplicate_ids
        )
        canonical = max(
            generations,
            key=lambda item: (
                int(item.get("sample_count", 0)),
                -int(item["pid"]),
            ),
        )
        recomputed.append({
            "rank": rank,
            "pid": int(canonical["pid"]),
            "pid_generations": [int(item["pid"]) for item in generations],
            "cpus_allowed_list": mask_history[0] if len(mask_history) == 1 else "",
            "observed_cpu_masks": mask_history,
            "current_sample_matches_assigned_singleton": exact,
            "affinity_violation_ever": not exact,
            "concurrent_duplicate_pid_ever": rank in derived_duplicate_ids,
        })
    return recomputed


def checked_run(case: dict, variant: str) -> tuple[dict, str]:
    stem = f"{case['name']}_p{case['ranks']}_{variant.lower()}"
    run_dir = RUN_ROOT / stem
    required = [run_dir / name for name in ("run.json", "returncode.txt", "cp2k.out")]
    if not all(path.is_file() for path in required):
        raise RuntimeError(f"missing result file in {run_dir}")
    if (run_dir / "returncode.txt").read_text().strip() != "0":
        raise RuntimeError(f"nonzero return code: {run_dir}")
    metadata = json.loads((run_dir / "run.json").read_text())
    if metadata.get("returncode") != 0 or metadata.get("variant") != variant:
        raise RuntimeError(f"metadata variant/result mismatch: {run_dir}")
    if metadata.get("case") != case["name"] or metadata.get("ranks") != case["ranks"]:
        raise RuntimeError(f"metadata case/rank mismatch: {run_dir}")
    input_path = (INPUT_ROOT / case["input"]).resolve()
    if metadata.get("input") != str(input_path):
        raise RuntimeError(f"input path mismatch: {run_dir}")
    if metadata.get("input_sha256") != sha256(input_path):
        raise RuntimeError(f"input hash mismatch: {run_dir}")
    if metadata.get("working_directory") != str(run_dir.resolve()):
        raise RuntimeError(f"working-directory isolation mismatch: {run_dir}")
    cp2k = Path(metadata.get("cp2k", ""))
    if not cp2k.is_file() or metadata.get("cp2k_sha256") != sha256(cp2k):
        raise RuntimeError(f"CP2K executable hash mismatch: {run_dir}")
    cp2k_lib = Path(metadata.get("cp2k_lib", ""))
    if not cp2k_lib.is_file() or metadata.get("cp2k_lib_sha256") != sha256(cp2k_lib):
        raise RuntimeError(f"CP2K shared-library hash mismatch: {run_dir}")
    affinity = metadata.get("affinity_proof")
    if not isinstance(affinity, list) or len(affinity) != case["ranks"]:
        raise RuntimeError(f"missing live rank-affinity proof: {run_dir}")
    schema = metadata.get("schema_version", metadata.get("schema"))
    if schema == 1:
        # Historical shared-mask records remain valid for numerical comparison
        # only. Their wall times are never scaling evidence.
        first_cpu, last_cpu = map(int, metadata["cpu_set"].split("-"))
        expected_cpus = set(range(first_cpu, last_cpu + 1))
        if any(
            set(item.get("cpus_allowed", [])) != expected_cpus
            or item.get("processor") not in expected_cpus
            for item in affinity
        ):
            raise RuntimeError(f"invalid historical affinity proof: {run_dir}")
        if not (run_dir / "cp2k.err").is_file():
            raise RuntimeError(f"missing historical stderr: {run_dir}")
        metadata["timing_classification"] = "legacy_timing_non_scaling"
    elif schema == 2:
        expected_cpus = parse_ordered_pe_list(str(metadata.get("ordered_pe_list", "")))
        if len(expected_cpus) != case["ranks"]:
            raise RuntimeError(f"ordered PE-list/rank mismatch: {run_dir}")
        try:
            recomputed_affinity = revalidated_rank_evidence(
                metadata, case["ranks"], expected_cpus
            )
        except RuntimeError as error:
            raise RuntimeError(f"rank-evidence revalidation failed: {run_dir}") from error
        if affinity != recomputed_affinity:
            raise RuntimeError(f"derived rank-affinity summary mismatch: {run_dir}")
        if any(
            not isinstance(item, dict)
            or not isinstance(item.get("rank"), int)
            or isinstance(item.get("rank"), bool)
            or not isinstance(item.get("pid"), int)
            or isinstance(item.get("pid"), bool)
            or not isinstance(item.get("pid_generations"), list)
            or len(item["pid_generations"]) != 1
            or any(
                not isinstance(pid, int) or isinstance(pid, bool)
                for pid in item["pid_generations"]
            )
            for item in affinity
        ):
            raise RuntimeError(f"invalid derived rank-affinity field types: {run_dir}")
        if [item.get("rank") for item in affinity] != list(range(case["ranks"])):
            raise RuntimeError(f"MPI rank ordering mismatch: {run_dir}")
        if any(
            item.get("cpus_allowed_list") != str(expected_cpus[index])
            or item.get("affinity_violation_ever") is not False
            or item.get("current_sample_matches_assigned_singleton") is not True
            for index, item in enumerate(affinity)
        ):
            raise RuntimeError(f"invalid singleton rank-affinity history: {run_dir}")
        if metadata.get("all_observed_rank_samples_match_ordered_pe_list") is not True:
            raise RuntimeError(f"sticky affinity gate failed: {run_dir}")
        if metadata.get("runtime_affinity_gate") is not True:
            raise RuntimeError(f"runtime affinity gate failed: {run_dir}")
        if metadata.get("cross_process_cpu_reservation_gate") is not True:
            raise RuntimeError(f"CPU reservation gate failed: {run_dir}")
        if metadata.get("live_compute_overlap_preflight_gate") is not True:
            raise RuntimeError(f"live CPU-overlap preflight failed: {run_dir}")
        if metadata.get("live_compute_overlap_preflight_owners") != []:
            raise RuntimeError(f"live CPU-overlap preflight owners recorded: {run_dir}")
        if metadata.get("live_compute_overlap_runtime_gate") is not True:
            raise RuntimeError(f"runtime live CPU-overlap gate failed: {run_dir}")
        if metadata.get("live_compute_overlap_runtime_samples") != []:
            raise RuntimeError(f"runtime live CPU-overlap was recorded: {run_dir}")
        if metadata.get("concurrent_duplicate_rank_processes_ever") is not False:
            raise RuntimeError(f"concurrently live duplicate MPI rank: {run_dir}")
        pid_generations = metadata.get("observed_cp2k_rank_pid_generations")
        if (
            not isinstance(pid_generations, list)
            or len(pid_generations) != case["ranks"]
            or any(
                item.get("pid_generations") != pid_generations[index]
                or not isinstance(pid_generations[index], list)
                or len(pid_generations[index]) != 1
                or any(
                    not isinstance(pid, int) or isinstance(pid, bool)
                    for pid in pid_generations[index]
                )
                for index, item in enumerate(affinity)
            )
        ):
            raise RuntimeError(f"invalid MPI rank PID-generation proof: {run_dir}")
        if metadata.get("observed_cp2k_process_generation_count") != case["ranks"]:
            raise RuntimeError(f"MPI rank process-generation count mismatch: {run_dir}")
        launcher = Path(metadata.get("mpi_launcher", ""))
        if not launcher.is_file() or metadata.get("mpi_launcher_sha256") != sha256(launcher):
            raise RuntimeError(f"MPI launcher hash mismatch: {run_dir}")
        pe_list = ",".join(str(cpu) for cpu in expected_cpus)
        expected_command = [
            str(launcher),
            "--map-by",
            f"pe-list={pe_list}:ordered",
            "--bind-to",
            "core",
            "--report-bindings",
            "-np",
            str(case["ranks"]),
            str(cp2k),
            "-i",
            str(input_path),
        ]
        if metadata.get("command") != expected_command:
            raise RuntimeError(f"MPI command/binding contract mismatch: {run_dir}")
        launcher_log = run_dir / str(metadata.get("launcher_log", ""))
        if (
            not launcher_log.is_file()
            or metadata.get("launcher_log_sha256") != sha256(launcher_log)
            or reported_binding_rank_ids(launcher_log.read_text(errors="replace"))
            != list(range(case["ranks"]))
            or metadata.get("reported_binding_rank_ids") != list(range(case["ranks"]))
        ):
            raise RuntimeError(f"launcher binding report mismatch: {run_dir}")
        if metadata.get("timing_classification") != "production_scaling_eligible":
            raise RuntimeError(f"schema-v2 timing is not scaling eligible: {run_dir}")
    else:
        raise RuntimeError(f"unsupported run schema in {run_dir}: {schema}")
    if metadata.get("output_sha256") != sha256(run_dir / "cp2k.out"):
        raise RuntimeError(f"output hash mismatch: {run_dir}")
    if schema == 1 and metadata.get("stderr_sha256") != sha256(run_dir / "cp2k.err"):
        raise RuntimeError(f"historical stderr hash mismatch: {run_dir}")
    text = (run_dir / "cp2k.out").read_text(errors="replace")
    if text.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"incomplete output: {run_dir}")
    modes = [tuple(map(int, match)) for match in MODE_RE.findall(text)]
    if not modes or any(nfull != case["nfull"] for _, _, nfull, _ in modes):
        raise RuntimeError(f"missing/mismatched partial-root mode marker: {run_dir}")
    return metadata, text


def observables(text: str) -> tuple[float, list[float], list[float]]:
    energies = [float(value) for value in re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
    )]
    forces = [tuple(map(float, row)) for row in re.findall(
        rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
        text,
        re.MULTILINE,
    )]
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    if not energies or not forces or not stress_blocks:
        raise RuntimeError("missing energy, force, or stress observable")
    stress_rows = re.findall(
        rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
        stress_blocks[-1],
        re.MULTILINE,
    )
    if len(stress_rows) != 3:
        raise RuntimeError("malformed stress block")
    force_values = [value for row in forces for value in row]
    stress_values = [float(value) for row in stress_rows for value in row]
    values = [energies[-1], *force_values, *stress_values]
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError("non-finite observable")
    return energies[-1], force_values, stress_values


def max_delta(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise RuntimeError("observable block lengths differ")
    return max(abs(a - b) for a, b in zip(left, right))


def star_residuals(case: dict, variant: str, text: str) -> tuple[float, float, float, int, int]:
    if variant == "STREAMED":
        matches = STREAM_RE.findall(text)
        if len(matches) != 1:
            raise RuntimeError("STREAMED run lacks exactly one selector marker")
        dense_full, streamed_peak, covariance, roundtrip = matches[0]
        dense_residual = float("nan")
        stream_residual = float(covariance)
        roundtrip_residual = float(roundtrip)
    elif variant == "QUALIFY":
        matches = QUALIFY_RE.findall(text)
        if not matches:
            raise RuntimeError("QUALIFY run lacks selector markers")
        dense_residual = max(float(match[1]) for match in matches)
        stream_residual = max(float(match[2]) for match in matches)
        roundtrip_residual = max(float(match[3]) for match in matches)
        covariance_delta = max(float(match[4]) for match in matches)
        if covariance_delta > MATRIX["gates"]["internal_covariance"]:
            raise RuntimeError(f"dense/stream covariance delta failed: {covariance_delta}")
        dense_values = {int(match[5]) for match in matches}
        stream_values = {int(match[6]) for match in matches}
        if len(dense_values) != 1 or len(stream_values) != 1:
            raise RuntimeError("memory counters changed during QUALIFY run")
        dense_full = str(next(iter(dense_values)))
        streamed_peak = str(next(iter(stream_values)))
    else:
        return float("nan"), float("nan"), float("nan"), 0, 0
    dense_full_i = int(dense_full)
    streamed_peak_i = int(streamed_peak)
    nspin = 2 if "UKS" in case["features"] else 1
    if dense_full_i * 3 != streamed_peak_i * nspin * case["nfull"]:
        raise RuntimeError("reported memory counters violate exact allocation formula")
    for value in (stream_residual, roundtrip_residual):
        if not math.isfinite(value) or value < 0.0:
            raise RuntimeError("invalid streamed residual")
    if stream_residual > MATRIX["gates"]["internal_covariance"]:
        raise RuntimeError(f"streamed covariance gate failed: {stream_residual}")
    if roundtrip_residual > MATRIX["gates"]["internal_roundtrip"]:
        raise RuntimeError(f"streamed roundtrip gate failed: {roundtrip_residual}")
    return dense_residual, stream_residual, roundtrip_residual, dense_full_i, streamed_peak_i


def main() -> int:
    rows = []
    executable_hashes = set()
    library_hashes = set()
    for case in expanded_cases():
        runs = {}
        for variant in ("DENSE", "STREAMED", "QUALIFY"):
            metadata, text = checked_run(case, variant)
            executable_hashes.add(metadata["cp2k_sha256"])
            library_hashes.add(metadata["cp2k_lib_sha256"])
            runs[variant] = (metadata, text, observables(text))
        if len(executable_hashes) != 1:
            raise RuntimeError("matrix used more than one CP2K executable")
        if len(library_hashes) != 1:
            raise RuntimeError("matrix used more than one CP2K shared library")
        dense_obs = runs["DENSE"][2]
        streamed_obs = runs["STREAMED"][2]
        qualify_obs = runs["QUALIFY"][2]
        d_energy = max(abs(dense_obs[0] - streamed_obs[0]), abs(dense_obs[0] - qualify_obs[0]))
        d_force = max(max_delta(dense_obs[1], streamed_obs[1]), max_delta(dense_obs[1], qualify_obs[1]))
        d_stress = max(max_delta(dense_obs[2], streamed_obs[2]), max_delta(dense_obs[2], qualify_obs[2]))
        gates = MATRIX["gates"]
        if d_energy > gates["external_energy_Ha"]:
            raise RuntimeError(f"energy gate failed for {case['name']}_p{case['ranks']}: {d_energy}")
        if d_force > gates["external_force_Ha_per_bohr"]:
            raise RuntimeError(f"force gate failed for {case['name']}_p{case['ranks']}: {d_force}")
        if d_stress > gates["external_stress_bar"]:
            raise RuntimeError(f"stress gate failed for {case['name']}_p{case['ranks']}: {d_stress}")
        _, stream_cov, stream_roundtrip, dense_full, streamed_peak = star_residuals(
            case, "STREAMED", runs["STREAMED"][1]
        )
        dense_cov, qualify_cov, qualify_roundtrip, qualify_dense, qualify_peak = star_residuals(
            case, "QUALIFY", runs["QUALIFY"][1]
        )
        if (dense_full, streamed_peak) != (qualify_dense, qualify_peak):
            raise RuntimeError("STREAMED/QUALIFY memory counters differ")
        rows.append({
            "case": case["name"],
            "ranks": case["ranks"],
            "features": ",".join(case["features"]),
            "nfull": case["nfull"],
            "external_dE_Ha": f"{d_energy:.16e}",
            "external_dForce_Ha_per_bohr": f"{d_force:.16e}",
            "external_dStress_bar": f"{d_stress:.16e}",
            "dense_covariance": f"{dense_cov:.16e}",
            "stream_covariance": f"{max(stream_cov, qualify_cov):.16e}",
            "stream_roundtrip": f"{max(stream_roundtrip, qualify_roundtrip):.16e}",
            "dense_full_complex": dense_full,
            "streamed_peak_complex": streamed_peak,
            "timing_classification": (
                "production_scaling_eligible"
                if all(
                    runs[variant][0].get("timing_classification")
                    == "production_scaling_eligible"
                    for variant in ("DENSE", "STREAMED", "QUALIFY")
                )
                else "legacy_timing_non_scaling"
            ),
            "status": "PASS",
        })
    expected = sum(len(case["ranks"]) for case in MATRIX["cases"])
    if len(rows) != expected:
        raise RuntimeError(f"matrix incomplete: {len(rows)} != {expected}")
    with SUMMARY_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)}/{expected} DENSE/STREAMED/QUALIFY triples")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        raise
