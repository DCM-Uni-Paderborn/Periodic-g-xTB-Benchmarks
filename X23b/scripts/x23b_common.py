#!/usr/bin/env python3
"""Shared safeguards for the additive X23b g-xTB workflow."""

from __future__ import annotations

import functools
import hashlib
import json
import os
import re
import subprocess
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping


PUBLISHED_METHODS = ("GFN1", "GFN2")
METHODS = (*PUBLISHED_METHODS, "GXTB")
BUSY_RETURN_CODE = 75
JOB_STAMP_NAME = "job_provenance.json"
GXTB_PROVENANCE_NAME = "build_provenance_gxtb.json"
DEFAULT_CAMPAIGN_MANIFEST = (
    Path(__file__).resolve().parents[2]
    / "campaigns"
    / "gxtb-pbc-v1-20260714"
    / "build_manifest.json"
)
CAMPAIGN_SCHEMA = 1
JOB_STAMP_SCHEMA = 2
CAMPAIGN_IDENTITY_FIELDS = (
    "schema",
    "campaign_id",
    "cp2k_executable_sha256",
    "cp2k_loaded_library_sha256",
    "cp2k_cmake_cache_sha256",
    "cp2k_embedded_source_revision",
    "cp2k_source_revision",
    "save_tblite_executable_sha256",
    "save_tblite_source_revision",
    "save_tblite_library_sha256",
    "save_tblite_cmake_cache_sha256",
    "dependency_lock_sha256",
)


@functools.lru_cache(maxsize=32)
def _sha256_cached(path: str, size: int, mtime_ns: int) -> str:
    del size, mtime_ns
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    """Return a cached content hash which is invalidated by size/mtime changes."""

    path = path.resolve(strict=True)
    stat = path.stat()
    return _sha256_cached(str(path), stat.st_size, stat.st_mtime_ns)


def _fingerprint(payload: Mapping[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def make_campaign_identity(
    *,
    campaign_id: str,
    cp2k_executable_sha256: str,
    cp2k_loaded_library_sha256: str,
    cp2k_cmake_cache_sha256: str,
    cp2k_embedded_source_revision: str,
    cp2k_source_revision: str,
    save_tblite_executable_sha256: str,
    save_tblite_source_revision: str,
    save_tblite_library_sha256: str,
    save_tblite_cmake_cache_sha256: str,
    dependency_lock_sha256: str,
) -> dict[str, object]:
    """Return the path-independent identity frozen for one GXTB campaign."""

    identity: dict[str, object] = {
        "schema": CAMPAIGN_SCHEMA,
        "campaign_id": campaign_id,
        "cp2k_executable_sha256": cp2k_executable_sha256.lower(),
        "cp2k_loaded_library_sha256": cp2k_loaded_library_sha256.lower(),
        "cp2k_cmake_cache_sha256": cp2k_cmake_cache_sha256.lower(),
        "cp2k_embedded_source_revision": cp2k_embedded_source_revision.lower(),
        "cp2k_source_revision": cp2k_source_revision.lower(),
        "save_tblite_executable_sha256": save_tblite_executable_sha256.lower(),
        "save_tblite_source_revision": save_tblite_source_revision.lower(),
        "save_tblite_library_sha256": save_tblite_library_sha256.lower(),
        "save_tblite_cmake_cache_sha256": save_tblite_cmake_cache_sha256.lower(),
        "dependency_lock_sha256": dependency_lock_sha256.lower(),
    }
    identity["fingerprint_sha256"] = _fingerprint(identity)
    validate_campaign_identity(identity)
    return identity


def validate_campaign_identity(identity: Mapping[str, object]) -> None:
    missing = [field for field in CAMPAIGN_IDENTITY_FIELDS if not identity.get(field)]
    if missing:
        raise ValueError("incomplete GXTB campaign identity: " + ", ".join(missing))
    if int(identity["schema"]) != CAMPAIGN_SCHEMA:
        raise ValueError(f"unsupported GXTB campaign schema: {identity['schema']}")
    for field in CAMPAIGN_IDENTITY_FIELDS:
        if field.endswith("_sha256") and not re.fullmatch(
            r"[0-9a-f]{64}", str(identity[field]).lower()
        ):
            raise ValueError(f"invalid GXTB campaign SHA256 field: {field}")
    if not re.fullmatch(r"[0-9a-f]{40}", str(identity["cp2k_source_revision"]).lower()):
        raise ValueError("invalid full CP2K source revision in campaign identity")
    if not re.fullmatch(
        r"[0-9a-f]{7,40}", str(identity["cp2k_embedded_source_revision"]).lower()
    ):
        raise ValueError("invalid embedded CP2K source revision in campaign identity")
    if not str(identity["cp2k_source_revision"]).lower().startswith(
        str(identity["cp2k_embedded_source_revision"]).lower()
    ):
        raise ValueError("embedded CP2K revision is not a prefix of the full source revision")
    if not re.fullmatch(r"[0-9a-f]{40}", str(identity["save_tblite_source_revision"]).lower()):
        raise ValueError("invalid save_tblite source revision in campaign identity")
    core = {field: identity[field] for field in CAMPAIGN_IDENTITY_FIELDS}
    expected = _fingerprint(core)
    if identity.get("fingerprint_sha256") != expected:
        raise ValueError("GXTB campaign fingerprint is internally inconsistent")


def load_campaign_identity(benchmark_root: Path) -> dict[str, object]:
    path = benchmark_root.resolve() / "data" / GXTB_PROVENANCE_NAME
    if not path.is_file():
        raise ValueError(f"GXTB campaign provenance is missing: {path}")
    payload = json.loads(path.read_text())
    identity = payload.get("campaign_identity")
    if not isinstance(identity, dict):
        raise ValueError(f"GXTB campaign identity is missing from {path}")
    validate_campaign_identity(identity)
    manifest_record = payload.get("campaign_manifest")
    if (
        not isinstance(manifest_record, dict)
        or not manifest_record.get("path")
        or not re.fullmatch(
            r"[0-9a-f]{64}", str(manifest_record.get("file_sha256", "")).lower()
        )
    ):
        raise ValueError(f"GXTB campaign manifest record is missing from {path}")
    manifest_path = _resolve_frozen_campaign_manifest(benchmark_root, manifest_record)
    manifest = json.loads(manifest_path.read_text())
    declared = _identity_from_manifest_declarations(manifest, manifest_path)
    observed = {
        field: identity[field]
        for field in (*CAMPAIGN_IDENTITY_FIELDS, "fingerprint_sha256")
    }
    if declared != observed:
        raise ValueError("current campaign manifest build identity differs from X23b provenance")
    return observed


def _resolve_frozen_campaign_manifest(
    benchmark_root: Path, manifest_record: Mapping[str, object]
) -> Path:
    """Resolve a frozen manifest by content after a benchmark checkout is relocated."""

    expected_sha256 = str(manifest_record["file_sha256"]).lower()
    recorded_path = Path(str(manifest_record["path"])).expanduser()
    campaign_id = str(manifest_record.get("campaign_id", "")).strip()
    candidates = [recorded_path]
    if campaign_id:
        candidates.append(
            benchmark_root.resolve().parent
            / "campaigns"
            / campaign_id
            / "build_manifest.json"
        )

    existing: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=True)
        except FileNotFoundError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        existing.append(resolved)
        if sha256_file(resolved) == expected_sha256:
            return resolved

    if existing:
        raise ValueError(
            "current campaign manifest fingerprint differs from the frozen X23b provenance"
        )
    raise ValueError(
        "frozen GXTB campaign manifest is missing at both its recorded and "
        "repository-relative locations"
    )


def require_gxtb_build_artifacts(
    *,
    cp2k: Path | None,
    cp2k_source: Path | None,
    save_tblite: Path | None,
    save_tblite_source: Path | None,
    campaign_manifest: Path | None,
) -> None:
    required = {
        "--cp2k": cp2k,
        "--cp2k-source": cp2k_source,
        "--save-tblite": save_tblite,
        "--save-tblite-source": save_tblite_source,
        "--campaign-manifest": campaign_manifest,
    }
    missing = [option for option, value in required.items() if value is None]
    if missing:
        raise ValueError("GXTB production requires " + ", ".join(missing))


def _identity_from_manifest_declarations(
    manifest: Mapping[str, object],
    manifest_path: Path,
    *,
    allowed_campaign_states: tuple[str, ...] = ("production_ready",),
) -> dict[str, object]:
    campaign_state = str(manifest.get("campaign_state", ""))
    if campaign_state not in allowed_campaign_states:
        if allowed_campaign_states == ("production_ready",):
            raise ValueError(
                f"GXTB campaign {manifest.get('campaign_id', '<unknown>')} is not "
                f"production_ready (state: {campaign_state or '<missing>'})"
            )
        allowed = ", ".join(allowed_campaign_states)
        raise ValueError(
            f"GXTB campaign {manifest.get('campaign_id', '<unknown>')} has state "
            f"{campaign_state or '<missing>'}; allowed state(s): {allowed}"
        )
    cp2k = manifest.get("cp2k")
    save = manifest.get("save_tblite")
    dependencies = manifest.get("fetched_dependencies")
    if not isinstance(cp2k, Mapping) or not isinstance(save, Mapping):
        raise ValueError(f"invalid CP2K/save_tblite records in {manifest_path}")
    if not isinstance(dependencies, Mapping) or not dependencies:
        raise ValueError(f"missing fetched-dependency lock in {manifest_path}")
    return make_campaign_identity(
        campaign_id=str(manifest.get("campaign_id", "")),
        cp2k_executable_sha256=str(cp2k.get("binary_sha256", "")),
        cp2k_loaded_library_sha256=str(cp2k.get("loaded_library_sha256", "")),
        cp2k_cmake_cache_sha256=str(cp2k.get("cmake_cache_sha256", "")),
        cp2k_embedded_source_revision=str(cp2k.get("reported_revision", "")),
        cp2k_source_revision=str(cp2k.get("revision", "")),
        save_tblite_executable_sha256=str(save.get("cli_sha256", "")),
        save_tblite_source_revision=str(save.get("revision", "")),
        save_tblite_library_sha256=str(save.get("static_library_sha256", "")),
        save_tblite_cmake_cache_sha256=str(save.get("cmake_cache_sha256", "")),
        dependency_lock_sha256=_fingerprint(dependencies),
    )


def thread_environment(threads: int) -> dict[str, str]:
    """Create the bounded-thread environment used by every X23b subprocess."""

    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": str(threads),
            "OPENBLAS_NUM_THREADS": str(threads),
            "MKL_NUM_THREADS": str(threads),
            "VECLIB_MAXIMUM_THREADS": "1",
            "OMP_WAIT_POLICY": "PASSIVE",
        }
    )
    return env


def _atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def job_identity(
    input_path: Path,
    cp2k: Path,
    method: str,
    phase: str,
    campaign_identity: Mapping[str, object] | None = None,
    extra_executables: Mapping[str, Path] | None = None,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> dict[str, object]:
    """Hash every artifact which determines whether a GXTB result is reusable."""

    if method == "GXTB":
        if campaign_identity is None:
            raise ValueError("GXTB job identity requires the frozen campaign identity")
        validate_campaign_identity(campaign_identity)
        campaign = {
            field: campaign_identity[field]
            for field in (*CAMPAIGN_IDENTITY_FIELDS, "fingerprint_sha256")
        }
    else:
        campaign = None
    identity: dict[str, object] = {
        "schema": JOB_STAMP_SCHEMA,
        "method": method,
        "phase": phase,
        "input": {
            "path": str(input_path.resolve()),
            "sha256": sha256_file(input_path),
        },
        "cp2k": {
            "path": str(cp2k.resolve(strict=True)),
            "sha256": sha256_file(cp2k),
        },
    }
    if campaign is not None:
        if identity["cp2k"]["sha256"] != campaign["cp2k_executable_sha256"]:
            raise ValueError(
                "CP2K executable does not match the frozen GXTB campaign fingerprint"
            )
        identity["campaign_identity"] = campaign
    if extra_executables:
        identity["extra_executables"] = {
            name: {
                "path": str(path.resolve(strict=True)),
                "sha256": sha256_file(path),
            }
            for name, path in sorted(extra_executables.items())
        }
    if protocol_identity is not None:
        identity["protocol_identity"] = dict(protocol_identity)
    if source_artifacts is not None:
        identity["source_artifacts"] = {
            role: {
                "path": str(path.resolve(strict=True)),
                "sha256": sha256_file(path),
            }
            for role, path in sorted(source_artifacts.items())
        }
    return identity


def write_job_stamp(
    run_dir: Path,
    input_path: Path,
    cp2k: Path,
    method: str,
    phase: str,
    status: str,
    *,
    details: Mapping[str, object] | None = None,
    campaign_identity: Mapping[str, object] | None = None,
    extra_executables: Mapping[str, Path] | None = None,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> Path:
    payload = job_identity(
        input_path,
        cp2k,
        method,
        phase,
        campaign_identity=campaign_identity,
        extra_executables=extra_executables,
        protocol_identity=protocol_identity,
        source_artifacts=source_artifacts,
    )
    payload["status"] = status
    payload["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    if details:
        payload["details"] = dict(details)
    path = run_dir / JOB_STAMP_NAME
    _atomic_json(path, payload)
    return path


def job_stamp_matches(
    run_dir: Path,
    input_path: Path,
    cp2k: Path,
    method: str,
    phase: str,
    *,
    campaign_identity: Mapping[str, object] | None = None,
    extra_executables: Mapping[str, Path] | None = None,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> tuple[bool, str]:
    """Check a GXTB stamp without imposing stamps on frozen published methods."""

    if method != "GXTB":
        return True, "legacy published method"
    path = run_dir / JOB_STAMP_NAME
    if not path.is_file():
        return False, f"missing {JOB_STAMP_NAME}"
    try:
        observed = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return False, f"invalid {JOB_STAMP_NAME}: {exc}"
    try:
        expected = job_identity(
            input_path,
            cp2k,
            method,
            phase,
            campaign_identity=campaign_identity,
            extra_executables=extra_executables,
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )
    except ValueError as exc:
        return False, str(exc)
    for key in (
        "schema",
        "method",
        "phase",
        "input",
        "cp2k",
        "campaign_identity",
        "extra_executables",
        "protocol_identity",
        "source_artifacts",
    ):
        if observed.get(key) != expected.get(key):
            return False, f"{key} fingerprint differs"
    return True, "matching input and executable fingerprints"


def recorded_job_stamp_matches(
    run_dir: Path,
    input_path: Path,
    method: str,
    phase: str,
    output: Path,
    *,
    campaign_identity: Mapping[str, object] | None = None,
    accepted_status_prefixes: tuple[str, ...] = ("converged",),
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> tuple[bool, str]:
    """Validate a completed result when the executable is no longer supplied."""

    if method != "GXTB":
        return True, "legacy published method"
    path = run_dir / JOB_STAMP_NAME
    if not path.is_file():
        return False, f"missing {JOB_STAMP_NAME}"
    try:
        stamp = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return False, f"invalid {JOB_STAMP_NAME}: {exc}"
    if campaign_identity is None:
        return False, "frozen campaign identity was not supplied"
    try:
        validate_campaign_identity(campaign_identity)
    except ValueError as exc:
        return False, str(exc)
    expected_campaign = {
        field: campaign_identity[field]
        for field in (*CAMPAIGN_IDENTITY_FIELDS, "fingerprint_sha256")
    }
    if (
        stamp.get("schema") != JOB_STAMP_SCHEMA
        or stamp.get("method") != method
        or stamp.get("phase") != phase
    ):
        return False, "stamp identity differs"
    if stamp.get("campaign_identity") != expected_campaign:
        return False, "campaign fingerprint differs"
    if stamp.get("protocol_identity") != (
        dict(protocol_identity) if protocol_identity is not None else None
    ):
        return False, "protocol identity differs"
    try:
        expected_sources = (
            {
                role: {
                    "path": str(path.resolve(strict=True)),
                    "sha256": sha256_file(path),
                }
                for role, path in sorted(source_artifacts.items())
            }
            if source_artifacts is not None
            else None
        )
    except (FileNotFoundError, OSError) as exc:
        return False, f"source artifact missing: {exc}"
    if stamp.get("source_artifacts") != expected_sources:
        return False, "source artifact fingerprint differs"
    recorded_input = stamp.get("input", {})
    recorded_cp2k = stamp.get("cp2k", {})
    if (
        not isinstance(recorded_input, dict)
        or recorded_input.get("path") != str(input_path.resolve(strict=True))
        or recorded_input.get("sha256") != sha256_file(input_path)
    ):
        return False, "input fingerprint differs"
    if not isinstance(recorded_cp2k, dict) or not recorded_cp2k.get("sha256"):
        return False, "CP2K fingerprint missing"
    if recorded_cp2k.get("sha256") != expected_campaign["cp2k_executable_sha256"]:
        return False, "CP2K fingerprint differs from the campaign"
    if not str(stamp.get("status", "")).startswith(accepted_status_prefixes):
        return False, "stamp does not record convergence"
    details = stamp.get("details", {})
    if (
        not isinstance(details, dict)
        or details.get("output") != str(output.resolve(strict=True))
        or details.get("output_sha256") != sha256_file(output)
    ):
        return False, "output fingerprint differs"
    return True, "recorded input, executable, and output fingerprints are consistent"


def _command_output(command: list[str], cwd: Path | None = None) -> str:
    process = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    return (process.stdout + process.stderr).strip()


def git_source_metadata(source: Path) -> dict[str, object]:
    source = source.resolve(strict=True)
    return {
        "path": str(source),
        "revision": _command_output(["git", "rev-parse", "HEAD"], source),
        "branch": _command_output(["git", "branch", "--show-current"], source),
        "status": _command_output(["git", "status", "--short"], source),
        "remotes": _command_output(["git", "remote", "-v"], source),
    }


def executable_metadata(executable: Path) -> dict[str, object]:
    executable = executable.resolve(strict=True)
    return {
        "path": str(executable),
        "executable_sha256": sha256_file(executable),
        "version_output": _command_output([str(executable), "--version"]),
    }


def _clean_git_source(source: Path, label: str) -> dict[str, object]:
    metadata = git_source_metadata(source)
    revision = str(metadata["revision"]).lower()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError(f"{label} source is not a valid Git checkout: {source}")
    if metadata["status"]:
        raise ValueError(f"{label} source checkout is dirty; commit or clean it before GXTB production")
    return metadata


def _campaign_from_manifest(
    *,
    cp2k: Path,
    cp2k_source: Path,
    save_tblite: Path,
    save_tblite_source: Path,
    campaign_manifest: Path,
    allowed_campaign_states: tuple[str, ...] = ("production_ready",),
) -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    manifest_path = campaign_manifest.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text())
    declared_identity = _identity_from_manifest_declarations(
        manifest,
        manifest_path,
        allowed_campaign_states=allowed_campaign_states,
    )
    cp2k_expected = manifest.get("cp2k")
    save_expected = manifest.get("save_tblite")
    dependencies = manifest.get("fetched_dependencies")
    if not isinstance(cp2k_expected, dict) or not isinstance(save_expected, dict):
        raise ValueError(f"invalid CP2K/save_tblite records in {manifest_path}")
    if not isinstance(dependencies, dict) or not dependencies:
        raise ValueError(f"missing fetched-dependency lock in {manifest_path}")

    cp2k_record = executable_metadata(cp2k)
    cp2k_source_record = _clean_git_source(cp2k_source, "CP2K")
    embedded = re.search(
        r"^\s*Source code revision\s+([0-9a-f]+)\s*$",
        str(cp2k_record["version_output"]),
        flags=re.I | re.M,
    )
    if embedded is None:
        raise ValueError("CP2K executable does not report an embedded Source code revision")
    embedded_revision = embedded.group(1).lower()
    source_revision = str(cp2k_source_record["revision"]).lower()
    if not source_revision.startswith(embedded_revision):
        raise ValueError(
            "CP2K executable/source revision mismatch: binary embeds "
            f"{embedded_revision}, source checkout is {source_revision}; "
            "perform a clean rebuild before GXTB production"
        )
    cp2k_declared = {
        "executable": str(cp2k_expected.get("binary_sha256", "")).lower(),
        "source": str(cp2k_expected.get("revision", "")).lower(),
        "embedded": str(cp2k_expected.get("reported_revision", "")).lower(),
    }
    cp2k_observed = {
        "executable": str(cp2k_record["executable_sha256"]).lower(),
        "source": source_revision,
        "embedded": embedded_revision,
    }
    if cp2k_declared != cp2k_observed:
        raise ValueError(
            "CP2K executable/source/reported revision differs from the campaign manifest: "
            f"declared={cp2k_declared}, observed={cp2k_observed}"
        )
    if cp2k_expected.get("source_clean") is not True:
        raise ValueError("campaign manifest does not certify a clean CP2K source")
    loaded_library = Path(str(cp2k_expected.get("loaded_library", ""))).resolve(strict=True)
    loaded_library_sha256 = sha256_file(loaded_library)
    if loaded_library_sha256 != str(cp2k_expected.get("loaded_library_sha256", "")).lower():
        raise ValueError("libcp2k artifact differs from the campaign manifest")
    cp2k_record.update(
        {
            "source": cp2k_source_record,
            "source_revision": source_revision,
            "embedded_source_revision": embedded_revision,
            "source_revision_match": True,
            "loaded_library_path": str(loaded_library),
            "loaded_library_sha256": loaded_library_sha256,
            "cmake_cache_sha256": str(cp2k_expected.get("cmake_cache_sha256", "")).lower(),
        }
    )

    save_record = executable_metadata(save_tblite)
    save_source_record = _clean_git_source(save_tblite_source, "save_tblite")
    save_source_revision = str(save_source_record["revision"]).lower()
    library = Path(str(save_expected.get("static_library", ""))).resolve(strict=True)
    if library.name != "libtblite.a":
        raise ValueError(f"expected the static save_tblite archive libtblite.a, got {library}")
    save_declared = {
        "executable": str(save_expected.get("cli_sha256", "")).lower(),
        "source": str(save_expected.get("revision", "")).lower(),
        "library": str(save_expected.get("static_library_sha256", "")).lower(),
    }
    save_observed = {
        "executable": str(save_record["executable_sha256"]).lower(),
        "source": save_source_revision,
        "library": sha256_file(library),
    }
    if save_declared != save_observed:
        raise ValueError(
            "save_tblite CLI/source/libtblite.a differs from the campaign manifest: "
            f"declared={save_declared}, observed={save_observed}"
        )
    if save_expected.get("source_clean") is not True:
        raise ValueError("campaign manifest does not certify a clean save_tblite source")
    reported_version = str(save_expected.get("reported_version", ""))
    if reported_version and reported_version not in str(save_record["version_output"]):
        raise ValueError("save_tblite CLI version differs from the campaign manifest")
    save_record.update(
        {
            "cli_sha256": save_record["executable_sha256"],
            "source": save_source_record,
            "source_revision": save_source_revision,
            "static_library_path": str(library),
            "static_library_sha256": save_observed["library"],
            "cmake_cache_sha256": str(save_expected.get("cmake_cache_sha256", "")).lower(),
            "static_library_role": (
                "explicit build-time archive fingerprint; no runtime-link inference is made"
            ),
        }
    )
    identity = make_campaign_identity(
        campaign_id=str(manifest.get("campaign_id", "")),
        cp2k_executable_sha256=str(cp2k_record["executable_sha256"]),
        cp2k_loaded_library_sha256=loaded_library_sha256,
        cp2k_cmake_cache_sha256=str(cp2k_record["cmake_cache_sha256"]),
        cp2k_embedded_source_revision=embedded_revision,
        cp2k_source_revision=source_revision,
        save_tblite_executable_sha256=str(save_record["executable_sha256"]),
        save_tblite_source_revision=str(save_record["source_revision"]),
        save_tblite_library_sha256=str(save_record["static_library_sha256"]),
        save_tblite_cmake_cache_sha256=str(save_record["cmake_cache_sha256"]),
        dependency_lock_sha256=_fingerprint(dependencies),
    )
    if identity != declared_identity:
        raise ValueError("observed GXTB build artifacts differ from the declared campaign identity")
    manifest_record = {
        "path": str(manifest_path),
        "file_sha256": sha256_file(manifest_path),
        "campaign_id": manifest["campaign_id"],
        "campaign_state": manifest["campaign_state"],
        "authority": "single source of truth for the frozen build artifacts",
    }
    return identity, cp2k_record, save_record, manifest_record


def validate_campaign_artifacts(
    *,
    cp2k: Path,
    cp2k_source: Path,
    save_tblite: Path,
    save_tblite_source: Path,
    campaign_manifest: Path,
    allowed_campaign_states: tuple[str, ...] = ("production_ready",),
) -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    """Validate all declared GXTB build artifacts against one campaign manifest.

    Production callers deliberately retain the default ``production_ready``
    state.  A pre-production validation gate may explicitly admit
    ``validation_in_progress`` without weakening the production boundary.
    """

    return _campaign_from_manifest(
        cp2k=cp2k,
        cp2k_source=cp2k_source,
        save_tblite=save_tblite,
        save_tblite_source=save_tblite_source,
        campaign_manifest=campaign_manifest,
        allowed_campaign_states=allowed_campaign_states,
    )


def declared_campaign_identity(
    campaign_manifest: Path,
    *,
    allowed_campaign_states: tuple[str, ...] = ("production_ready",),
) -> tuple[dict[str, object], str]:
    """Return the path-independent declared identity and current campaign state."""

    manifest_path = campaign_manifest.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text())
    identity = _identity_from_manifest_declarations(
        manifest,
        manifest_path,
        allowed_campaign_states=allowed_campaign_states,
    )
    return identity, str(manifest["campaign_state"])


def _completed(path: Path, *markers: str) -> bool:
    if not path.is_file():
        return False
    text = path.read_text(errors="ignore")
    return all(marker in text for marker in markers)


def _completed_gxtb(
    output: Path,
    input_path: Path,
    phase: str,
    campaign_identity: Mapping[str, object],
    *markers: str,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> bool:
    if not _completed(output, *markers) or not input_path.is_file():
        return False
    valid, _ = recorded_job_stamp_matches(
        output.parent,
        input_path,
        "GXTB",
        phase,
        output,
        campaign_identity=campaign_identity,
        protocol_identity=protocol_identity,
        source_artifacts=source_artifacts,
    )
    return valid


def _project_input(run_dir: Path, project: str) -> Path:
    return run_dir / f"{project.replace('-', '_')}.inp"


def _merge_workflow_paths(current: dict[str, object], updates: Mapping[str, object]) -> None:
    for key, value in updates.items():
        if isinstance(value, Mapping):
            nested = current.get(key)
            if not isinstance(nested, dict):
                nested = {}
                current[key] = nested
            _merge_workflow_paths(nested, value)
        elif isinstance(value, Path):
            current[key] = str(value.resolve())
        else:
            current[key] = value


def update_gxtb_provenance(
    benchmark_root: Path,
    *,
    cp2k: Path | None = None,
    cp2k_source: Path | None = None,
    save_tblite: Path | None = None,
    save_tblite_source: Path | None = None,
    campaign_manifest: Path | None = None,
    workflow_paths: Mapping[str, object] | None = None,
    publication_tables_updated: bool | None = None,
) -> Path:
    """Atomically maintain GXTB-only build and coverage provenance.

    The filename is fixed deliberately: this helper can never overwrite the
    frozen ``build_provenance.json`` for GFN1/GFN2.
    """

    benchmark_root = benchmark_root.resolve()
    path = benchmark_root / "data" / GXTB_PROVENANCE_NAME
    payload: dict[str, object]
    if path.is_file():
        payload = json.loads(path.read_text())
        if payload.get("method") not in (None, "GXTB"):
            raise ValueError(f"refusing to update non-GXTB provenance: {path}")
    else:
        payload = {}

    artifact_values = (cp2k, cp2k_source, save_tblite, save_tblite_source, campaign_manifest)
    candidate_identity: dict[str, object] | None = None
    candidate_cp2k: dict[str, object] | None = None
    candidate_save: dict[str, object] | None = None
    candidate_manifest: dict[str, object] | None = None
    if any(value is not None for value in artifact_values):
        require_gxtb_build_artifacts(
            cp2k=cp2k,
            cp2k_source=cp2k_source,
            save_tblite=save_tblite,
            save_tblite_source=save_tblite_source,
            campaign_manifest=campaign_manifest,
        )
        candidate_identity, candidate_cp2k, candidate_save, candidate_manifest = _campaign_from_manifest(
            cp2k=cp2k,  # type: ignore[arg-type]
            cp2k_source=cp2k_source,  # type: ignore[arg-type]
            save_tblite=save_tblite,  # type: ignore[arg-type]
            save_tblite_source=save_tblite_source,  # type: ignore[arg-type]
            campaign_manifest=campaign_manifest,  # type: ignore[arg-type]
        )

    frozen_identity = payload.get("campaign_identity")
    if frozen_identity is not None:
        if not isinstance(frozen_identity, dict):
            raise ValueError("invalid frozen GXTB campaign identity")
        validate_campaign_identity(frozen_identity)
        if candidate_identity is not None and candidate_identity != frozen_identity:
            raise ValueError(
                "GXTB campaign build fingerprint differs from the frozen campaign; "
                "mixed CP2K/save_tblite artifacts are forbidden"
            )
        campaign_identity = {
            field: frozen_identity[field]
            for field in (*CAMPAIGN_IDENTITY_FIELDS, "fingerprint_sha256")
        }
    elif candidate_identity is not None:
        campaign_identity = candidate_identity
        payload["campaign_identity"] = candidate_identity
    else:
        raise ValueError(
            "GXTB campaign identity has not been established; provide CP2K/source, "
            "save_tblite CLI/source, and the campaign manifest"
        )

    frozen_manifest = payload.get("campaign_manifest")
    if frozen_identity is not None and frozen_manifest is None:
        raise ValueError("frozen GXTB campaign manifest record is missing")
    if frozen_manifest is not None:
        if (
            not isinstance(frozen_manifest, dict)
            or not frozen_manifest.get("path")
            or not re.fullmatch(
                r"[0-9a-f]{64}", str(frozen_manifest.get("file_sha256", "")).lower()
            )
        ):
            raise ValueError("invalid frozen GXTB campaign manifest record")
        frozen_manifest_path = Path(str(frozen_manifest["path"])).resolve(strict=True)
        frozen_manifest_sha = str(frozen_manifest["file_sha256"]).lower()
        if sha256_file(frozen_manifest_path) != frozen_manifest_sha:
            raise ValueError(
                "current campaign manifest fingerprint differs from the frozen X23b provenance"
            )
        if candidate_manifest is not None and (
            Path(str(candidate_manifest["path"])).resolve() != frozen_manifest_path
            or str(candidate_manifest["file_sha256"]).lower() != frozen_manifest_sha
        ):
            raise ValueError(
                "GXTB campaign manifest path/fingerprint differs from the frozen campaign"
            )

    payload.update({"benchmark": "X23b", "method": "GXTB"})
    payload["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    protocol = payload.setdefault("protocol", {})
    if not isinstance(protocol, dict):
        raise ValueError("invalid protocol provenance record")
    protocol.setdefault("model", "METHOD GXTB")
    protocol.setdefault("scc_mixer", "SCC_MIXER TBLITE (native save_tblite complete-Fock potential DIIS)")
    protocol.setdefault(
        "fixed_reference_diagnostic",
        "optional frozen-crystal Gamma/k111/k222/k333 single points; excluded from production defaults",
    )
    protocol["kpoint_mesh_contract"] = (
        "GXTB uses the same SPGLIB-reduced CP2K MACDONALD meshes as GFN1/GFN2 "
        "(SYMMETRY T, FULL_GRID F, SYMMETRY_BACKEND SPGLIB, "
        "SYMMETRY_REDUCTION_METHOD SPGLIB); CP2K expands to and folds from the complete mesh"
    )
    protocol.setdefault(
        "cell_optimization",
        "native Bloch SPGLIB-reduced k222 KEEP_ANGLES with an explicitly recorded source policy",
    )
    protocol["reported_lattice_energy"] = (
        "native Bloch SPGLIB-reduced k333 single point on the matching final GXTB k222 geometry"
    )
    protocol["energy_convergence_check"] = (
        "native Bloch SPGLIB-reduced k444 single point on the matching final GXTB k222 geometry"
    )
    protocol["legacy_full_grid_policy"] = (
        "pre-SPGLIB GXTB full-grid outputs are diagnostics only and are never accepted as production results"
    )
    protocol["resume_policy"] = (
        "GXTB skips and collectors require matching input/output hashes and the complete frozen "
        "CP2K/save_tblite campaign fingerprint; BUSY and STALE_OUTPUT are failures"
    )
    separation = payload.setdefault("separation", {})
    if not isinstance(separation, dict):
        raise ValueError("invalid separation provenance record")
    separation.setdefault("reuses_gfn1_or_gfn2_monomer_energies", False)
    separation.setdefault("reuses_gfn1_or_gfn2_geometries", False)
    separation.setdefault("overwrites_published_gfn1_or_gfn2_data", False)
    separation.setdefault("gxtb_analysis_default", "X23b/data/gxtb_staging")

    if candidate_cp2k is not None and candidate_save is not None and candidate_manifest is not None:
        payload["cp2k"] = candidate_cp2k
        payload["save_tblite"] = candidate_save
        if frozen_manifest is None:
            payload["campaign_manifest"] = candidate_manifest
    cp2k_record = payload.get("cp2k")
    save_record = payload.get("save_tblite")
    if not isinstance(cp2k_record, dict) or not isinstance(save_record, dict):
        raise ValueError("GXTB campaign build records are missing")
    record_checks = {
        "CP2K executable": (
            cp2k_record.get("executable_sha256"),
            campaign_identity["cp2k_executable_sha256"],
        ),
        "CP2K loaded library": (
            cp2k_record.get("loaded_library_sha256"),
            campaign_identity["cp2k_loaded_library_sha256"],
        ),
        "CP2K CMake cache": (
            cp2k_record.get("cmake_cache_sha256"),
            campaign_identity["cp2k_cmake_cache_sha256"],
        ),
        "CP2K embedded source": (
            cp2k_record.get("embedded_source_revision"),
            campaign_identity["cp2k_embedded_source_revision"],
        ),
        "CP2K source": (
            cp2k_record.get("source_revision"),
            campaign_identity["cp2k_source_revision"],
        ),
        "save_tblite executable": (
            save_record.get("executable_sha256"),
            campaign_identity["save_tblite_executable_sha256"],
        ),
        "save_tblite source": (
            save_record.get("source_revision"),
            campaign_identity["save_tblite_source_revision"],
        ),
        "save_tblite static library": (
            save_record.get("static_library_sha256"),
            campaign_identity["save_tblite_library_sha256"],
        ),
        "save_tblite CMake cache": (
            save_record.get("cmake_cache_sha256"),
            campaign_identity["save_tblite_cmake_cache_sha256"],
        ),
    }
    for label, (observed, expected_value) in record_checks.items():
        if observed != expected_value:
            raise ValueError(f"{label} record differs from the frozen GXTB campaign")

    paths = payload.setdefault("workflow_paths", {})
    if not isinstance(paths, dict):
        raise ValueError("invalid workflow_paths provenance record")
    if workflow_paths:
        _merge_workflow_paths(paths, workflow_paths)
    source_policy = str(paths.get("k222_source_policy", "gamma_cellopt_restart"))
    if source_policy not in ("gamma_cellopt_restart", "experimental_reference"):
        raise ValueError(f"unsupported frozen k222 source policy: {source_policy}")
    protocol["k222_source_policy"] = source_policy
    protocol["cell_optimization"] = (
        "native Bloch SPGLIB-reduced k222 KEEP_ANGLES, seeded by the parsed and stamped "
        "experimental shifted-k222 ENERGY_FORCE preflight on the frozen X23 reference structure"
        if source_policy == "experimental_reference"
        else "native Bloch SPGLIB-reduced k222 KEEP_ANGLES, seeded only by the matching GXTB Gamma restart"
    )

    metadata_path = benchmark_root / "data" / "metadata.json"
    systems = []
    if metadata_path.is_file():
        systems = [str(row["id"]) for row in json.loads(metadata_path.read_text())["systems"]]
    expected = len(systems)

    def count_local(phase: str, suffix: str, markers: tuple[str, ...]) -> int:
        count = 0
        for system in systems:
            stem = f"{system}_{suffix}"
            output = benchmark_root / "runs" / phase / "GXTB" / stem / f"{stem}.out"
            input_path = benchmark_root / "inputs" / phase / "GXTB" / f"{stem}.inp"
            if _completed_gxtb(
                output,
                input_path,
                f"x23b_{phase}",
                campaign_identity,
                *markers,
            ):
                count += 1
        return count

    validation = payload.setdefault("validation", {})
    if not isinstance(validation, dict):
        raise ValueError("invalid validation provenance record")
    validation.update(
        {
            "gas_optimizations_expected": expected,
            "gas_optimizations_completed": count_local(
                "molecule_geoopt", "GXTB_mol_geoopt", ("PROGRAM ENDED", "GEOMETRY OPTIMIZATION COMPLETED")
            ),
            "gamma_cell_optimizations_expected": expected,
            "gamma_cell_optimizations_completed": count_local(
                "cellopt_gamma", "GXTB_gamma_cellopt", ("PROGRAM ENDED", "GEOMETRY OPTIMIZATION COMPLETED")
            ),
        }
    )

    cellopt_root_value = paths.get("k222_cellopt_root")
    cellopt_root = Path(str(cellopt_root_value)) if cellopt_root_value else None
    cellopt_records: dict[str, dict[str, str]] = {}
    if cellopt_root is not None:
        cellopt_manifest = cellopt_root.resolve().parent / "x23b_k222_cellopt_manifest.csv"
        if cellopt_manifest.is_file():
            with cellopt_manifest.open(newline="") as handle:
                for row in csv.DictReader(handle):
                    if row.get("method") == "GXTB" and row.get("system"):
                        cellopt_records[str(row["system"])] = row

    def cellopt_completed(system: str) -> bool:
        row = cellopt_records.get(system)
        if row is None:
            if source_policy != "gamma_cellopt_restart" or cellopt_root is None:
                return False
            legacy_run_dir = cellopt_root / "GXTB" / system / "k222_cellopt_keep_angles"
            return _completed_gxtb(
                legacy_run_dir / "cp2k.out",
                _project_input(
                    legacy_run_dir,
                    f"{system}_GXTB_k222_cellopt_keep_angles",
                ),
                "x23b_k222_cellopt",
                campaign_identity,
                "PROGRAM ENDED",
                "GEOMETRY OPTIMIZATION COMPLETED",
            )
        input_path = Path(str(row.get("input", "")))
        run_dir = Path(str(row.get("run_dir", "")))
        policy = str(row.get("source_policy", ""))
        schema = int(row.get("schema", "1") or "1")
        protocol_identity = None
        source_artifacts = None
        if schema >= 2:
            variant = str(row.get("variant", ""))
            lineage = str(row.get("lineage", ""))
            if not policy or not variant or not lineage:
                return False
            protocol_identity = {
                "source_policy": policy,
                "variant": variant,
                "lineage": lineage,
            }
            if policy == "experimental_reference":
                source_artifacts = {
                    "reference_input": Path(str(row.get("source_path", ""))),
                    "reference_structure": Path(str(row.get("structure_path", ""))),
                    "preflight_input": Path(str(row.get("preflight_input", ""))),
                    "preflight_output": Path(str(row.get("preflight_output", ""))),
                }
                if row.get("preflight_stamp"):
                    source_artifacts["preflight_stamp"] = Path(str(row["preflight_stamp"]))
            elif policy == "gamma_cellopt_restart":
                source_artifacts = {"gamma_restart": Path(str(row.get("source_restart", "")))}
            else:
                return False
        return _completed_gxtb(
            run_dir / "cp2k.out",
            input_path,
            "x23b_k222_cellopt",
            campaign_identity,
            "PROGRAM ENDED",
            "GEOMETRY OPTIMIZATION COMPLETED",
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )

    validation["k222_cell_optimizations_expected"] = expected
    validation["k222_cell_optimizations_completed"] = (
        sum(cellopt_completed(system) for system in systems)
        if cellopt_root is not None
        else 0
    )

    preflight_root_value = paths.get("experimental_k222_preflight_root")
    preflight_root = Path(str(preflight_root_value)) if preflight_root_value else None
    preflight_records: dict[str, dict[str, object]] = {}
    if preflight_root is not None:
        preflight_manifest = preflight_root / "experimental_k222_preflight_manifest.json"
        if preflight_manifest.is_file():
            try:
                preflight_payload = json.loads(preflight_manifest.read_text())
            except (json.JSONDecodeError, OSError):
                preflight_payload = {}
            if (
                preflight_payload.get("phase") == "x23b_experimental_k222_preflight"
                and preflight_payload.get("campaign_identity") == campaign_identity
            ):
                for row in preflight_payload.get("systems", []):
                    if isinstance(row, dict) and row.get("system"):
                        preflight_records[str(row["system"])] = row

    preflight_protocol = {
        "source_policy": "experimental_reference",
        "variant": "experimental_k222_preflight",
        "mesh": "MACDONALD 2 2 2 0.25 0.25 0.25",
        "symmetry": "SPGLIB reduced",
        "run_type": "ENERGY_FORCE",
        "stress": "ANALYTICAL GPa",
    }

    def preflight_completed(system: str) -> bool:
        row = preflight_records.get(system)
        if row is None:
            return False
        return _completed_gxtb(
            Path(str(row.get("output", ""))),
            Path(str(row.get("input", ""))),
            "x23b_experimental_k222_preflight",
            campaign_identity,
            "PROGRAM ENDED",
            "ENERGY| Total FORCE_EVAL",
            "FORCES| Atomic forces",
            "STRESS| Analytical stress tensor",
            protocol_identity=preflight_protocol,
            source_artifacts={
                "reference_input": Path(str(row.get("source_input", ""))),
                "reference_structure": Path(str(row.get("structure_path", ""))),
            },
        )

    validation["experimental_k222_preflight_expected"] = expected
    validation["experimental_k222_preflight_completed"] = (
        sum(preflight_completed(system) for system in systems)
        if preflight_root is not None
        else 0
    )

    final_roots = paths.get("final_single_point_roots", {})
    if not isinstance(final_roots, dict):
        final_roots = {}
    for mesh in (3, 4):
        mesh_id = f"k{mesh}{mesh}{mesh}"
        root_value = final_roots.get(mesh_id)
        root = Path(str(root_value)) if root_value else None
        final_records: dict[str, dict[str, str]] = {}
        if root is not None and (root / "manifest.csv").is_file():
            with (root / "manifest.csv").open(newline="") as handle:
                for row in csv.DictReader(handle):
                    if row.get("method") == "GXTB" and row.get("system"):
                        final_records[str(row["system"])] = row

        def final_completed(system: str) -> bool:
            if root is None:
                return False
            row = final_records.get(system)
            if row is None or not row.get("source_protocol_identity"):
                if source_policy != "gamma_cellopt_restart":
                    return False
                legacy_run = root / "GXTB" / system / f"{mesh_id}_sp_on_k222"
                return _completed_gxtb(
                    legacy_run / "cp2k.out",
                    _project_input(legacy_run, f"{system}_GXTB_{mesh_id}_sp_on_k222"),
                    f"x23b_final_{mesh_id}_on_k222",
                    campaign_identity,
                    "PROGRAM ENDED",
                    "ENERGY| Total FORCE_EVAL",
                )
            try:
                source_protocol = json.loads(row["source_protocol_identity"])
                source_payload = json.loads(row["source_artifacts"])
                if not isinstance(source_protocol, dict) or not isinstance(source_payload, dict):
                    return False
                if (
                    source_protocol.get("source_policy") != row.get("source_policy")
                    or source_protocol.get("variant") != row.get("source_variant")
                ):
                    return False
                target_protocol = {
                    "source_policy": row["source_policy"],
                    "source_variant": row["source_variant"],
                    "source_protocol_identity": source_protocol,
                    "target_mesh": mesh_id,
                }
                source_artifacts = {
                    name: Path(path) for name, path in source_payload.items()
                }
                source_artifacts.update(
                    {
                        "cellopt_input": Path(row["source_input"]),
                        "cellopt_output": Path(row["source_run_dir"]) / "cp2k.out",
                        "cellopt_restart": Path(row["source_restart"]),
                    }
                )
                input_path = Path(row["input"])
                run_dir = Path(row["run_dir"])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                return False
            return _completed_gxtb(
                run_dir / "cp2k.out",
                input_path,
                f"x23b_final_{mesh_id}_on_k222",
                campaign_identity,
                "PROGRAM ENDED",
                "ENERGY| Total FORCE_EVAL",
                protocol_identity=target_protocol,
                source_artifacts=source_artifacts,
            )

        validation[f"{mesh_id}_single_points_expected"] = expected
        validation[f"{mesh_id}_single_points_completed"] = (
            sum(final_completed(system) for system in systems)
            if root is not None
            else 0
        )

    for mesh_id in ("gamma", "k111", "k222", "k333"):
        validation[f"fixed_reference_{mesh_id}_single_points_expected"] = expected
        validation[f"fixed_reference_{mesh_id}_single_points_completed"] = sum(
            _completed_gxtb(
                benchmark_root
                / "runs"
                / "crystal_sp"
                / mesh_id
                / "GXTB"
                / f"{system}_GXTB_{mesh_id}_sp"
                / f"{system}_GXTB_{mesh_id}_sp.out",
                benchmark_root
                / "inputs"
                / "crystal_sp"
                / mesh_id
                / "GXTB"
                / f"{system}_GXTB_{mesh_id}_sp.inp",
                f"x23b_fixed_reference_{mesh_id}",
                campaign_identity,
                "PROGRAM ENDED",
                "ENERGY| Total FORCE_EVAL",
            )
            for system in systems
        )

    if publication_tables_updated is not None:
        validation["publication_tables_updated"] = publication_tables_updated
    source_gate = (
        "experimental_k222_preflight"
        if source_policy == "experimental_reference"
        else "gamma_cell_optimizations"
    )
    main_keys = (
        "gas_optimizations",
        source_gate,
        "k222_cell_optimizations",
        "k333_single_points",
        "k444_single_points",
    )
    complete_main = expected > 0 and all(
        int(validation.get(f"{key}_completed", 0)) == int(validation.get(f"{key}_expected", expected))
        for key in main_keys
    )
    any_completed = any(
        key.endswith("_completed") and isinstance(value, int) and value > 0
        for key, value in validation.items()
    )
    payload["status"] = (
        "production_complete"
        if complete_main
        else "production_in_progress"
        if any_completed
        else "runner_ready_production_not_started"
    )
    _atomic_json(path, payload)
    return path


def validate_method_input(text: str, method: str) -> None:
    """Reject inputs which could silently run a different model or mixer."""

    if method not in METHODS:
        raise ValueError(f"unsupported X23b method: {method}")
    tblite_methods = re.findall(r"^\s*METHOD\s+(GFN1|GFN2|GXTB)\s*$", text, flags=re.I | re.M)
    if not tblite_methods or tblite_methods[-1].upper() != method:
        raise ValueError(f"input does not select METHOD {method}")
    if method != "GXTB":
        return
    mixers = [value.upper() for value in re.findall(r"^\s*SCC_MIXER\s+(\S+)\s*$", text, flags=re.I | re.M)]
    if mixers != ["TBLITE"]:
        raise ValueError("GXTB production inputs require exactly one 'SCC_MIXER TBLITE'")
    if re.search(r"^\s*&TBLITE_MIXER\b", text, flags=re.I | re.M):
        raise ValueError("GXTB production inputs must not override the native save_tblite Fock-DIIS mixer")
    if re.search(r"^\s*&KPOINTS\b", text, flags=re.I | re.M):
        required = (
            re.search(r"^\s*SYMMETRY\s+T\s*$", text, flags=re.I | re.M),
            re.search(r"^\s*FULL_GRID\s+F\s*$", text, flags=re.I | re.M),
            re.search(r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$", text, flags=re.I | re.M),
            re.search(r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$", text, flags=re.I | re.M),
        )
        if not all(required):
            raise ValueError(
                "GXTB k-point production inputs require the SPGLIB-reduced mesh contract "
                "(SYMMETRY T, FULL_GRID F, SPGLIB backend and reduction method)"
            )


def _numbered_restart_step(path: Path) -> int:
    match = re.search(r"-1_(\d+)\.restart$", path.name)
    return int(match.group(1)) if match else -1


def final_restart(run_dir: Path) -> Path | None:
    """Return the restart needed to seed the next protocol phase."""

    unnumbered = list(run_dir.glob("*-1.restart"))
    if unnumbered:
        return max(unnumbered, key=lambda path: path.stat().st_mtime_ns)
    numbered = [path for path in run_dir.glob("*-1_*.restart") if _numbered_restart_step(path) >= 0]
    if numbered:
        return max(numbered, key=lambda path: (_numbered_restart_step(path), path.stat().st_mtime_ns))
    return None


def prune_gxtb_transients(run_dir: Path, *, keep_final_restart: bool = True) -> dict[str, object]:
    """Prune a *completed* GXTB run without touching reusable/audit files.

    Completion is deliberately checked by each phase-specific caller before
    entering here.  The path guard is an additional defence against applying
    the operation to the frozen GFN1/GFN2 production tree.
    """

    run_dir = run_dir.resolve()
    if "GXTB" not in run_dir.parts:
        raise ValueError(f"refusing to prune a run outside a GXTB directory: {run_dir}")

    retained_restart = final_restart(run_dir) if keep_final_restart else None
    candidates: set[Path] = set()
    patterns = (
        "*-RESTART.kp*",
        "*-RESTART.wfn.bak-*",
        "*-RESTART.wfn.bak_*",
        "*.restart.bak-*",
        "*.restart.bak_*",
        "*-pos-*.xyz.bak-*",
        "*-cell-*.cell.bak-*",
    )
    for pattern in patterns:
        candidates.update(path for path in run_dir.rglob(pattern) if path.is_file())
    for path in run_dir.rglob("*-1_*.restart"):
        if path.is_file() and _numbered_restart_step(path) >= 0:
            candidates.add(path)
    if not keep_final_restart:
        candidates.update(path for path in run_dir.rglob("*-1.restart") if path.is_file())
    if retained_restart is not None:
        candidates.discard(retained_restart)

    deleted: list[dict[str, object]] = []
    for path in sorted(candidates):
        try:
            relative = path.relative_to(run_dir).as_posix()
        except ValueError as exc:  # pragma: no cover - rglob cannot escape
            raise ValueError(f"transient escaped the GXTB run directory: {path}") from exc
        size = path.stat().st_size
        path.unlink()
        deleted.append({"path": relative, "bytes": size})

    report: dict[str, object] = {
        "policy": "completed GXTB jobs only",
        "deleted": deleted,
        "bytes_deleted": sum(int(row["bytes"]) for row in deleted),
        "retained_final_restart": (
            retained_restart.relative_to(run_dir).as_posix() if retained_restart is not None else None
        ),
        "retained_current_wfn": sorted(
            path.relative_to(run_dir).as_posix()
            for path in run_dir.rglob("*-RESTART.wfn")
            if path.is_file()
        ),
    }
    (run_dir / "prune_transients.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
