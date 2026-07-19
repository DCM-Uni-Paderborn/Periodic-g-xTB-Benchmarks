#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: monitor_qualified_mixed_mae.sh ROOT REFERENCE_CSV PAPER_JSON REQUIRED_BINARY_SHA256 STATUS_DIR [INTERVAL_SECONDS]

Repeatedly evaluates the densest available, provenance-qualified same-mesh
DMC-ICE13 MAE. Set ONCE=1 for one evaluation without polling.
EOF
  exit 2
}

[[ $# -ge 5 && $# -le 6 ]] || usage

root=$1
reference=$2
paper=$3
required_binary=$(printf '%s' "$4" | tr '[:upper:]' '[:lower:]')
status_dir=$5
interval=${6:-60}
tool_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
tool=$tool_dir/dmc_mixed_mae.py
meshes=${MESHES:-9,8,7,6,5,4,3,2,1}

[[ $required_binary =~ ^[0-9a-f]{64}$ ]] || {
  printf 'invalid required executable SHA-256: %s\n' "$required_binary" >&2
  exit 2
}
[[ $interval =~ ^[1-9][0-9]*$ ]] || {
  printf 'invalid polling interval: %s\n' "$interval" >&2
  exit 2
}
[[ -x $tool ]] || {
  printf 'missing executable evaluator: %s\n' "$tool" >&2
  exit 2
}

latest=$status_dir/dmc_mixed_qualified_latest.tsv
history=$status_dir/dmc_mixed_qualified_history.tsv
signature_file=$status_dir/dmc_mixed_qualified_signature.sha256
readiness=$status_dir/dmc_mixed_qualified.status
failure_signature_file=$status_dir/dmc_mixed_qualified_failure.sha256
mkdir -p "$status_dir"

hash_stdin() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum
  else
    shasum -a 256
  fi
}

timestamp_now() {
  date -u '+%Y-%m-%dT%H:%M:%SZ'
}

phase_pattern='^(II|III|IV|VI|VII|VIII|IX|XI|XIII|XIV|XV|XVII)$'

while true; do
  temporary=$(mktemp "$status_dir/dmc_mixed_qualified.XXXXXX")
  stderr_file=$temporary.stderr
  if "$tool" "$root" "$reference" --meshes "$meshes" \
       --paper-json "$paper" --require-binary-sha256 "$required_binary" \
       >"$temporary" 2>"$stderr_file"; then
    signature=$(awk -F '\t' -v pattern="$phase_pattern" \
      '$1 ~ pattern {print $1 ":" $2}' "$temporary" | hash_stdin | awk '{print $1}')
    previous=$(awk 'NR == 1 {print $1}' "$signature_file" 2>/dev/null || true)
    if [[ $signature != "$previous" ]]; then
      timestamp=$(timestamp_now)
      current=$(awk -F '\t' '$1 == "mixed_mae_kj_mol" {print $2}' "$temporary")
      paper_value=$(awk -F '\t' '$1 == "paper_comparator_mae_kj_mol" {print $2}' "$temporary")
      paper_same_mesh=$(awk -F '\t' '$1 == "paper_comparator_all_same_mesh" {print $2}' "$temporary")
      improvement=$(awk -F '\t' '$1 == "mae_improvement_kj_mol" {print $2}' "$temporary")
      percent=$(awk -F '\t' '$1 == "mae_improvement_percent" {print $2}' "$temporary")
      mesh_vector=$(awk -F '\t' -v pattern="$phase_pattern" \
        '$1 ~ pattern {printf "%s%s:%s", separator, $1, $2; separator=","} END {print ""}' \
        "$temporary")
      mv "$temporary" "$latest"
      printf '%s  %s\n' "$signature" "$timestamp" >"$signature_file"
      if [[ ! -s $history ]]; then
        printf 'timestamp\tmixed_mae_kj_mol\tpaper_comparator_mae_kj_mol\tpaper_comparator_all_same_mesh\tmae_improvement_kj_mol\tmae_improvement_percent\tmesh_vector\n' \
          >"$history"
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$timestamp" "$current" "$paper_value" "$paper_same_mesh" \
        "$improvement" "$percent" "$mesh_vector" >>"$history"
      printf 'updated timestamp=%s mixed_mae=%s paper_comparator_mae=%s same_mesh=%s mesh_vector=%s\n' \
        "$timestamp" "$current" "$paper_value" "$paper_same_mesh" "$mesh_vector"
    else
      rm -f "$temporary"
    fi
    {
      printf 'status=READY\n'
      printf 'required_binary_sha256=%s\n' "$required_binary"
      printf 'updated=%s\n' "$(timestamp_now)"
    } >"$readiness"
    rm -f "$stderr_file"
    result=0
  else
    reason=$(tail -n 1 "$stderr_file")
    failure_signature=$(printf '%s' "$reason" | hash_stdin | awk '{print $1}')
    previous_failure=$(awk 'NR == 1 {print $1}' "$failure_signature_file" 2>/dev/null || true)
    if [[ $failure_signature != "$previous_failure" ]]; then
      printf 'not_ready timestamp=%s reason=%s\n' \
        "$(timestamp_now)" "$reason" >&2
      printf '%s\n' "$failure_signature" >"$failure_signature_file"
    fi
    {
      printf 'status=NOT_READY\n'
      printf 'required_binary_sha256=%s\n' "$required_binary"
      printf 'reason=%s\n' "$reason"
      printf 'checked=%s\n' "$(timestamp_now)"
    } >"$readiness"
    rm -f "$temporary" "$stderr_file"
    result=1
  fi
  if [[ ${ONCE:-0} == 1 ]]; then
    exit "$result"
  fi
  sleep "$interval"
done
