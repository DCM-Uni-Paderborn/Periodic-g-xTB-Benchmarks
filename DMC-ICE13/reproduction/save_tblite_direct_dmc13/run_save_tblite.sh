#!/usr/bin/env bash
set -euo pipefail

: "${TBLITE_EXE:?Set TBLITE_EXE to the save_tblite CLI executable}"

root=$(cd "$(dirname "$0")" && pwd)
meshes=${MESHES:-"1 2 3 4"}
accuracy=${ACCURACY:-0.1}
iterations=${ITERATIONS:-300}
result_root=${RESULT_ROOT:-"$root/results/recalculated"}
read -r -a phases <<< "${PHASES:-Ih II III IV VI VII VIII IX XI XIII XIV XV XVII}"

mkdir -p "$result_root"
sha256_files() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$@"
  else
    shasum -a 256 "$@"
  fi
}
executable_hash=$(sha256_files "$TBLITE_EXE" | awk '{print $1}')
{
  printf 'executable=%s\n' "$TBLITE_EXE"
  printf 'executable_sha256=%s\n' "$executable_hash"
  printf 'meshes=%s\n' "$meshes"
  printf 'phases=%s\n' "${phases[*]}"
  printf 'accuracy=%s\n' "$accuracy"
  printf 'iterations=%s\n' "$iterations"
  printf 'gradients=%s\n' "${GRADIENTS:-0}"
  "$TBLITE_EXE" --version 2>&1 || true
} > "$result_root/run_metadata.txt"

qualified_existing_result() {
  local output=$1
  local input=$2
  [[ -s "$output/result.json" ]] || return 1
  [[ -s "$output/process.out" ]] || return 1
  [[ -f "$output/exit_status" ]] || return 1
  [[ $(tr -d '[:space:]' < "$output/exit_status") == 0 ]] || return 1
  [[ -f "$output/binary.sha256" ]] || return 1
  [[ -f "$output/input.sha256" ]] || return 1
  [[ $(awk 'NR == 1 {print $1}' "$output/binary.sha256") == "$executable_hash" ]] || return 1
  [[ $(awk 'NR == 1 {print $1}' "$output/input.sha256") == $(sha256_files "$input" | awk '{print $1}') ]] || return 1
  grep -q 'total energy' "$output/process.out" || return 1
  grep -q 'JSON dump of results written' "$output/process.out" || return 1
}

for mesh in $meshes; do
  mesh_id="k${mesh}${mesh}${mesh}"
  for phase in "${phases[@]}"; do
    input="$root/structures/$mesh_id/$phase/POSCAR"
    output="$result_root/$mesh_id/$phase"
    mkdir -p "$output"
    if [[ ${SKIP_EXISTING:-0} == 1 ]] && qualified_existing_result "$output" "$input"; then
      continue
    fi
    sha256_files "$TBLITE_EXE" > "$output/binary.sha256"
    sha256_files "$input" > "$output/input.sha256"
    rm -f "$output/result.json" "$output/gradient.txt" "$output/exit_status"
    set +e
    if [[ ${GRADIENTS:-0} == 1 ]]; then
      "$TBLITE_EXE" run --method gxtb --acc "$accuracy" --iterations "$iterations" \
        --no-restart --json "$output/result.json" --grad "$output/gradient.txt" "$input" \
        > "$output/process.out" 2>&1
    else
      "$TBLITE_EXE" run --method gxtb --acc "$accuracy" --iterations "$iterations" \
        --no-restart --json "$output/result.json" "$input" \
        > "$output/process.out" 2>&1
    fi
    status=$?
    set -e
    if [[ $status -eq 0 && ! -s "$output/result.json" ]]; then
      status=66
    fi
    printf '%s\n' "$status" > "$output/exit_status"
    if [[ $status -ne 0 ]]; then
      printf 'save_tblite CLI failed for %s/%s with exit status %s\n' \
        "$mesh_id" "$phase" "$status" >&2
      exit "$status"
    fi
    sha256_files \
      "$input" "$TBLITE_EXE" "$output/result.json" "$output/process.out" \
      "$output/exit_status" > "$output/SHA256SUMS"
  done
done
