#!/usr/bin/env bash
set -euo pipefail

root=${DMC_ROOT:-/home/kuehne88/work/gxtb-native-bvk-20260718}
binary=${CP2K_BINARY:-/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp}
launcher=${PINNED_LAUNCHER:-$root/launch_pinned_cp2k.sh}
gamma_oracle_controller=${GAMMA_ORACLE_CONTROLLER:-$root/run_gamma_supercell_oracle.sh}
expected_binary=${REQUIRED_BINARY_SHA256:-b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f}
archive_root="$root/runs/pre-strict-adaptive-20260719"
decision_root="$root/status/strict-adaptive-completion"
log="$root/status/strict-adaptive-completion.log"
status_file="$root/status/strict-adaptive-completion.status"
minimum_available_gib=${MINIMUM_AVAILABLE_GIB:-400}
convergence_threshold=${CONVERGENCE_THRESHOLD:-0.10}
maximum_mesh=${MAXIMUM_MESH:-12}

mkdir -p "$archive_root" "$decision_root" "$root/status"

available_gib() {
  if [[ -n ${AVAILABLE_GIB_OVERRIDE:-} ]]; then
    printf '%s\n' "$AVAILABLE_GIB_OVERRIDE"
    return
  fi
  awk '/^MemAvailable:/{printf "%d\n", $2/1024/1024}' /proc/meminfo
}

mesh_dir() {
  local mesh=$1
  printf 'k%s%s%s-reduced\n' "$mesh" "$mesh" "$mesh"
}

wait_for_memory() {
  while (( $(available_gib) < minimum_available_gib )); do
    printf '%s waiting available_GiB=%s required_GiB=%s\n' \
      "$(date --iso-8601=seconds)" "$(available_gib)" \
      "$minimum_available_gib" >>"$log"
    sleep 30
  done
}

wait_for_preceding_controllers() {
  local pattern
  pattern='run_(ih8_priority|k888_final_pair|dmc8_remaining_when_safe|response_corrected_endpoints|vii_k222_symmetry_gate)\.sh|validate_response_corrected_endpoints\.sh'
  while pgrep -u "$USER" -f "$pattern" >/dev/null; do
    printf '%s waiting for preceding controllers\n' \
      "$(date --iso-8601=seconds)" >>"$log"
    sleep 30
  done
}

is_exact_complete() {
  local result=$1 digest
  [[ -f "$result/exit_status" ]] || return 1
  [[ "$(tr -d '\n' <"$result/exit_status")" == 0 ]] || return 1
  grep -q 'PROGRAM ENDED AT' "$result/cp2k.out" 2>/dev/null || return 1
  digest=$(awk 'NR == 1 {print $1}' "$result/binary.sha256" 2>/dev/null || true)
  [[ "$digest" == "$expected_binary" ]]
}

ensure_input() {
  local mesh=$1 phase=$2 target source_mesh source provenance
  target="$root/inputs/$(mesh_dir "$mesh")/$phase/input.inp"
  [[ -r "$target" ]] && return
  source_mesh=$((mesh - 1))
  while (( source_mesh > 0 )); do
    source="$root/inputs/$(mesh_dir "$source_mesh")/$phase/input.inp"
    [[ -r "$source" ]] && break
    source_mesh=$((source_mesh - 1))
  done
  [[ -r "$source" ]]
  provenance="${target%/*}/mesh_rewrite.json"
  python3 "$root/tools/build_native_mesh_input.py" \
    "$source" "$target" "$mesh" --provenance "$provenance"
  printf '%s generated input mesh=%s phase=%s source_mesh=%s\n' \
    "$(date --iso-8601=seconds)" "$mesh" "$phase" "$source_mesh" >>"$log"
}

run_one() {
  local mesh=$1 phase=$2 cpu=$3 md input result archive
  md=$(mesh_dir "$mesh")
  input="$root/inputs/$md/$phase/input.inp"
  result="$root/runs/$md/$phase"
  ensure_input "$mesh" "$phase"
  if is_exact_complete "$result"; then
    printf '%s already complete mesh=%s phase=%s\n' \
      "$(date --iso-8601=seconds)" "$mesh" "$phase" >>"$log"
    return
  fi
  if [[ -d "$result" ]]; then
    archive="$archive_root/$md/$phase-$(date +%Y%m%dT%H%M%S)"
    mkdir -p "${archive%/*}"
    mv "$result" "$archive"
    printf '%s archived prior result mesh=%s phase=%s destination=%s\n' \
      "$(date --iso-8601=seconds)" "$mesh" "$phase" "$archive" >>"$log"
  fi
  mkdir -p "$result"
  printf '%s launching mesh=%s phase=%s cpu=%s available_GiB=%s\n' \
    "$(date --iso-8601=seconds)" "$mesh" "$phase" "$cpu" \
    "$(available_gib)" >>"$log"
  "$launcher" "dmc-k${mesh}${mesh}${mesh}-strict-$phase" \
    "$cpu" "$binary" "$input" "$result" >"$result/launcher.log" 2>&1
  is_exact_complete "$result"
  printf '%s complete mesh=%s phase=%s\n' \
    "$(date --iso-8601=seconds)" "$mesh" "$phase" >>"$log"
}

run_batch() {
  local specification mesh phase cpu
  local -a pids=()
  wait_for_memory
  printf '%s batch start specifications=%s available_GiB=%s\n' \
    "$(date --iso-8601=seconds)" "$*" "$(available_gib)" >>"$log"
  for specification in "$@"; do
    IFS=: read -r mesh phase cpu <<<"$specification"
    run_one "$mesh" "$phase" "$cpu" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
  printf '%s batch complete specifications=%s\n' \
    "$(date --iso-8601=seconds)" "$*" >>"$log"
}

run_phase_round() {
  local mesh=$1
  shift
  local -a cpus
  if (( mesh >= 8 )); then
    cpus=(83)
  elif (( mesh == 7 )); then
    cpus=(83 79)
  elif (( mesh == 6 )); then
    # Archived peak-memory data bound a six-job batch below the 400 GiB
    # reserve while keeping every calculation on its own physical CPU.
    cpus=(83 81 79 77 75 73)
  else
    # The complete 4^3 and 5^3 phase sets peak at about 72 and 249 GiB,
    # respectively, so all twelve non-reference phases can run together.
    cpus=(83 81 79 77 75 73 71 69 67 65 63 61)
  fi
  local -a batch=()
  local phase cpu index=0
  for phase in "$@"; do
    cpu=${cpus[$((index % ${#cpus[@]}))]}
    batch+=("$mesh:$phase:$cpu")
    index=$((index + 1))
    if (( ${#batch[@]} == ${#cpus[@]} )); then
      run_batch "${batch[@]}"
      batch=()
    fi
  done
  if (( ${#batch[@]} )); then
    run_batch "${batch[@]}"
  fi
}

run_reference_and_phase_round() {
  local mesh=$1
  shift
  local -a cpus=(85 83 81 79 77 75 73 71 69 67 65 63 61)
  local -a specifications=("$mesh:Ih:${cpus[0]}")
  local phase index=1
  for phase in "$@"; do
    specifications+=("$mesh:$phase:${cpus[$index]}")
    index=$((index + 1))
  done
  run_batch "${specifications[@]}"
}

pair_status() {
  local phase=$1 previous=$2 current=$3 output=$4
  set +e
  python3 "$root/tools/dmc_phase_convergence.py" \
    "$root" "$previous" "$current" "$phase" \
    --threshold "$convergence_threshold" \
    --require-binary-sha256 "$expected_binary" \
    >"$output" 2>&1
  local rc=$?
  set -e
  return "$rc"
}

test "$(sha256sum "$binary" | awk '{print $1}')" = "$expected_binary"
printf '%s strict adaptive completion waiting\n' \
  "$(date --iso-8601=seconds)" >>"$log"
wait_for_preceding_controllers

# Preserve the memory-safe priority order requested for the 8^3 pair.  A
# completed, hash-qualified VII result is retained; otherwise it runs alone.
# Only after that result exists may the equally qualified Ih reference run.
run_one 8 VII 91
run_one 8 Ih 85
test -x "$gamma_oracle_controller"
"$gamma_oracle_controller"

phases=(II III IV VI VII VIII IX XI XIII XIV XV XVII)

# Qualify the complete 4^3 boundary explicitly.  The archived response-
# corrected 1^3--4^3 series proves that no earlier adjacent pair passes.  The
# current controller still rechecks every 4^3 result against the exact binary
# hash before it can participate in the adaptive sequence.
run_one 4 Ih 85
run_phase_round 4 "${phases[@]}"

# Every phase needs 5^3 once because none passed through 3^3->4^3.  The
# archived peak-memory envelope for Ih plus all twelve phase calculations is
# about 266 GiB, so the complete 5^3 set can use thirteen disjoint physical
# CPUs concurrently while retaining the 400-GiB pre-launch availability gate.
# Denser rounds are generated exclusively from the unresolved set of the
# preceding adjacent pair.
run_reference_and_phase_round 5 "${phases[@]}"

incomplete=0
pending=()
for phase in "${phases[@]}"; do
  phase_dir="$decision_root/$phase"
  mkdir -p "$phase_dir"
  if pair_status "$phase" 4 5 "$phase_dir/k4-k5.tsv"; then
    continue
  else
    rc=$?
  fi
  if (( rc == 1 )); then
    pending+=("$phase")
  else
    incomplete=$((incomplete + 1))
  fi
done

previous=5
for ((mesh = 6; mesh <= maximum_mesh; mesh++)); do
  (( ${#pending[@]} )) || break
  run_one "$mesh" Ih 85
  run_phase_round "$mesh" "${pending[@]}"
  next_pending=()
  for phase in "${pending[@]}"; do
    phase_dir="$decision_root/$phase"
    if pair_status "$phase" "$previous" "$mesh" \
      "$phase_dir/k${previous}-k${mesh}.tsv"; then
      continue
    else
      rc=$?
    fi
    if (( rc == 1 )); then
      next_pending+=("$phase")
    else
      incomplete=$((incomplete + 1))
    fi
  done
  if (( ${#next_pending[@]} )); then
    pending=("${next_pending[@]}")
  else
    pending=()
  fi
  previous=$mesh
done
unresolved=${#pending[@]}

set +e
python3 "$root/tools/select_adaptive_endpoints.py" \
  "$root" "$root/tools/dmc_ice13_relative_energies.csv" \
  --meshes "$(seq -s, 4 "$maximum_mesh")" \
  --threshold "$convergence_threshold" \
  --require-binary-sha256 "$expected_binary" \
  --output-json "$decision_root/adaptive_endpoints.json" \
  --output-csv "$decision_root/adaptive_endpoints.csv" \
  >"$decision_root/adaptive_endpoints.stdout.json" \
  2>"$decision_root/adaptive_endpoints.stderr"
selector_rc=$?
set -e

{
  printf 'required_binary_sha256=%s\n' "$expected_binary"
  printf 'threshold_kj_mol_per_water=%s\n' "$convergence_threshold"
  printf 'maximum_mesh=%s\n' "$maximum_mesh"
  printf 'unresolved_phases=%s\n' "$unresolved"
  printf 'incomplete_checks=%s\n' "$incomplete"
  printf 'selector_rc=%s\n' "$selector_rc"
  if (( incomplete != 0 || selector_rc == 2 )); then
    printf 'status=INCOMPLETE\n'
  elif (( unresolved != 0 || selector_rc == 1 )); then
    printf 'status=UNRESOLVED\n'
  else
    printf 'status=PASS\n'
  fi
} >"$status_file"

printf '%s strict adaptive completion ended unresolved=%s incomplete=%s selector_rc=%s\n' \
  "$(date --iso-8601=seconds)" "$unresolved" "$incomplete" "$selector_rc" >>"$log"
if (( incomplete != 0 || selector_rc == 2 )); then
  exit 2
fi
if (( unresolved != 0 || selector_rc == 1 )); then
  exit 1
fi
