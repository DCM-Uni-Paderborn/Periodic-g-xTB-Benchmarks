#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-native-bvk-20260718
binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
qualified_binary=b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f
cpu=137

log() {
  printf '%s %s\n' "$(date --iso-8601=seconds)" "$*"
}

validate_result() {
  local result_dir=$1
  [[ $(cat "$result_dir/exit_status") == 0 ]]
  grep -q 'PROGRAM ENDED AT' "$result_dir/cp2k.out"
  grep -q "$qualified_binary" "$result_dir/binary.sha256"
  [[ $(sha256sum "$binary" | awk '{print $1}') == "$qualified_binary" ]]
}

wait_for_memory() {
  local required_kib=398458880 available_kib
  while :; do
    available_kib=$(awk '/^MemAvailable:/{print $2}' /proc/meminfo)
    if (( available_kib >= required_kib )); then
      log "memory gate passed available_kib=$available_kib required_kib=$required_kib"
      return
    fi
    log "waiting for memory available_kib=$available_kib required_kib=$required_kib"
    sleep 30
  done
}

run_endpoint() {
  local job_id=$1 input=$2 result_dir=$3
  [[ -s $input ]]
  [[ ! -e $result_dir ]]
  wait_for_memory
  log "launching job_id=$job_id cpu=$cpu"
  bash "$root/launch_pinned_cp2k.sh" "$job_id" "$cpu" "$binary" "$input" "$result_dir"
  validate_result "$result_dir"
  log "qualified endpoint complete job_id=$job_id"
}

[[ $(sha256sum "$binary" | awk '{print $1}') == "$qualified_binary" ]]

run_endpoint dmc-k888-qualified-XIV \
  "$root/inputs/k888-reduced/XIV/input.inp" \
  "$root/runs/k888-reduced/XIV"

run_endpoint dmc-k777-qualified-XI \
  "$root/inputs/k777-reduced/XI/input.inp" \
  "$root/runs/k777-reduced/XI"

log 'qualified endpoint chain below 9x9x9 complete'
