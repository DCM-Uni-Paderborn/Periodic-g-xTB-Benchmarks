#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-native-bvk-20260718
campaign="$root/validation/xvii-full-derivatives-current"
binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
qualified_binary=b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f
cpu=153

log() {
  printf '%s %s\n' "$(date --iso-8601=seconds)" "$*"
}

live_cp2k_pids() {
  local pid state comm uid self_uid
  self_uid=$(id -u)
  for status in /proc/[0-9]*/status; do
    [[ -r $status ]] || continue
    pid=${status#/proc/}
    pid=${pid%/status}
    state=$(awk '/^State:/{print $2}' "$status")
    [[ $state == Z ]] && continue
    uid=$(awk '/^Uid:/{print $2}' "$status")
    [[ $uid == "$self_uid" ]] || continue
    comm=$(awk '/^Name:/{print $2}' "$status")
    [[ $comm == cp2k.psmp || $comm == cp2k.popt || $comm == cp2k.ssmp || $comm == cp2k.sopt ]] || continue
    printf '%s\n' "$pid"
  done
}

wait_for_dense_chain() {
  while pgrep -u "$USER" -f "$root/tools/continue_dmc_below_9.sh" >/dev/null; do
    sleep 30
  done

  # Never replace a failed or interrupted DMC endpoint chain with derivative
  # work.  Both remaining sub-9x9x9 prerequisites must be complete with the
  # frozen executable before this campaign is allowed to consume the node.
  validate_result "$root/runs/k888-reduced/XIV"
  validate_result "$root/runs/k777-reduced/XI"

  while [[ -n $(live_cp2k_pids) ]]; do
    log "waiting for existing CP2K process ids=$(live_cp2k_pids | tr '\n' ',')"
    sleep 30
  done
}

wait_for_memory() {
  local required_kib=100663296 available_kib
  while :; do
    available_kib=$(awk '/^MemAvailable:/{print $2}' /proc/meminfo)
    if (( available_kib >= required_kib )); then
      return
    fi
    log "waiting for memory available_kib=$available_kib required_kib=$required_kib"
    sleep 30
  done
}

validate_result() {
  local result_dir=$1
  [[ $(cat "$result_dir/exit_status") == 0 ]]
  grep -q 'PROGRAM ENDED AT' "$result_dir/cp2k.out"
  grep -q "$qualified_binary" "$result_dir/binary.sha256"
  [[ $(sha256sum "$binary" | awk '{print $1}') == "$qualified_binary" ]]
}

run_case() {
  local name=$1 input="$campaign/inputs/$1/input.inp" result="$campaign/runs/$1"
  [[ -s $input ]]
  if [[ -e $result ]]; then
    validate_result "$result"
    log "qualified $name already present"
    return
  fi
  wait_for_memory
  log "launching $name on singleton CPU $cpu"
  bash "$root/launch_pinned_cp2k.sh" "gxtb-xvii-$name" "$cpu" "$binary" "$input" "$result"
  validate_result "$result"
  log "qualified $name complete"
}

[[ $(sha256sum "$binary" | awk '{print $1}') == "$qualified_binary" ]]
wait_for_dense_chain
run_case full
run_case reduced
run_case force-plus
run_case force-minus
run_case strain-plus
run_case strain-minus
log 'current-binary ice-XVII derivative requalification complete'
