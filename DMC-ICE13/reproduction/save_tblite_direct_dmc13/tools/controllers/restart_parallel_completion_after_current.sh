#!/usr/bin/env bash
set -euo pipefail

root=${DMC_ROOT:-/home/kuehne88/work/gxtb-native-bvk-20260718}
binary=${CP2K_BINARY:-/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp}
required_binary=${REQUIRED_BINARY_SHA256:-b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f}
old_controller_pid=${OLD_CONTROLLER_PID:?OLD_CONTROLLER_PID is required}
current_result=${CURRENT_RESULT:-$root/runs/k444-reduced/Ih}
log="$root/status/parallel-controller-restart.log"

recorded_digest() {
  awk 'NR == 1 {print $1}' "$1" 2>/dev/null || true
}

exact_complete() {
  [[ -f "$current_result/exit_status" ]] || return 1
  [[ "$(tr -d '\n' <"$current_result/exit_status")" == 0 ]] || return 1
  [[ "$(recorded_digest "$current_result/binary.sha256")" == "$required_binary" ]] || return 1
  grep -q 'PROGRAM ENDED AT' "$current_result/cp2k.out" 2>/dev/null
}

test "$(sha256sum "$binary" | awk '{print $1}')" = "$required_binary"
printf '%s waiting for current calculation pid=%s result=%s\n' \
  "$(date --iso-8601=seconds)" "$old_controller_pid" "$current_result" >>"$log"

while ! exact_complete; do
  if [[ -f "$current_result/exit_status" ]] && \
     [[ "$(tr -d '\n' <"$current_result/exit_status")" != 0 ]]; then
    printf '%s current calculation failed; resuming original controller\n' \
      "$(date --iso-8601=seconds)" >>"$log"
    kill -CONT "$old_controller_pid" 2>/dev/null || true
    exit 2
  fi
  sleep 5
done

# Do not replace the scheduler until every live CP2K or MPI process owned by
# this account has left.  Zombie entries are non-live and deliberately ignored.
while true; do
  live=0
  for proc in /proc/[0-9]*; do
    [[ $(stat -c %U "$proc" 2>/dev/null || true) == "$USER" ]] || continue
    comm=$(cat "$proc/comm" 2>/dev/null || true)
    case "$comm" in
      cp2k.psmp|cp2k.popt|mpirun|mpiexec|orterun) ;;
      *) continue ;;
    esac
    state=$(awk '/^State:/{print $2}' "$proc/status" 2>/dev/null || true)
    [[ $state == Z ]] || live=1
  done
  (( live == 0 )) && break
  sleep 2
done

if [[ -d /proc/$old_controller_pid ]]; then
  kill -KILL "$old_controller_pid"
fi

out="$root/status/strict-adaptive-completion.controller.out"
nohup setsid env \
  DMC_ROOT="$root" \
  CP2K_BINARY="$binary" \
  REQUIRED_BINARY_SHA256="$required_binary" \
  CONVERGENCE_THRESHOLD=0.10 \
  MINIMUM_AVAILABLE_GIB=400 \
  bash "$root/run_strict_adaptive_completion.sh" \
  >"$out" 2>&1 </dev/null &
new_pid=$!
printf '%s\n' "$new_pid" >"$root/status/strict-adaptive-completion.controller.pid"
printf '%s launched parallel completion controller pid=%s\n' \
  "$(date --iso-8601=seconds)" "$new_pid" >>"$log"
