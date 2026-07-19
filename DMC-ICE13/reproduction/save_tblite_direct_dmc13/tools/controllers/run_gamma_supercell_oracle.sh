#!/usr/bin/env bash
set -euo pipefail

root=${DMC_ROOT:-/home/kuehne88/work/gxtb-native-bvk-20260718}
binary=${CP2K_BINARY:-/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp}
launcher=${PINNED_LAUNCHER:-$root/launch_pinned_cp2k.sh}
required_binary=${REQUIRED_BINARY_SHA256:-b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f}
input=${GAMMA_ORACLE_INPUT:-$root/validation/explicit-gamma-supercell-k222-XVII/input.inp}
required_input=${GAMMA_ORACLE_INPUT_SHA256:-bfd43e54957647f6ad6b8df8c13f2b9dbc34cde41038ea55b4ec5c63c8abbec1}
result="$root/validation/explicit-gamma-supercell-k222-XVII/result"
archive_root="$root/runs/pre-explicit-gamma-oracle-20260719"
status_file="$root/status/explicit-gamma-supercell-k222-XVII.status"
log="$root/status/explicit-gamma-supercell-k222-XVII.log"

mkdir -p "$archive_root" "$root/status" "${result%/*}"

recorded_digest() {
  awk 'NR == 1 {print $1}' "$1" 2>/dev/null || true
}

exact_complete() {
  local directory=$1
  [[ -f "$directory/exit_status" ]] || return 1
  [[ "$(tr -d '\n' <"$directory/exit_status")" == 0 ]] || return 1
  [[ "$(recorded_digest "$directory/binary.sha256")" == "$required_binary" ]] || return 1
  grep -q 'PROGRAM ENDED AT' "$directory/cp2k.out" 2>/dev/null
}

test "$(sha256sum "$binary" | awk '{print $1}')" = "$required_binary"
test "$(sha256sum "$input" | awk '{print $1}')" = "$required_input"

# The large 8^3 phase/reference pair has priority.  The oracle is forbidden
# from starting until both exact same-build outputs are present.
exact_complete "$root/runs/k888-reduced/VII"
exact_complete "$root/runs/k888-reduced/Ih"

if exact_complete "$result" && \
   [[ "$(recorded_digest "$result/input.sha256")" == "$required_input" ]]; then
  printf 'status=PASS\nrequired_binary_sha256=%s\ninput_sha256=%s\n' \
    "$required_binary" "$required_input" >"$status_file"
  exit 0
fi

if [[ -d "$result" ]]; then
  archived="$archive_root/result-$(date +%Y%m%dT%H%M%S)"
  mv "$result" "$archived"
  printf '%s archived prior result destination=%s\n' \
    "$(date --iso-8601=seconds)" "$archived" >>"$log"
fi
mkdir -p "$result"

while (( ${AVAILABLE_GIB_OVERRIDE:-$(awk '/^MemAvailable:/{printf "%d\n", $2/1024/1024}' /proc/meminfo)} < ${MINIMUM_AVAILABLE_GIB:-400} )); do
  printf '%s waiting for 400 GiB available memory\n' \
    "$(date --iso-8601=seconds)" >>"$log"
  sleep 30
done

printf '%s launching explicit Gamma-supercell oracle\n' \
  "$(date --iso-8601=seconds)" >>"$log"
"$launcher" dmc-k222-XVII-gamma-supercell 91 "$binary" "$input" "$result" \
  >"$result/launcher.log" 2>&1

exact_complete "$result"
[[ "$(recorded_digest "$result/input.sha256")" == "$required_input" ]]
printf 'status=PASS\nrequired_binary_sha256=%s\ninput_sha256=%s\n' \
  "$required_binary" "$required_input" >"$status_file"
printf '%s explicit Gamma-supercell oracle complete\n' \
  "$(date --iso-8601=seconds)" >>"$log"
