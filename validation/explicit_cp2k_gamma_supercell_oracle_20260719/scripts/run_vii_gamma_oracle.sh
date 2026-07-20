#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-native-bvk-20260718
campaign="$root/validation/vii-gamma-supercell-oracle"
binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
binary_sha256=b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f
input_sha256=6481b70e1a437e92247649e06ff90798358f0606dedf60c07b11e4ecab557a21
minimum_available_kib=134217728
cpu=169

input="$campaign/input.inp"
result="$campaign/run"
[[ $(sha256sum "$binary" | awk '{print $1}') == "$binary_sha256" ]]
[[ $(sha256sum "$input" | awk '{print $1}') == "$input_sha256" ]]
[[ ! -e $result ]]

while :; do
  available_kib=$(awk '/^MemAvailable:/{print $2}' /proc/meminfo)
  if (( available_kib >= minimum_available_kib )); then
    break
  fi
  sleep 30
done

bash "$root/launch_pinned_cp2k.sh" \
  dmc-vii-k222-gamma-oracle \
  "$cpu" \
  "$binary" \
  "$input" \
  "$result"

[[ $(cat "$result/exit_status") == 0 ]]
grep -q 'PROGRAM ENDED AT' "$result/cp2k.out"
grep -q "$binary_sha256" "$result/binary.sha256"
