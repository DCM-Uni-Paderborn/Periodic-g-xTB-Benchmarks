#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  printf 'usage: %s CPU WORK_ROOT\n' "$0" >&2
  exit 64
fi

cpu=$1
work=$2
binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
launcher=/home/kuehne88/work/gxtb-native-bvk-20260718/launch_pinned_cp2k.sh
expected_binary=b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f
expected_launcher=ee99ad6085ba5dd78b5fadd2f17bf630d3838eaae18da6a818684ad780fedeb6
minimum_available_gib=228

tasks=(
  derivative:CH4_gxtb_gamma_force_stress.inp
  derivative:CH4_gxtb_kp_111_force_stress.inp
  derivative:CH4_gxtb_kp_full_222_force_stress.inp
  derivative:CH4_gxtb_kp_k290_force_stress.inp
  derivative:CH4_gxtb_kp_spglib_111_force_stress.inp
  derivative:CH4_gxtb_kp_spglib_222_force_stress.inp
  derivative:H2O_gxtb_gamma_force_stress.inp
  derivative:H2O_gxtb_kp_gamma_force_stress.inp
  derivative:H2_gxtb_gamma_supercell_311_force_stress.inp
  derivative:H2_gxtb_kp_311_force_stress.inp
  derivative:H2_gxtb_kp_311_tr_force_stress.inp
  partial:gxtb_1d_native_gamma_centered_k211.inp
  partial:gxtb_1d_native_k211.inp
  partial:gxtb_1d_supercell_k111.inp
  partial:gxtb_1d_x_k211_force_stress.inp
  partial:gxtb_1d_x_k211_force_stress_full.inp
  partial:gxtb_1d_x_k211_force_stress_spglib.inp
  partial:gxtb_2d_native_gamma_centered_k212.inp
  partial:gxtb_2d_native_k212.inp
  partial:gxtb_2d_supercell_k111.inp
  partial:gxtb_2d_xz_k212_force_stress.inp
  partial:gxtb_2d_xz_k212_force_stress_full.inp
  partial:gxtb_2d_xz_k212_force_stress_spglib.inp
)

if [[ ! $cpu =~ ^[0-9]+$ ]]; then
  printf 'invalid CPU: %s\n' "$cpu" >&2
  exit 64
fi
test "$(sha256sum "$binary" | awk '{print $1}')" = "$expected_binary"
test "$(sha256sum "$launcher" | awk '{print $1}')" = "$expected_launcher"
test -d "$work/inputs/derivative"
test -d "$work/inputs/partial"
test ! -e "$work/results"
mkdir -p "$work/results"

{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'cp2k_sha256=%s\n' "$expected_binary"
  printf 'launcher_sha256=%s\n' "$expected_launcher"
  printf 'cpu=%s\nthreads=1\n' "$cpu"
} > "$work/run_identity.txt"

for task in "${tasks[@]}"; do
  suite=${task%%:*}
  filename=${task#*:}
  stem=${filename%.inp}
  input=$work/inputs/$suite/$filename
  result=$work/results/$suite/$stem
  test -r "$input"
  test ! -e "$result"

  while true; do
    available_gib=$(awk '/^MemAvailable:/{printf "%d\n", $2/1024/1024}' /proc/meminfo)
    live_cp2k=$(ps -eo pid=,stat=,rss=,args= | awk '$2 !~ /^Z/ && /[c]p2k.psmp/ {print}')
    if [[ -z $live_cp2k && $available_gib -ge $minimum_available_gib ]]; then
      break
    fi
    {
      printf '%s waiting task=%s available_GiB=%s required_GiB=%s\n' \
        "$(date --iso-8601=seconds)" "$task" "$available_gib" "$minimum_available_gib"
      [[ -z $live_cp2k ]] || printf '%s\n' "$live_cp2k"
    } >> "$work/controller.log"
    sleep 15
  done

  mkdir -p "$(dirname "$result")"
  printf '%s launching task=%s cpu=%s available_GiB=%s\n' \
    "$(date --iso-8601=seconds)" "$task" "$cpu" "$available_gib" \
    >> "$work/controller.log"
  "$launcher" "lowk-rerun-${suite}-${stem}" "$cpu" "$binary" "$input" "$result" \
    > "$result-launcher.log" 2>&1
  test "$(tr -d '[:space:]' < "$result/exit_status")" = 0
  grep -q 'PROGRAM ENDED AT' "$result/cp2k.out"
  test "$(awk 'NR == 1 {print $1}' "$result/binary.sha256")" = "$expected_binary"
  test "$(awk 'NR == 1 {print $1}' "$result/input.sha256")" = \
    "$(sha256sum "$input" | awk '{print $1}')"
  printf '%s complete task=%s\n' "$(date --iso-8601=seconds)" "$task" \
    >> "$work/controller.log"
done

printf 'completed_utc=%s\nstatus=PASS\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  > "$work/STATUS"
