#!/usr/bin/env bash
set -euo pipefail

work=/home/kuehne88/work/gxtb-final-lowk-derivatives-20260719
binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
launcher=/home/kuehne88/work/gxtb-native-bvk-20260718/launch_pinned_cp2k.sh
expected_binary=b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f
expected_launcher=fd65cd0cff356ae9eaa8e556b88e1750c5193cc0c3a0b37f76983e7579d9405c
status_file=$work/controller.exit_status
log=$work/controller.log
minimum_available_gib=150
cpus=(60 61 62 63 64 65 66 67 68 69 70 71)

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

mkdir -p "$work/results"
rm -f "$status_file"
trap 'status=$?; printf "%s\n" "$status" > "$status_file"; trap - EXIT; exit "$status"' EXIT

available_gib() {
  awk '/^MemAvailable:/{printf "%d\n", $2/1024/1024}' /proc/meminfo
}

is_exact_complete() {
  local result=$1 input=$2 recorded_binary recorded_input actual_input
  [[ -f $result/exit_status ]] || return 1
  [[ $(tr -d '[:space:]' < "$result/exit_status") == 0 ]] || return 1
  grep -q 'PROGRAM ENDED AT' "$result/cp2k.out" 2>/dev/null || return 1
  recorded_binary=$(awk 'NR == 1 {print $1}' "$result/binary.sha256" 2>/dev/null || true)
  recorded_input=$(awk 'NR == 1 {print $1}' "$result/input.sha256" 2>/dev/null || true)
  actual_input=$(sha256sum "$input" | awk '{print $1}')
  [[ $recorded_binary == "$expected_binary" ]]
  [[ $recorded_input == "$actual_input" ]]
}

run_one() {
  local task=$1 cpu=$2 suite filename stem input result status
  suite=${task%%:*}
  filename=${task#*:}
  stem=${filename%.inp}
  input=$work/inputs/$suite/$filename
  result=$work/results/$suite/$stem
  [[ -r $input ]]
  if is_exact_complete "$result" "$input"; then
    printf '%s already_complete task=%s\n' "$(date --iso-8601=seconds)" "$task" >> "$log"
    return
  fi
  rm -rf "$result"
  mkdir -p "$(dirname "$result")"
  while (( $(available_gib) < minimum_available_gib )); do
    printf '%s waiting task=%s available_GiB=%s required_GiB=%s\n' \
      "$(date --iso-8601=seconds)" "$task" "$(available_gib)" \
      "$minimum_available_gib" >> "$log"
    sleep 30
  done
  while true; do
    printf '%s launching task=%s cpu=%s available_GiB=%s\n' \
      "$(date --iso-8601=seconds)" "$task" "$cpu" "$(available_gib)" >> "$log"
    set +e
    "$launcher" "lowk-${suite}-${stem}" "$cpu" "$binary" "$input" "$result" \
      > "$result-launcher.log" 2>&1
    status=$?
    set -e
    if [[ $status -eq 75 ]]; then
      sleep 30
      continue
    fi
    [[ $status -eq 0 ]]
    break
  done
  is_exact_complete "$result" "$input"
  printf '%s complete task=%s\n' "$(date --iso-8601=seconds)" "$task" >> "$log"
}

test "$(sha256sum "$binary" | awk '{print $1}')" = "$expected_binary"
test "$(sha256sum "$launcher" | awk '{print $1}')" = "$expected_launcher"
printf 'cp2k_sha256=%s\nlauncher_sha256=%s\ncpus=%s\nthreads=1\n' \
  "$expected_binary" "$expected_launcher" "$(IFS=,; echo "${cpus[*]}")" \
  > "$work/run_identity.txt"

for ((offset=0; offset<${#tasks[@]}; offset+=${#cpus[@]})); do
  pids=()
  for ((slot=0; slot<${#cpus[@]}; slot++)); do
    index=$((offset + slot))
    (( index < ${#tasks[@]} )) || break
    run_one "${tasks[$index]}" "${cpus[$slot]}" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
done

find "$work/results" -type f -name affinity_preexec.txt -print0 \
  | sort -z | xargs -0 sha256sum > "$work/affinity_proofs.sha256"
find "$work/results" -type f -name cp2k.out -print0 \
  | sort -z | xargs -0 sha256sum > "$work/outputs.sha256"
