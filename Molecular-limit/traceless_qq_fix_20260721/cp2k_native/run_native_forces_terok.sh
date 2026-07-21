#!/usr/bin/env bash
set -euo pipefail

cpu=${CPU:-42}
build_root=/home/kuehne88/work/gxtb-multipole-limit-fix-20260721
work=${WORK:-/home/kuehne88/work/h2o-molecular-limit-native-20260721}
binary=$build_root/cp2k-relink/cp2k.psmp
provider=$build_root/save_tblite-install
launcher=/home/kuehne88/work/gxtb-native-bvk-20260718/launch_pinned_cp2k.sh
template=$work/input_templates/H2O_gxtb_molecular_limit_native_forces.inp
raw=$work/raw
status=$work/status
source_commit=fad7fe4b188f99794d7c047d5b710667c3a2ce84
boxes=(8 10 12 15 20 30 40 50 60 80 100 150 200 250)

test "$(cat "$build_root/build.status")" = complete
test -x "$binary"
expected_binary=$(sha256sum "$binary" | awk '{print $1}')
expected_libcp2k=$(sha256sum "$build_root/cp2k-relink/libcp2k.so.2026.2" | awk '{print $1}')
expected_libtblite=$(sha256sum "$provider/lib/libtblite.a" | awk '{print $1}')
mkdir -p "$raw" "$status"

export PATH="$provider/bin:$PATH"
export LD_LIBRARY_PATH="$build_root/cp2k-relink:$provider/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

snapshot_resources() {
  local case_name=$1 snapshot=$status/${case_name}.prelaunch.txt
  {
    printf 'utc='; date -u +%Y-%m-%dT%H:%M:%SZ
    awk '/^(MemTotal|MemFree|MemAvailable|SwapTotal|SwapFree):/{print}' /proc/meminfo
    printf 'candidate_peak_budget_kib=%s\n' $((2 * 1024 * 1024))
    printf 'minimum_reserve_kib=%s\n' $((128 * 1024 * 1024))
    ps -eo pid=,ppid=,uid=,rss=,stat=,psr=,comm=,args= --sort=-rss
    printf 'reservations:\n'
    find /tmp -maxdepth 2 -type f -path '/tmp/gxtb-cpu-reservations-*/*.reservation' \
      -print -exec sed -n '1p' {} \; 2>/dev/null || true
    printf 'live_cp2k_affinities:\n'
    for proc in /proc/[0-9]*; do
      pid=${proc##*/}
      comm=$(cat "$proc/comm" 2>/dev/null || true)
      case "$comm" in
        cp2k.psmp|cp2k.popt|mpirun|mpiexec|orterun) ;;
        *) continue ;;
      esac
      allowed=$(awk '/^Cpus_allowed_list:/{print $2}' "$proc/status" 2>/dev/null || true)
      printf '%s %s %s\n' "$pid" "$comm" "$allowed"
    done
  } > "$snapshot"
  available=$(awk '/^MemAvailable:/{print $2}' /proc/meminfo)
  (( available - 2 * 1024 * 1024 >= 128 * 1024 * 1024 ))
}

render_input() {
  local case_name=$1 periodicity=$2 box=$3 destination=$4
  sed \
    -e "s/\${PROJECT}/$case_name/g" \
    -e "s/\${L}/$box/g" \
    -e "s/\${PBC}/$periodicity/g" \
    "$template" > "$destination"
}

run_case() {
  local case_name=$1 periodicity=$2 box=$3 directory
  directory=$raw/$case_name
  mkdir -p "$directory"
  render_input "$case_name" "$periodicity" "$box" "$directory/input.inp"
  snapshot_resources "$case_name"
  bash "$launcher" "h2o-native-$case_name" "$cpu" "$binary" \
    "$directory/input.inp" "$directory"
  test "$(cat "$directory/exit_status")" = 0
  grep -q 'PROGRAM ENDED AT' "$directory/cp2k.out"
  grep -q 'FORCES| Atomic forces \[hartree/bohr\]' "$directory/cp2k.out"
  test "$(awk '{print $1}' "$directory/binary.sha256")" = "$expected_binary"
  test "$(awk -F= '/allowed=/{for(i=1;i<=NF;i++) if($i ~ /^42$/) found=1} END{print found+0}' "$directory/affinity_preexec.txt")" = 1
  printf 'qualified\n' > "$directory/qualification_status"
}

{
  printf 'started_utc='; date -u +%Y-%m-%dT%H:%M:%SZ
  printf 'host='; hostname
  printf 'cpu=%s\n' "$cpu"
  printf 'cp2k_binary=%s\n' "$binary"
  printf 'cp2k_sha256=%s\n' "$expected_binary"
  printf 'libcp2k_sha256=%s\n' "$expected_libcp2k"
  printf 'libtblite_sha256=%s\n' "$expected_libtblite"
  printf 'cp2k_source_commit='; git -C /home/kuehne88/work/gxtb-final-clean-20260718/cp2k rev-parse HEAD
  printf 'save_tblite_source_commit=%s\n' "$source_commit"
  printf 'reference_cli=disabled\n'
} > "$work/provenance.env"

run_case H2O_0D NONE 20
for box in "${boxes[@]}"; do
  printf -v tag '%02d' "$box"
  run_case "H2O_L${tag}" XYZ "$box"
done

printf 'completed_utc=' >> "$work/provenance.env"
date -u +%Y-%m-%dT%H:%M:%SZ >> "$work/provenance.env"
printf 'complete\n' > "$work/run.status"
