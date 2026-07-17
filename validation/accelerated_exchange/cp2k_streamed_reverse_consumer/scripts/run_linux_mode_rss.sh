#!/usr/bin/env bash
set -u

if [[ $# -ne 1 ]]; then
  echo "usage: $0 CPU_FOR_DENSE,CPU_FOR_STREAMED" >&2
  exit 64
fi
IFS=, read -r -a reserved_cpus <<< "$1"
if [[ ${#reserved_cpus[@]} -ne 2 ]] || [[ ! ${reserved_cpus[0]} =~ ^[0-9]+$ ]] ||
  [[ ! ${reserved_cpus[1]} =~ ^[0-9]+$ ]] || [[ ${reserved_cpus[0]} == "${reserved_cpus[1]}" ]]; then
  echo "two distinct literal CPUs are required" >&2
  exit 65
fi

root=/home/kuehne88/work/cp2k-gxtb-streamed-reverse-consumer-20260716
envdir=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
runner="$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/run_with_rss.py"
input="$root/source/tests/xTB/regtest-tblite-gxtb-spglib/Si_prim_gxtb_kp_shifted_spglib.inp"
base="$root/evidence/linux_mode_rss"

mkdir -p "$base"
exec 9> "$base/.run_linux_mode_rss.lock"
if ! flock -n 9; then
  echo "another RSS-mode writer already owns $base" >&2
  exit 75
fi
cpu_lock_root="/tmp/periodic-gxtb-cpu-reservations-$(id -u)"
mkdir -p "$cpu_lock_root"
reservation_lock_fds=()
for cpu in "${reserved_cpus[@]}"; do
  exec {lock_fd}>> "$cpu_lock_root/cpu-${cpu}.lock"
  if ! flock -n "$lock_fd"; then
    echo "logical CPU $cpu is already reserved by another production launcher" >&2
    exit 76
  fi
  printf '{"cpu":%s,"pid":%s,"source":"run_linux_mode_rss.sh"}\n' \
    "$cpu" "$$" >&"$lock_fd"
  reservation_lock_fds+=("$lock_fd")
done

export PATH="$envdir/bin:$PATH"
export LD_LIBRARY_PATH="$root/build/src:$envdir/lib:${LD_LIBRARY_PATH:-}"
export CP2K_DATA_DIR="$root/source/data"
export CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE=2
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export LSAN_OPTIONS=detect_leaks=0
for inherited_mca_key in ${!OMPI_MCA_@} ${!PRTE_MCA_@}; do
  unset "$inherited_mca_key"
done
ulimit -c 0

run_mode() {
  local mode=$1 cpu=$2
  local run="$base/${mode,,}_p1"
  mkdir -p "$run"
  cd "$run" || exit 90
  rm -f PASS FAIL run.out launcher.log rusage.json returncode.txt \
    selected_evidence.txt SHA256SUMS
  export CP2K_GXTB_EXCHANGE_GRADIENT_MODE="$mode"
  printf '%s\n' "$mode" > mode.txt
  printf '%s\n' "$cpu" > ordered_pe_list.txt
  sha256sum "$input" > input.sha256
  set +e
  python3 "$runner" --mpi-ranks 1 --ordered-pe-list "$cpu" \
    --cp2k "$root/build/bin/cp2k.pdbg" --launcher-log launcher.log \
    rusage.json timeout 1200 \
    "$envdir/bin/mpiexec" --map-by "pe-list=${cpu}:ordered" --bind-to core \
    --report-bindings -n 1 "$root/build/bin/cp2k.pdbg" \
    -i "$input" -o run.out > launcher.log 2>&1
  local rc=$?
  set -e
  printf '%s\n' "$rc" > returncode.txt
  grep -E "ENERGY\||Total energy:|GXTB STREAMED-REVERSE|PROGRAM ENDED" run.out \
    > selected_evidence.txt 2> /dev/null || true
  if [[ $rc -eq 0 ]] && grep -q "PROGRAM ENDED" run.out 2> /dev/null; then
    : > PASS
  else
    : > FAIL
    if [[ $rc -eq 0 ]]; then
      rc=98
      printf '%s\n' "$rc" > returncode.txt
    fi
  fi
  sha256sum mode.txt ordered_pe_list.txt input.sha256 launcher.log rusage.json returncode.txt \
    > SHA256SUMS 2> /dev/null || true
  [[ -f run.out ]] && sha256sum run.out >> SHA256SUMS
  [[ -f PASS ]] && sha256sum PASS >> SHA256SUMS
  [[ -f FAIL ]] && sha256sum FAIL >> SHA256SUMS
  return "$rc"
}

run_mode DENSE "${reserved_cpus[0]}" &
pid_dense=$!
run_mode STREAMED "${reserved_cpus[1]}" &
pid_streamed=$!
status=0
wait "$pid_dense" || status=1
wait "$pid_streamed" || status=1
date --iso-8601=seconds > "$base/finished.txt"
exit "$status"
