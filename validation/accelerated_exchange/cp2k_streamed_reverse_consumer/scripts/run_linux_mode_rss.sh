#!/usr/bin/env bash
set -u

root=/home/kuehne88/work/cp2k-gxtb-streamed-reverse-consumer-20260716
envdir=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
runner="$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/run_with_rss.py"
input="$root/source/tests/xTB/regtest-tblite-gxtb-spglib/Si_prim_gxtb_kp_shifted_spglib.inp"
base="$root/evidence/linux_mode_rss"

export PATH="$envdir/bin:$PATH"
export LD_LIBRARY_PATH="$root/build/src:$envdir/lib:${LD_LIBRARY_PATH:-}"
export CP2K_DATA_DIR="$root/source/data"
export CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE=2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export LSAN_OPTIONS=detect_leaks=0
export OMPI_MCA_hwloc_base_binding_policy=none
ulimit -c 0

run_mode() {
  local mode=$1 cpu=$2
  local run="$base/${mode,,}_p1"
  mkdir -p "$run"
  cd "$run" || exit 90
  export CP2K_GXTB_EXCHANGE_GRADIENT_MODE="$mode"
  printf '%s\n' "$mode" > mode.txt
  printf '%s\n' "$cpu" > cpu_set.txt
  sha256sum "$input" > input.sha256
  set +e
  python3 "$runner" rusage.json timeout 1200 taskset -c "$cpu" \
    "$envdir/bin/mpiexec" --bind-to none -n 1 "$root/build/bin/cp2k.pdbg" \
    -i "$input" -o run.out > launcher.log 2>&1
  local rc=$?
  set -e
  printf '%s\n' "$rc" > returncode.txt
  grep -E "ENERGY\||Total energy:|GXTB STREAMED-REVERSE|PROGRAM ENDED" run.out \
    > selected_evidence.txt 2> /dev/null || true
  sha256sum mode.txt cpu_set.txt input.sha256 launcher.log rusage.json returncode.txt \
    > SHA256SUMS 2> /dev/null || true
  [[ -f run.out ]] && sha256sum run.out >> SHA256SUMS
  if [[ $rc -eq 0 ]] && grep -q "PROGRAM ENDED" run.out 2> /dev/null; then
    : > PASS
  else
    : > FAIL
  fi
}

mkdir -p "$base"
run_mode DENSE 50 &
pid_dense=$!
run_mode STREAMED 51 &
pid_streamed=$!
status=0
wait "$pid_dense" || status=1
wait "$pid_streamed" || status=1
date --iso-8601=seconds > "$base/finished.txt"
exit "$status"
