#!/usr/bin/env bash
set -u

root=/home/kuehne88/work/cp2k-gxtb-streamed-reverse-consumer-20260716
envdir=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
runner="$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/run_with_rss.py"
matrix="$root/evidence/linux_matrix"

export PATH="$envdir/bin:$PATH"
export LD_LIBRARY_PATH="$root/build/src:$envdir/lib:${LD_LIBRARY_PATH:-}"
export CP2K_DATA_DIR="$root/source/data"
export CP2K_GXTB_EXCHANGE_GRADIENT_MODE=QUALIFY
export CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE=2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export LSAN_OPTIONS=detect_leaks=0
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_rmaps_base_oversubscribe=1
ulimit -c 0

names=(
  k290_rks_3d_fd
  shifted_spglib_rks_3d
  tr_rks_3d_fd
  fullmesh_rks_3d_fd
  fullmesh_uks_3d_fd
  rks_1d_fd
  spglib_rks_2d_fd
)
inputs=(
  "$root/source/tests/xTB/regtest-tblite-gxtb/CH4_gxtb_kp_k290_force_stress.inp"
  "$root/source/tests/xTB/regtest-tblite-gxtb-spglib/Si_prim_gxtb_kp_shifted_spglib.inp"
  "$root/source/tests/xTB/regtest-tblite-gxtb/H2_gxtb_kp_311_tr_force_stress.inp"
  "$root/source/tests/xTB/regtest-tblite-gxtb/H2_gxtb_kp_311_force_stress.inp"
  "$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/inputs/O2_gxtb_uks_kp_311_force_stress.inp"
  "$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/inputs/Ar2_chain_gxtb_x_kp_debug.inp"
  "$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/inputs/Ar_layer_gxtb_yz_kp_debug.inp"
)

run_one() {
  local name=$1 input=$2 ranks=$3 first_cpu=$4 last_cpu=$5
  local run="$matrix/${name}_p${ranks}"
  mkdir -p "$run"
  cd "$run" || exit 90
  printf '%s\n' "$input" > input.path
  printf '%s-%s\n' "$first_cpu" "$last_cpu" > cpu_set.txt
  printf '%s\n' "$ranks" > mpi_ranks.txt
  sha256sum "$input" > input.sha256
  date --iso-8601=seconds > started.txt
  : > RUNNING

  set +e
  python3 "$runner" rusage.json timeout 1200 taskset -c "${first_cpu}-${last_cpu}" \
    "$envdir/bin/mpiexec" --bind-to none -n "$ranks" \
    "$root/build/bin/cp2k.pdbg" -i "$input" -o run.out > launcher.log 2>&1
  local rc=$?
  set -e
  printf '%s\n' "$rc" > returncode.txt
  date --iso-8601=seconds > finished.txt
  if [[ -f run.out ]]; then
    grep -E "GXTB-QUALIFICATION_ONLY STREAMED-REVERSE|GXTB STREAMED-REVERSE|DEBUG.*Sum of differences|PROGRAM ENDED" \
      run.out > selected_evidence.txt || true
  fi
  sha256sum input.path input.sha256 cpu_set.txt mpi_ranks.txt launcher.log rusage.json \
    returncode.txt started.txt finished.txt > SHA256SUMS 2> /dev/null || true
  [[ -f run.out ]] && sha256sum run.out >> SHA256SUMS
  rm -f RUNNING
  if [[ $rc -eq 0 ]] && grep -q "PROGRAM ENDED" run.out 2> /dev/null &&
    grep -q "GXTB-QUALIFICATION_ONLY STREAMED-REVERSE" run.out 2> /dev/null; then
    : > PASS
  else
    : > FAIL
  fi
}

mkdir -p "$matrix"
next_cpu=0
pids=()
for index in "${!names[@]}"; do
  for ranks in 1 2 4; do
    first_cpu=$next_cpu
    last_cpu=$((next_cpu + ranks - 1))
    next_cpu=$((last_cpu + 1))
    run_one "${names[$index]}" "${inputs[$index]}" "$ranks" "$first_cpu" "$last_cpu" &
    pids+=("$!")
  done
done

overall=0
for pid in "${pids[@]}"; do
  wait "$pid" || overall=1
done

find "$matrix" -mindepth 2 -maxdepth 2 -type f \( -name PASS -o -name FAIL \) -print | sort \
  > "$matrix/status.txt"
date --iso-8601=seconds > "$matrix/finished.txt"
exit "$overall"
