#!/usr/bin/env bash
set -u

test_root=/home/kuehne88/work/codex-gxtb-cp2k-block-helpers-20260716
baseline_root=/home/kuehne88/work/codex-gxtb-exchange-cache-cp2k-20260716
run_root=${test_root}/runtime_block_helpers_baseline
binary=${baseline_root}/build/bin/cp2k.psmp
mkdir -p "${run_root}"

launch() {
  name=$1
  input=$2
  cores=$3
  run=${run_root}/${name}
  mkdir -p "${run}"
  cp "${input}" "${run}/input.inp"
  sha256sum "${run}/input.inp" "${binary}" >"${run}/SHA256SUMS.initial"
  nohup bash -c '
    run=$1
    binary=$2
    data=$3
    cores=$4
    cd "${run}" || exit 97
    export OMP_NUM_THREADS=4
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export BLIS_NUM_THREADS=1
    export CP2K_DATA_DIR=${data}
    taskset -c "${cores}" "${binary}" -i input.inp -o cp2k.out
    rc=$?
    printf "%s\n" "${rc}" >returncode.txt
    sha256sum input.inp cp2k.out returncode.txt >SHA256SUMS.final
    exit "${rc}"
  ' runtime-job "${run}" "${binary}" "${baseline_root}/data" "${cores}" \
    >"${run}/launcher.log" 2>&1 </dev/null &
  printf "%s\n" "$!" >"${run}/launcher.pid"
}

launch time_reversal_311 \
  "${test_root}/tests/xTB/regtest-tblite-gxtb/H2_gxtb_kp_311_tr.inp" 200-203
launch k290_222_force_stress \
  "${test_root}/tests/xTB/regtest-tblite-gxtb/CH4_gxtb_kp_k290.inp" 204-207
launch spglib_shifted_222_force_stress \
  "${test_root}/tests/xTB/regtest-tblite-gxtb-spglib/Si_prim_gxtb_kp_shifted_spglib.inp" 208-211
