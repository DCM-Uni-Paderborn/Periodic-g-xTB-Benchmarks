#!/usr/bin/env bash
set -u

test_root=/home/kuehne88/work/codex-gxtb-cp2k-block-helpers-20260716
baseline_root=/home/kuehne88/work/codex-gxtb-exchange-cache-cp2k-20260716

launch() {
  run_root=$1
  binary=$2
  data=$3
  name=$4
  input=$5
  cores=$6
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
  ' runtime-job "${run}" "${binary}" "${data}" "${cores}" \
    >"${run}/launcher.log" 2>&1 </dev/null &
  printf "%s\n" "$!" >"${run}/launcher.pid"
}

new_runs=${test_root}/runtime_block_helpers_derivatives
baseline_runs=${test_root}/runtime_block_helpers_derivatives_baseline
mkdir -p "${new_runs}" "${baseline_runs}"

h2=${test_root}/tests/xTB/regtest-tblite-gxtb/H2_gxtb_kp_311_tr_force_stress.inp
ch4=${test_root}/tests/xTB/regtest-tblite-gxtb/CH4_gxtb_kp_k290_force_stress.inp

launch "${new_runs}" "${test_root}/build/bin/cp2k.psmp" "${test_root}/data" \
  time_reversal_311_force_stress "${h2}" 200-203
launch "${new_runs}" "${test_root}/build/bin/cp2k.psmp" "${test_root}/data" \
  k290_222_force_stress_debug "${ch4}" 204-207
launch "${baseline_runs}" "${baseline_root}/build/bin/cp2k.psmp" "${baseline_root}/data" \
  time_reversal_311_force_stress "${h2}" 208-211
launch "${baseline_runs}" "${baseline_root}/build/bin/cp2k.psmp" "${baseline_root}/data" \
  k290_222_force_stress_debug "${ch4}" 212-215
