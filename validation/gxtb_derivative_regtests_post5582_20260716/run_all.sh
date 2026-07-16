#!/usr/bin/env bash
set -u

campaign_root=$(cd "$(dirname "$0")" && pwd)
cp2k_root=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714
cp2k_bin=${cp2k_root}/build/cp2k/bin/cp2k.psmp
mpirun_bin=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env/bin/mpirun
export CP2K_DATA_DIR=${cp2k_root}/data
export OMP_NUM_THREADS=1

mkdir -p "${campaign_root}/runs"
: > "${campaign_root}/launched.tsv"

for input in "${campaign_root}"/inputs/*.inp; do
  case_name=$(basename "${input}" .inp)
  run_dir=${campaign_root}/runs/${case_name}
  mkdir -p "${run_dir}"
  cp "${input}" "${run_dir}/input.inp"
  (
    cd "${run_dir}" || exit 97
    sha256sum input.inp "${cp2k_bin}" > SHA256SUMS.initial
    "${mpirun_bin}" --bind-to none -np 4 "${cp2k_bin}" -i input.inp -o cp2k.out > launcher.log 2>&1
    rc=$?
    printf '%s\n' "${rc}" > returncode.txt
    sha256sum input.inp cp2k.out launcher.log returncode.txt > SHA256SUMS.final
    exit "${rc}"
  ) &
  child=$!
  printf '%s\t%s\n' "${case_name}" "${child}" >> "${campaign_root}/launched.tsv"
done

wait
date -u +%Y-%m-%dT%H:%M:%SZ > "${campaign_root}/completed_utc.txt"
