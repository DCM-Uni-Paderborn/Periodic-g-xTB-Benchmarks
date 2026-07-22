#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  printf 'usage: %s CASE INPUT CPU_BASE\n' "$0" >&2
  exit 2
fi

CASE_NAME=$1
INPUT_SOURCE=$2
CPU_BASE=$3

ROOT=/home/kuehne88/work/gxtb-compact-partial-combination-20260722
ENVROOT=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
DEPS=/home/kuehne88/work/gxtb-symmetry-fused-combination-20260721/save_tblite-install
BINARY=${ROOT}/cp2k-build-r1/bin/cp2k.psmp
RESULT=${ROOT}/results/${CASE_NAME}
PROOF=${RESULT}/provenance/affinity

export PATH=${ENVROOT}/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=${ROOT}/provider-install/lib:${DEPS}/lib:${ENVROOT}/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "${RESULT}" "${PROOF}"
cp "${INPUT_SOURCE}" "${RESULT}/input.inp"

{
  date -u +'%Y-%m-%dT%H:%M:%SZ'
  hostname
  uname -a
} > "${RESULT}/provenance/host.txt"

env | LC_ALL=C sort > "${RESULT}/provenance/environment.txt"
grep -E '^(MemTotal|MemAvailable):' /proc/meminfo \
  > "${RESULT}/provenance/prelaunch-memory.txt"
ps -e -o pid=,ppid=,user=,stat=,psr=,rss=,etimes=,comm=,args= \
  > "${RESULT}/provenance/prelaunch-all-rss.txt"
sha256sum \
  "${BINARY}" \
  "${ROOT}/provider-install/lib/libtblite.a" \
  "${RESULT}/input.inp" \
  "${ROOT}/cp2k-src/src/cp_control_utils.F" \
  "${ROOT}/cp2k-src/src/tblite_interface.F" \
  "${ROOT}/provider-src/src/tblite/cp2k_compat.f90" \
  > "${RESULT}/provenance/pre-run-sha256.txt"

cd "${RESULT}"
"${ENVROOT}/bin/mpirun" --bind-to none -np 2 \
  "${ROOT}/mpi_singleton_launch.sh" "${CPU_BASE}" "${PROOF}" \
  "${BINARY}" -i input.inp -o cp2k.out \
  > mpi-launch.out 2>&1

sha256sum cp2k.out mpi-launch.out > provenance/output-sha256.txt
