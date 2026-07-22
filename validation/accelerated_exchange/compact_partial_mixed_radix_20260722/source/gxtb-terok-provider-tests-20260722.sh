#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/kuehne88/work/gxtb-compact-partial-combination-20260722
ENVROOT=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
BUILD=${ROOT}/provider-build-r3

export PATH=${ENVROOT}/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=${ROOT}/provider-install/lib:${ENVROOT}/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "${ROOT}/logs" \
  "${ROOT}/provenance/provider-exchange-preexec" \
  "${ROOT}/provenance/provider-gxtb-preexec"

grep -E '^(MemTotal|MemAvailable):' /proc/meminfo \
  > "${ROOT}/provenance/provider-tests-prelaunch-memory.txt"
ps -e -o pid=,ppid=,user=,stat=,psr=,rss=,etimes=,comm=,args= \
  > "${ROOT}/provenance/provider-tests-prelaunch-all-rss.txt"

cd "${BUILD}/test/unit"
"${ROOT}/mpi_singleton_launch.sh" 142 \
  "${ROOT}/provenance/provider-exchange-preexec" \
  "${BUILD}/test/unit/tblite-tester" exchange \
  > "${ROOT}/logs/provider-exchange-direct.log" 2>&1

"${ROOT}/mpi_singleton_launch.sh" 142 \
  "${ROOT}/provenance/provider-gxtb-preexec" \
  "${BUILD}/test/unit/tblite-tester" gxtb \
  > "${ROOT}/logs/provider-gxtb-direct.log" 2>&1

sha256sum \
  "${ROOT}/logs/provider-exchange-direct.log" \
  "${ROOT}/logs/provider-gxtb-direct.log" \
  > "${ROOT}/provenance/provider-direct-output-sha256.txt"
