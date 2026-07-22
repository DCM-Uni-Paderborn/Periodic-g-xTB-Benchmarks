#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/kuehne88/work/gxtb-compact-partial-combination-20260722
ENVROOT=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
DEPENDENCY_PREFIX=/home/kuehne88/work/gxtb-symmetry-fused-combination-20260721/save_tblite-install

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PATH="${ENVROOT}/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="${ENVROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

"${ENVROOT}/bin/cmake" \
  -S "${ROOT}/provider-src" \
  -B "${ROOT}/provider-build-r3" \
  -G Ninja \
  -DCMAKE_MAKE_PROGRAM="${ENVROOT}/bin/ninja" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${ENVROOT}/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_Fortran_COMPILER="${ENVROOT}/bin/x86_64-conda-linux-gnu-gfortran" \
  -DCMAKE_INSTALL_PREFIX="${ROOT}/provider-install" \
  -DCMAKE_PREFIX_PATH="${DEPENDENCY_PREFIX}" \
  -DWITH_TESTS=ON \
  -DWITH_DDX=OFF \
  2>&1 | tee "${ROOT}/logs/provider-configure.log"

"${ENVROOT}/bin/ninja" -C "${ROOT}/provider-build-r3" -j1 \
  2>&1 | tee "${ROOT}/logs/provider-build.log"

"${ENVROOT}/bin/ctest" \
  --test-dir "${ROOT}/provider-build-r3" \
  -R "^tblite/(exchange|gxtb)$" \
  --output-on-failure \
  2>&1 | tee "${ROOT}/logs/provider-tests.log"

"${ENVROOT}/bin/cmake" --install "${ROOT}/provider-build-r3" \
  2>&1 | tee "${ROOT}/logs/provider-install.log"
