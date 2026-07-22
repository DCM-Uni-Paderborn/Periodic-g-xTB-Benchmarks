#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/kuehne88/work/gxtb-compact-partial-combination-20260722
ENVROOT=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
DEPS=/home/kuehne88/work/gxtb-symmetry-fused-combination-20260721/save_tblite-install
SOURCE=${ROOT}/cp2k-src
BUILD=${ROOT}/cp2k-build-r1
INSTALL=${ROOT}/cp2k-install

export PATH=${ENVROOT}/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=${ROOT}/provider-install/lib:${DEPS}/lib:${ENVROOT}/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "${BUILD}" "${INSTALL}" "${ROOT}/logs" "${ROOT}/provenance"

"${ENVROOT}/bin/cmake" -S "${SOURCE}" -B "${BUILD}" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${INSTALL}" \
  -DCMAKE_MAKE_PROGRAM="${ENVROOT}/bin/ninja" \
  -DCMAKE_C_COMPILER="${ENVROOT}/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="${ENVROOT}/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="${ENVROOT}/bin/mpifort" \
  -DCMAKE_PREFIX_PATH="${ROOT}/provider-install;${DEPS};${ENVROOT}" \
  -DCP2K_USE_MPI=ON \
  -DCP2K_USE_OPENMP=ON \
  -DCP2K_USE_SCALAPACK=ON \
  -DCP2K_USE_TBLITE=ON \
  -DCP2K_TBLITE_PROVIDER=SAVE \
  -DCP2K_USE_SPGLIB=ON \
  -DCP2K_USE_FFTW3=OFF \
  -DCP2K_USE_LIBXC=OFF \
  -DCP2K_USE_DFTD4=OFF \
  -DCP2K_USE_ACCEL=NONE \
  -DCP2K_WITH_GPU=NONE \
  -DDBCSR_DIR="${ENVROOT}/lib/cmake/dbcsr" \
  -DSpglib_DIR="${ENVROOT}/lib/cmake/Spglib" \
  -Dtblite_DIR="${ROOT}/provider-install/lib/cmake/tblite" \
  -Dmctc-lib_DIR="${DEPS}/lib/cmake/mctc-lib" \
  -Dtoml-f_DIR="${DEPS}/lib/cmake/toml-f" \
  -Ds-dftd3_DIR="${DEPS}/lib/cmake/s-dftd3" \
  -Ddftd4_DIR="${DEPS}/lib/cmake/dftd4" \
  -Dmulticharge_DIR="${DEPS}/lib/cmake/multicharge" \
  2>&1 | tee "${ROOT}/logs/cp2k-configure.log"

"${ENVROOT}/bin/cmake" --build "${BUILD}" --target cp2k.psmp -j1 \
  2>&1 | tee "${ROOT}/logs/cp2k-build.log"

sha256sum \
  "${BUILD}/bin/cp2k.psmp" \
  "${ROOT}/provider-install/lib/libtblite.a" \
  > "${ROOT}/provenance/qualified-binary-sha256.txt"
