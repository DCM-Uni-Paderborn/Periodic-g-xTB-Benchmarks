# Reproduction commands

The campaign used clean out-of-tree builds. Shell continuations below are only
for readability; the hashes in `build_manifest.json` identify the exact
artifacts used for all accepted calculations.

## save_tblite

```sh
cmake -S /tmp/save_tblite_cp2k \
  -B /tmp/save_tblite_1449feb_bench_build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/tmp/save_tblite_1449feb_bench_install \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-16 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DBUILD_SHARED_LIBS=OFF -DWITH_OpenMP=ON \
  -DWITH_DDX=OFF -DWITH_HDF5=OFF -DWITH_TREXIO=OFF \
  -DWITH_TESTS=ON -Dtblite-dependency-method=fetch \
  -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build /tmp/save_tblite_1449feb_bench_build -j 8
cmake --install /tmp/save_tblite_1449feb_bench_build
```

The full CTest run passed 88/93 wrapper tests. The g-XTB wrapper passed. Three
individual finite-difference cases exceed machine-precision-derived thresholds
only narrowly: `2.24436e-9` versus `2.22045e-9`, `2.49000e-10` versus
`2.22045e-10`, and `2.87804e-10` versus `2.22045e-10`. The remaining C-API and
XTBML failures are exclusively DDX cases in a build configured with DDX off.
The campaign contains no solvent calculations.

## CP2K

```sh
cmake -S /tmp/cp2k_gxtb \
  -B /tmp/cp2k_gxtb_18d37c_bench_build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-16 \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_PREFIX_PATH='/tmp/save_tblite_1449feb_bench_install;/opt/homebrew/opt/openblas;/opt/homebrew/opt/scalapack;/opt/homebrew/opt/spglib;/opt/homebrew/opt/dbcsr;/opt/homebrew/opt/open-mpi' \
  -DCP2K_USE_MPI=ON -DCP2K_USE_TBLITE=ON \
  -DCP2K_TBLITE_PROVIDER=SAVE -DCP2K_USE_SPGLIB=ON \
  -DCP2K_BLAS_VENDOR=OpenBLAS -DCP2K_SCALAPACK_VENDOR=GENERIC \
  -DCP2K_USE_FFTW3=OFF -DCP2K_USE_LIBXC=OFF \
  -DCP2K_USE_DFTD4=OFF -DCP2K_USE_HDF5=OFF \
  -DCP2K_USE_LIBINT2=OFF -DCP2K_USE_ELPA=OFF \
  -DCP2K_WITH_GPU=NONE
cmake --build /tmp/cp2k_gxtb_18d37c_bench_build \
  -j 8 --target cp2k.psmp
```

The benchmark environment pins all nested threading:

```sh
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_WAIT_POLICY=PASSIVE
```

The exact build passed all 49 core g-XTB CP2K matchers and all 9 SPGLIB
g-XTB matchers. The CLI-reference matchers require
`/tmp/save_tblite_1449feb_bench_install/bin` in `PATH`.
