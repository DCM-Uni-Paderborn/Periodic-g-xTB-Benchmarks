#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-acp-stream-20260722
envroot=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
reuse_prefix=/home/kuehne88/work/gxtb-range-fft-fix-20260722/save_tblite-install-r2
cpu=141
candidate_kib=$((16*1024*1024))
native_peak_kib=$((100*1024*1024))
minimum_margin_kib=$((128*1024*1024))
provider_build="$root/save_tblite-build-r3"
provider_install="$root/save_tblite-install-r3"
cp2k_build="$root/cp2k-build-r3"
log="$root/logs/build-r3.log"
exit_file="$root/logs/build-r3.exit"

trap 'rc=$?; printf "%s\n" "$rc" > "$exit_file"' EXIT

export PATH="$envroot/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$envroot/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export BLIS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

printf 'pid=%s expected_cpu=%s\n' "$$" "$cpu" > "$root/provenance/build-r3-affinity-preexec.txt"
grep -E '^(Name|State|Cpus_allowed|Cpus_allowed_list):' /proc/$$/status \
  >> "$root/provenance/build-r3-affinity-preexec.txt"
if [[ $(awk '/^Cpus_allowed_list:/ {print $2}' /proc/$$/status) != "$cpu" ]]; then
  printf 'build launcher is not bound to singleton CPU %s\n' "$cpu" >&2
  exit 74
fi

grep -E '^(MemAvailable|MemTotal):' /proc/meminfo \
  > "$root/provenance/prelaunch-memory-r3.txt"
ps -e -o pid=,ppid=,sid=,psr=,rss=,vsz=,nlwp=,etimes=,stat=,comm=,args= --sort=-rss \
  > "$root/provenance/prelaunch-all-rss-r3.tsv"

mem_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
remaining_kib=0
live_cp2k=0
while read -r pid; do
  [[ -r /proc/$pid/status ]] || continue
  state=$(awk '/^State:/ {print $2}' /proc/$pid/status)
  rss=$(awk '/^VmRSS:/ {print $2}' /proc/$pid/status)
  [[ $state != Z && ${rss:-0} -gt 0 ]] || continue
  allowance=$((native_peak_kib > rss ? native_peak_kib-rss : 0))
  remaining_kib=$((remaining_kib+allowance))
  live_cp2k=$((live_cp2k+1))
done < <(pgrep -x cp2k.psmp || true)
margin_kib=$((mem_kib-remaining_kib-candidate_kib))
{
  printf 'mem_available_kib=%s\n' "$mem_kib"
  printf 'live_cp2k=%s\n' "$live_cp2k"
  printf 'remaining_growth_allowance_kib=%s\n' "$remaining_kib"
  printf 'candidate_build_peak_kib=%s\n' "$candidate_kib"
  printf 'computed_margin_kib=%s\n' "$margin_kib"
  printf 'minimum_margin_kib=%s\n' "$minimum_margin_kib"
} > "$root/provenance/prelaunch-budget-r3.txt"
if ((margin_kib < minimum_margin_kib)); then
  printf 'build launch margin is below 128 GiB\n' >&2
  exit 75
fi

git -C "$root/save_tblite-src" diff --check
git -C "$root/cp2k" diff --check

for library in dftd4 mctc-lib multicharge s-dftd3 toml-f; do
  test -f "$reuse_prefix/lib/lib${library}.a"
  test -d "$reuse_prefix/lib/cmake/$library"
done
sha256sum "$reuse_prefix/lib"/lib{dftd4,mctc-lib,multicharge,s-dftd3,toml-f}.a \
  > "$root/provenance/reused-dependencies-r3.sha256"
grep -E '^(CMAKE_(C|Fortran)_COMPILER|CMAKE_BUILD_TYPE|CMAKE_INSTALL_PREFIX)' \
  /home/kuehne88/work/gxtb-range-fft-fix-20260722/save_tblite-build-r2/CMakeCache.txt \
  > "$root/provenance/reused-dependencies-r3-cmake.txt"

"$envroot/bin/cmake" \
  -S "$root/save_tblite-src" \
  -B "$provider_build" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$envroot/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_Fortran_COMPILER="$envroot/bin/x86_64-conda-linux-gnu-gfortran" \
  -DCMAKE_INSTALL_PREFIX="$provider_install" \
  "-DCMAKE_PREFIX_PATH=$reuse_prefix;$envroot" \
  -DBUILD_SHARED_LIBS=OFF \
  -DWITH_API=ON \
  -DWITH_BLAS=ON \
  -DWITH_DDX=OFF \
  -DWITH_HDF5=OFF \
  -DWITH_JSON=OFF \
  -DWITH_OpenMP=ON \
  -DWITH_TESTS=ON \
  -DWITH_TREXIO=OFF \
  > "$log" 2>&1
"$envroot/bin/cmake" --build "$provider_build" --parallel 1 >> "$log" 2>&1
"$envroot/bin/ctest" --test-dir "$provider_build" --output-on-failure -V \
  -R '^tblite/(acp|gxtb)$' > "$root/logs/provider-acp-gxtb-release-r3.log" 2>&1
"$envroot/bin/cmake" --install "$provider_build" >> "$log" 2>&1

"$envroot/bin/cmake" \
  -S "$root/cp2k" \
  -B "$cp2k_build" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$envroot/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$envroot/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="$envroot/bin/mpifort" \
  "-DCMAKE_PREFIX_PATH=$provider_install;$reuse_prefix;$envroot" \
  -DBUILD_SHARED_LIBS=ON \
  -DCP2K_USE_MPI=ON \
  -DCP2K_USE_OPENMP=ON \
  -DCP2K_USE_SCALAPACK=ON \
  -DCP2K_USE_SPGLIB=ON \
  -DCP2K_USE_TBLITE=ON \
  -DCP2K_TBLITE_PROVIDER=SAVE \
  -DCP2K_TBLITE_REVISION=acp-stream-r3 \
  -DCP2K_USE_ACCEL=NONE \
  -DCP2K_WITH_GPU=NONE \
  -DCP2K_USE_FFTW3=OFF \
  -DCP2K_USE_DFTD4=OFF \
  -DCP2K_USE_HDF5=OFF \
  -DCP2K_USE_LIBINT2=OFF \
  -DCP2K_USE_LIBXC=OFF \
  >> "$log" 2>&1
"$envroot/bin/cmake" --build "$cp2k_build" --target cp2k.psmp --parallel 1 \
  >> "$log" 2>&1

sha256sum \
  "$provider_install/lib/libtblite.a" \
  "$cp2k_build/bin/cp2k.psmp" \
  "$cp2k_build/src/libcp2k.so.2026.2" \
  "$root/provenance/provider.patch" \
  "$root/provenance/cp2k.patch" \
  > "$root/provenance/build-r3.sha256"
ldd "$cp2k_build/bin/cp2k.psmp" > "$root/provenance/build-r3-ldd.txt"
git -C "$root/save_tblite-src" status --short \
  > "$root/provenance/provider-status-r3.txt"
git -C "$root/cp2k" status --short \
  > "$root/provenance/cp2k-status-r3.txt"
