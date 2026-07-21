#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  printf 'usage: %s CPU SOURCE_DIR RUNTIME_DIR DEPENDENCY_DIR RESULT_ROOT\n' "$0" >&2
  exit 64
fi

cpu=$1
source_dir=$2
runtime=$3
deps=$4
root=$5
expected_commit=9322cc6d43be5099f2e7dc4866abc45b8387835a

if [[ ! $cpu =~ ^[0-9]+$ ]]; then
  printf 'invalid CPU: %s\n' "$cpu" >&2
  exit 64
fi
test -d "$source_dir/.git"
test -x "$runtime/bin/cmake"
test "$(git -C "$source_dir" rev-parse HEAD)" = "$expected_commit"
test -z "$(git -C "$source_dir" status --porcelain)"

build=$root/build
logs=$root/logs
provenance=$root/provenance
reservation_root="/tmp/gxtb-cpu-reservations-${USER}"
mkdir -p "$build" "$logs" "$provenance" "$reservation_root"

exec 9>"$reservation_root/.lock"
flock -x 9
for reservation_file in "$reservation_root"/*.reservation; do
  [[ -e $reservation_file ]] || continue
  read -r reserved_pid reserved_cpu < "$reservation_file" || true
  if [[ -z ${reserved_pid:-} || ! -d /proc/$reserved_pid ]]; then
    rm -f "$reservation_file"
    continue
  fi
  state=$(awk '/^State:/{print $2}' "/proc/$reserved_pid/status" 2>/dev/null || true)
  if [[ $state == Z ]]; then
    rm -f "$reservation_file"
    continue
  fi
  if [[ $reserved_cpu == "$cpu" ]]; then
    printf 'CPU %s is reserved by PID %s\n' "$cpu" "$reserved_pid" >&2
    exit 75
  fi
done
reservation="$reservation_root/provider-head-linux-${expected_commit:0:8}.reservation"
printf '%s %s\n' "$$" "$cpu" > "$reservation"
flock -u 9
trap 'rm -f "$reservation"' EXIT INT TERM

export PATH="$runtime/bin:$PATH"
export LD_LIBRARY_PATH="$runtime/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$runtime/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

taskset -c "$cpu" bash -c '
  set -euo pipefail
  expected_cpu=$1
  proof=$2
  shift 2
  allowed=$(awk "/^Cpus_allowed_list:/{print \$2}" "/proc/$$/status")
  {
    printf "pid=%s expected_cpu=%s allowed=%s\\n" "$$" "$expected_cpu" "$allowed"
    grep -E "^(Name|State|Cpus_allowed|Cpus_allowed_list):" "/proc/$$/status"
    env | grep -E "^(OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|BLIS_NUM_THREADS|VECLIB_MAXIMUM_THREADS|NUMEXPR_NUM_THREADS)=" | sort
  } > "$proof"
  [[ $allowed == "$expected_cpu" ]]
  exec "$@"
' bash "$cpu" "$provenance/affinity_preexec.txt" "$runtime/bin/cmake" \
  -S "$source_dir" -B "$build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_SKIP_RPATH=ON \
  -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=ON -DWITH_TESTS=ON \
  -DWITH_DDX=OFF -DWITH_HDF5=OFF -DWITH_JSON=OFF -DWITH_TREXIO=OFF \
  -DCMAKE_PREFIX_PATH="$runtime" \
  -DCMAKE_C_COMPILER="$runtime/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_Fortran_COMPILER="$runtime/bin/x86_64-conda-linux-gnu-gfortran" \
  -DFETCHCONTENT_SOURCE_DIR_DFTD4="$deps/dftd4" \
  -DFETCHCONTENT_SOURCE_DIR_MCTC-LIB="$deps/mctc-lib" \
  -DFETCHCONTENT_SOURCE_DIR_MSTORE="$deps/mstore" \
  -DFETCHCONTENT_SOURCE_DIR_MULTICHARGE="$deps/multicharge" \
  -DFETCHCONTENT_SOURCE_DIR_S-DFTD3="$deps/s-dftd3" \
  -DFETCHCONTENT_SOURCE_DIR_TEST-DRIVE="$deps/test-drive" \
  -DFETCHCONTENT_SOURCE_DIR_TOML-F="$deps/toml-f" \
  > "$logs/configure.log" 2>&1

taskset -c "$cpu" "$runtime/bin/cmake" --build "$build" -j1 \
  > "$logs/build.log" 2>&1

set +e
taskset -c "$cpu" "$runtime/bin/ctest" --test-dir "$build" \
  --output-on-failure --timeout 300 -j1 > "$logs/ctest-full.log" 2>&1
ctest_rc=$?
set -e
printf '%s\n' "$ctest_rc" > "$logs/ctest-full.rc"

taskset -c "$cpu" "$runtime/bin/ctest" --test-dir "$build" \
  --output-on-failure --timeout 300 -j1 \
  -R '(^pbc$|^wignerseitz$|s-dftd3/periodic|tblite/(acp|exchange|gxtb|hamiltonian|integral-trafo|mixer|q-vszp|wavefunction-restart|wignerseitz))' \
  > "$logs/ctest-focused.log" 2>&1

{
  printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'source_commit=%s\n' "$(git -C "$source_dir" rev-parse HEAD)"
  printf 'source_branch=%s\n' "$(git -C "$source_dir" branch --show-current)"
  printf 'ctest_full_rc=%s\n' "$ctest_rc"
  printf 'ctest_focused_rc=0\n'
} > "$root/STATUS"

exit "$ctest_rc"
