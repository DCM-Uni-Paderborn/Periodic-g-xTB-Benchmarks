#!/usr/bin/env bash
set -u

root=/home/kuehne88/work/gxtb-partial-root-final-35e7942-20260716/fault-hooks-release
inputs=/home/kuehne88/work/gxtb-partial-root-final-35e7942-20260716/inputs
build=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_35e7942_release_20260716
provider=/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942_install_release
env_prefix=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
dep=/home/kuehne88/work/codex-gxtb-acpfix-20260714/install/save_tblite
source=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_clean_20260716
binary="$build/bin/cp2k.psmp"

export PATH="$env_prefix/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$build/src:$provider/lib:$env_prefix/lib:$dep/lib"
export CP2K_DATA_DIR="$source/data"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export CP2K_GXTB_EXCHANGE_STREAM_MODE=KGROUP_PARTIAL_ROOT
export CP2K_GXTB_EXCHANGE_GRADIENT_MODE=QUALIFY
export CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION=1
export CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE=3
export CP2K_GXTB_EXCHANGE_TRANSFORM_MODE=SEPARABLE

mkdir -p "$root"
printf 'case\tselector\tnp\texit_code\ttimed_out\tdiagnostic\tpass\n' > "$root/status.tsv"

run_case() {
  name=$1
  selector=$2
  np=$3
  input=$4
  diagnostic=$5
  case_dir="$root/$name"
  mkdir -p "$case_dir"
  cp "$inputs/$input" "$case_dir/input.inp"
  export CP2K_GXTB_PARTIAL_QUALIFY_INJECT="$selector"
  (
    cd "$case_dir" || exit 99
    timeout --signal=TERM --kill-after=10s 180s taskset -c 160-175 \
      mpiexec -n "$np" "$binary" -i input.inp > run.log 2>&1
  )
  rc=$?
  timed_out=0
  if test "$rc" -eq 124 -o "$rc" -eq 137; then timed_out=1; fi
  found=0
  if grep -Fq "$diagnostic" "$case_dir/run.log"; then found=1; fi
  pass=0
  if test "$rc" -ne 0 -a "$timed_out" -eq 0 -a "$found" -eq 1; then pass=1; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$name" "$selector" "$np" "$rc" "$timed_out" "$found" "$pass" >> "$root/status.tsv"
}

run_case nonleader_failure NONLEADER_FAILURE 4 ch4_k222_pgs2.inp \
  'save_tblite KGROUP_PARTIAL_ROOT initialization failed'
run_case antihermitian_reverse ANTIHERMITIAN_REVERSE 2 ch4_k222_default.inp \
  'Folded KGROUP_PARTIAL_ROOT overlap adjoint is not Hermitian'
run_case nonfinite_forward_result NONFINITE_FORWARD_RESULT 2 ch4_k222_default.inp \
  'save_tblite KGROUP_PARTIAL_ROOT result failed'
run_case unknown_selector BOGUS 2 ch4_k222_default.inp \
  'Unknown CP2K_GXTB_PARTIAL_QUALIFY_INJECT value'
run_case truncated_selector XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 2 \
  ch4_k222_default.inp 'CP2K_GXTB_PARTIAL_QUALIFY_INJECT is truncated or unreadable'

unset CP2K_GXTB_PARTIAL_QUALIFY_INJECT
sha256sum "$binary" "$provider/lib/libtblite.a" "$source/src/tblite_interface.F" \
  "$inputs"/*.inp "$root"/*/run.log > "$root/SHA256SUMS"
awk -F '\t' 'NR > 1 && $7 != 1 {failed = 1} END {exit failed}' "$root/status.tsv"
