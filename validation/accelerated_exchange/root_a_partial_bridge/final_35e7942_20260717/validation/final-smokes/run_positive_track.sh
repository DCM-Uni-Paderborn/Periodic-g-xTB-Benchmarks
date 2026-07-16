#!/usr/bin/env bash
set -u

tag=${1:?usage: run_positive_track.sh release|debug}
root_base=/home/kuehne88/work/gxtb-partial-root-final-35e7942-20260716
inputs="$root_base/inputs"
source=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_clean_20260716
env_prefix=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
dep=/home/kuehne88/work/codex-gxtb-acpfix-20260714/install/save_tblite

case "$tag" in
  release)
    build=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_35e7942_release_20260716
    provider=/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942_install_release
    binary="$build/bin/cp2k.psmp"
    cpu_set=144-159
    ;;
  debug)
    build=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_35e7942_debug_20260716
    provider=/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942_install_debug
    binary="$build/bin/cp2k.pdbg"
    cpu_set=176-191
    export LSAN_OPTIONS=detect_leaks=0
    ;;
  *)
    printf 'unknown build tag: %s\n' "$tag" >&2
    exit 2
    ;;
esac

root="$root_base/positive-smokes-$tag"
export PATH="$env_prefix/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$build/src:$provider/lib:$env_prefix/lib:$dep/lib"
export CP2K_DATA_DIR="$source/data"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export GOTO_NUM_THREADS=1
unset CP2K_GXTB_PARTIAL_QUALIFY_INJECT

mkdir -p "$root"
printf 'case\tmode\tnp\texit_code\ttimed_out\tprogram_ended\tforward_qual\treverse_qual\tpass\n' \
  > "$root/status.tsv"

run_case() {
  name=$1
  mode=$2
  np=$3
  input=$4
  case_dir="$root/$name"
  mkdir -p "$case_dir"
  cp "$inputs/$input" "$case_dir/input.inp"
  if test "$mode" = partial; then
    export CP2K_GXTB_EXCHANGE_STREAM_MODE=KGROUP_PARTIAL_ROOT
    export CP2K_GXTB_EXCHANGE_GRADIENT_MODE=QUALIFY
    export CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION=1
    export CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE=3
    export CP2K_GXTB_EXCHANGE_TRANSFORM_MODE=SEPARABLE
  else
    export CP2K_GXTB_EXCHANGE_STREAM_MODE=LEGACY
    export CP2K_GXTB_EXCHANGE_GRADIENT_MODE=DENSE
    export CP2K_GXTB_EXCHANGE_TRANSFORM_MODE=DENSE
    unset CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION
    unset CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE
  fi
  (
    cd "$case_dir" || exit 99
    timeout --signal=TERM --kill-after=10s 600s taskset -c "$cpu_set" \
      mpiexec -n "$np" "$binary" -i input.inp > run.log 2>&1
  )
  rc=$?
  timed_out=0
  if test "$rc" -eq 124 -o "$rc" -eq 137; then timed_out=1; fi
  ended=0
  forward=0
  reverse=0
  if grep -aFq 'PROGRAM ENDED' "$case_dir/run.log"; then ended=1; fi
  if grep -aFq 'GXTB-QUALIFICATION_ONLY KGROUP-PARTIAL-ROOT iter=' "$case_dir/run.log"; then forward=1; fi
  if grep -aFq 'GXTB-QUALIFICATION_ONLY KGROUP-PARTIAL-ROOT-REVERSE' "$case_dir/run.log"; then reverse=1; fi
  pass=0
  if test "$rc" -eq 0 -a "$timed_out" -eq 0 -a "$ended" -eq 1; then
    if test "$mode" = dense; then
      pass=1
    elif test "$input" = ch4_gamma_explicit.inp -a "$forward" -eq 1 -a "$reverse" -eq 0; then
      pass=1
    elif test "$forward" -eq 1 -a "$reverse" -eq 1; then
      pass=1
    fi
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$name" "$mode" "$np" "$rc" "$timed_out" "$ended" "$forward" "$reverse" "$pass" \
    >> "$root/status.tsv"
}

run_case k222_mode6_p1 partial 1 ch4_k222_default.inp
run_case k222_mode6_p2 partial 2 ch4_k222_default.inp
run_case k222_mode6_p4 partial 4 ch4_k222_default.inp
run_case gamma_explicit_mode6_p2 partial 2 ch4_gamma_explicit.inp
run_case gamma_explicit_dense_p2 dense 2 ch4_gamma_explicit.inp
run_case gamma_implicit_dense_p2 dense 2 ch4_gamma_implicit.inp

python3 "$root_base/compare.py" "$root/k222_mode6_p2/run.log" "$root/k222_mode6_p1/run.log" \
  > "$root/k222_p2_vs_p1.tsv"
python3 "$root_base/compare.py" "$root/k222_mode6_p4/run.log" "$root/k222_mode6_p1/run.log" \
  > "$root/k222_p4_vs_p1.tsv"
python3 "$root_base/compare.py" "$root/gamma_explicit_mode6_p2/run.log" \
  "$root/gamma_explicit_dense_p2/run.log" > "$root/gamma_mode6_vs_explicit_dense.tsv"
python3 "$root_base/compare.py" "$root/gamma_explicit_mode6_p2/run.log" \
  "$root/gamma_implicit_dense_p2/run.log" > "$root/gamma_mode6_vs_implicit_dense.tsv"

source_line=$(grep -n 'CALL tb_gxtb_kpoint_exchange_gradient_dense(qs_env, tb, gradient, sigma)' \
  "$source/src/tblite_interface.F" | tail -n 1)
reverse_count=$(grep -ac 'GXTB KGROUP-PARTIAL-ROOT-REVERSE' \
  "$root/gamma_explicit_mode6_p2/run.log")
force_count=$(grep -ac '^ FORCES|' "$root/gamma_explicit_mode6_p2/run.log")
stress_count=$(grep -ac '^ STRESS| Analytical stress tensor' \
  "$root/gamma_explicit_mode6_p2/run.log")
printf 'source_sha256\tsource_fallback_line\tpartial_reverse_markers\tforce_rows\tstress_blocks\n' \
  > "$root/gamma_reverse_dense_fallback.tsv"
printf '%s\t%s\t%s\t%s\t%s\n' \
  "$(sha256sum "$source/src/tblite_interface.F" | awk '{print $1}')" \
  "$source_line" "$reverse_count" "$force_count" "$stress_count" \
  >> "$root/gamma_reverse_dense_fallback.tsv"

sha256sum "$binary" "$provider/lib/libtblite.a" "$source/src/tblite_interface.F" \
  "$inputs"/*.inp "$root"/*/run.log "$root"/*.tsv > "$root/SHA256SUMS"
awk -F '\t' 'NR > 1 && $9 != 1 {failed = 1} END {exit failed}' "$root/status.tsv"
