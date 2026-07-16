#!/usr/bin/env bash
set -u

root_base=/home/kuehne88/work/gxtb-partial-root-final-35e7942-20260716
root="$root_base/gamma-ultratight-dense-release"
inputs="$root_base/inputs"
source=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_clean_20260716
build=/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_35e7942_release_20260716
provider=/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942_install_release
env_prefix=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
dep=/home/kuehne88/work/codex-gxtb-acpfix-20260714/install/save_tblite
binary="$build/bin/cp2k.psmp"

export PATH="$env_prefix/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$build/src:$provider/lib:$env_prefix/lib:$dep/lib"
export CP2K_DATA_DIR="$source/data"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 GOTO_NUM_THREADS=1
export CP2K_GXTB_EXCHANGE_STREAM_MODE=LEGACY
export CP2K_GXTB_EXCHANGE_GRADIENT_MODE=DENSE
export CP2K_GXTB_EXCHANGE_TRANSFORM_MODE=DENSE
unset CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION
unset CP2K_GXTB_PARTIAL_QUALIFY_INJECT

mkdir -p "$root"
printf 'case\texit_code\ttimed_out\tprogram_ended\tpass\n' > "$root/status.tsv"
for spec in \
  explicit:ch4_gamma_explicit_real_ultratight.inp \
  implicit:ch4_gamma_implicit_ultratight.inp; do
  name=${spec%%:*}
  input=${spec#*:}
  case_dir="$root/$name"
  mkdir -p "$case_dir"
  cp "$inputs/$input" "$case_dir/input.inp"
  (
    cd "$case_dir" || exit 99
    timeout --signal=TERM --kill-after=10s 600s taskset -c 160-175 \
      mpiexec -n 2 "$binary" -i input.inp > run.log 2>&1
  )
  rc=$?
  timed_out=0
  if test "$rc" -eq 124 -o "$rc" -eq 137; then timed_out=1; fi
  ended=0
  if grep -aFq 'PROGRAM ENDED' "$case_dir/run.log"; then ended=1; fi
  pass=0
  if test "$rc" -eq 0 -a "$timed_out" -eq 0 -a "$ended" -eq 1; then pass=1; fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$rc" "$timed_out" "$ended" "$pass" \
    >> "$root/status.tsv"
done

python3 "$root_base/compare.py" "$root/explicit/run.log" "$root/implicit/run.log" \
  > "$root/explicit_dense_vs_implicit_dense.tsv"
sha256sum "$binary" "$provider/lib/libtblite.a" "$source/src/tblite_interface.F" \
  "$inputs/ch4_gamma_explicit_real_ultratight.inp" \
  "$inputs/ch4_gamma_implicit_ultratight.inp" "$root"/*/run.log "$root"/*.tsv \
  > "$root/SHA256SUMS"
awk -F '\t' 'NR > 1 && $5 != 1 {failed = 1} END {exit failed}' "$root/status.tsv"
