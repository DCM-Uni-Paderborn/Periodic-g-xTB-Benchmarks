#!/usr/bin/env bash
set -uo pipefail
root=/home/kuehne88/work/gxtb-dmc-cell-dilation-post5582-20260716T0837Z
run="$root/component_smokes/no_acp"
cp2k=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714/build/cp2k/bin/cp2k.psmp
mpirun=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env/bin/mpirun
tblite=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714/install/save_tblite/bin/tblite
param=/home/kuehne88/gx5582/gxtb_no_acp.toml
mkdir -p "$run"
"$tblite" param "$param" --output "$root/component_variants/validated_gxtb_no_acp.toml" > "$root/component_variants/gxtb_no_acp.toml.validation.stdout" 2> "$root/component_variants/gxtb_no_acp.toml.validation.stderr"
param_rc=$?
cp "$root/component_campaigns/no_acp/inputs/dmc_Ih_GXTB_pbc_s3_k222.inp" "$run/input.inp"
date -u +%Y-%m-%dT%H:%M:%SZ > "$run/started_utc.txt"
(
  cd "$run" || exit 98
  { sha256sum input.inp; sha256sum "$param"; sha256sum "$cp2k"; } > SHA256SUMS.initial
  OMP_NUM_THREADS=1 taskset -c 220-223 "$mpirun" --bind-to none -np 4 "$cp2k" -i input.inp -o cp2k.out > launcher.log 2>&1
  rc=$?
  printf '%s\n' "$rc" > returncode.txt
  date -u +%Y-%m-%dT%H:%M:%SZ > completed_utc.txt
  find . -maxdepth 1 -type f ! -name SHA256SUMS.final -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS.final
  energy=$(awk '/ENERGY\| Total FORCE_EVAL/ {value=$NF} END {print value}' cp2k.out)
  ended=$(grep -c 'PROGRAM ENDED AT' cp2k.out || true)
  printf 'variant\tparam_rc\treturncode\tprogram_ended\ttotal_energy_eh\nno_acp\t%s\t%s\t%s\t%s\n' "$param_rc" "$rc" "$ended" "$energy" > smoke_summary.tsv
  exit "$rc"
)
