#!/usr/bin/env bash
set -uo pipefail

root=/home/kuehne88/work/gxtb-dmc-cell-dilation-post5582-20260716T0837Z
cp2k=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714/build/cp2k/bin/cp2k.psmp
mpirun=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env/bin/mpirun
case_name=dmc_Ih_GXTB_pbc_s3_k222
variants=(no_exchange frozen_qvszp no_anisotropic_multipole)

mkdir -p "$root/component_smokes"
: > "$root/component_smokes/smoke_summary.tsv"
printf 'variant\treturncode\tprogram_ended\ttotal_energy_eh\n' >> "$root/component_smokes/smoke_summary.tsv"

for variant in "${variants[@]}"; do
  run="$root/component_smokes/$variant"
  mkdir -p "$run"
  cp "$root/component_campaigns/$variant/inputs/$case_name.inp" "$run/input.inp"
  (
    cd "$run" || exit 98
    date -u +%Y-%m-%dT%H:%M:%SZ > started_utc.txt
    param=$(awk '$1 == "PARAM" {print $2; exit}' input.inp)
    {
      sha256sum input.inp
      sha256sum "$param"
      sha256sum "$cp2k"
    } > SHA256SUMS.initial
    OMP_NUM_THREADS=1 taskset -c 220-223 "$mpirun" --bind-to none -np 4 "$cp2k" -i input.inp -o cp2k.out > launcher.log 2>&1
    rc=$?
    printf '%s\n' "$rc" > returncode.txt
    date -u +%Y-%m-%dT%H:%M:%SZ > completed_utc.txt
    find . -maxdepth 1 -type f ! -name SHA256SUMS.final -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS.final
    energy=$(awk '/ENERGY\| Total FORCE_EVAL/ {value=$NF} END {print value}' cp2k.out)
    ended=$(grep -c 'PROGRAM ENDED AT' cp2k.out || true)
    printf '%s\t%s\t%s\t%s\n' "$variant" "$rc" "$ended" "$energy" >> "$root/component_smokes/smoke_summary.tsv"
  )
done

(
  cd "$root/component_smokes" || exit 98
  find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS
)
