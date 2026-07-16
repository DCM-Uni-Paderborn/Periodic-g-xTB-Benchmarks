#!/usr/bin/env bash
set -u

root=$(cd "$(dirname "$0")" && pwd)
cp2k_root=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714
cp2k=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714/build/cp2k/bin/cp2k.psmp
mpirun=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env/bin/mpirun

export CP2K_DATA_DIR="$cp2k_root/data"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

mkdir -p "$root/runs"
: > "$root/launched.tsv"

for src in "$root"/inputs/*.inp; do
  id=$(basename "$src" .inp)
  run="$root/runs/$id"
  mkdir -p "$run"
  cp "$src" "$run/input.inp"
  sha256sum "$run/input.inp" "$cp2k" > "$run/SHA256SUMS.initial"
  (
    cd "$run" || exit 97
    "$mpirun" --bind-to none -np 4 "$cp2k" -i input.inp -o cp2k.out > launcher.log 2>&1
    rc=$?
    printf '%s\n' "$rc" > returncode.txt
    sha256sum input.inp cp2k.out launcher.log returncode.txt > SHA256SUMS.final
    exit "$rc"
  ) &
  printf '%s\t%s\n' "$id" "$!" | tee -a "$root/launched.tsv"
done

wait
date -u +%FT%TZ > "$root/completed_utc.txt"
