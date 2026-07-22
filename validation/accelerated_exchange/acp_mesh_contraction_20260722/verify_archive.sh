#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")" && pwd)
cd "$root"

sha256sum -c SHA256SUMS

while IFS= read -r -d '' output; do
  grep -q 'PROGRAM ENDED AT' "$output"
done < <(find local linux/results/final-r4 -type f -name cp2k.out \
  ! -path '*/invalid-selector/*' -print0)

for output in \
  local/ch4-streamed/cp2k.out \
  local/ch4-3x3x3-streamed/cp2k.out \
  local/permanent-symmetry-regression/production/cp2k.out \
  local/si-streamed/cp2k.out \
  linux/results/final-r4/ch4-2x2x2-streamed/cp2k.out \
  linux/results/final-r4/ch4-3x3x3-streamed/cp2k.out \
  linux/results/final-r4/si-shifted-streamed/cp2k.out; do
  grep -q 'GXTB-ACP-MESH STREAMED nFull=' "$output"
  grep -q 'fullStorage=0' "$output"
  grep -q 'GXTB-ACP-MESH SPARSE-REVERSE projectorImages=' "$output"
  grep -q 'fullDifferenceSet=0' "$output"
done

for output in \
  local/ch4-qualify/cp2k.out \
  local/ch4-3x3x3-qualify/cp2k.out \
  local/permanent-symmetry-regression/qualify/cp2k.out \
  local/si-qualify/cp2k.out \
  linux/results/final-r4/ch4-2x2x2-qualify/cp2k.out \
  linux/results/final-r4/ch4-3x3x3-qualify/cp2k.out \
  linux/results/final-r4/si-shifted-qualify/cp2k.out; do
  grep -q 'ACP_SPARSE_RESPONSE_QUALIFY residual=' "$output"
  grep -q 'GXTB-QUALIFICATION_ONLY ACP-SPARSE-REVERSE' "$output"
done

for control in local/invalid-selector linux/results/final-r4/invalid-selector; do
  test "$(tr -d '[:space:]' < "$control/exit-code.txt")" != 0
  grep -q 'Unknown CP2K_GXTB_ACP_MESH_CONTRACTION value' "$control/cp2k.out"
done

printf 'ACP archive verification passed\n'
