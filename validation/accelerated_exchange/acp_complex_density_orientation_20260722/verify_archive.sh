#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c SHA256SUMS
else
  shasum -a 256 -c SHA256SUMS
fi

grep -q 'residual=  3.9903E-05' local/pre_gate/cp2k.out

for output in \
  local/diagnostic/swapped-qualify/cp2k.out \
  local/diagnostic/two-dimensional-qualify/cp2k.out \
  linux/cases/acp-h2-1d-3x1x1-fixed/cp2k.out \
  linux/cases/acp-h2-2d-3x1x3-fixed/cp2k.out; do
  grep -q 'ACP_SPARSE_RESPONSE_QUALIFY residual=  0.0000E+00' "$output"
  grep -q 'GXTB-QUALIFICATION_ONLY ACP-SPARSE-REVERSE' "$output"
  grep -q 'PROGRAM ENDED AT' "$output"
done

grep -q 'residual=  5.5511E-17' local/diagnostic/swapped-qualify/cp2k.out
grep -q 'residual=  6.6613E-16' local/diagnostic/two-dimensional-qualify/cp2k.out
grep -q 'residual=  1.1102E-16' linux/cases/acp-h2-1d-3x1x1-fixed/cp2k.out
grep -q 'residual=  5.5511E-16' linux/cases/acp-h2-2d-3x1x3-fixed/cp2k.out

for case_dir in linux/cases/acp-h2-1d-3x1x1-fixed linux/cases/acp-h2-2d-3x1x3-fixed; do
  grep -q '^df6466552a495e94e710174dbf468ec1765f1a4690ef955f9129ef983f35790b ' \
    "$case_dir/provenance/pre-run-sha256.txt"
  grep -q '^fe210c64a4c4fa6897668a8657dd234046143c55bda2f7c9279d24108f2f152a ' \
    "$case_dir/provenance/pre-run-sha256.txt"
  grep -q $'Cpus_allowed_list:\t90' "$case_dir/provenance/affinity/affinity-rank0.txt"
  grep -q $'Cpus_allowed_list:\t91' "$case_dir/provenance/affinity/affinity-rank1.txt"
done

grep -q '100% tests passed, 0 tests failed out of 2' linux/build/provider-tests.log
grep -q 'Signed-off-by: Thomas D. Kühne <tkuehne@cp2k.org>' source/provider/0001-*.patch
grep -q 'Signed-off-by: Thomas D. Kühne <tkuehne@cp2k.org>' source/provider/0002-*.patch
grep -q 'Signed-off-by: Thomas D. Kühne <tkuehne@cp2k.org>' source/cp2k/0001-*.patch

printf '%s\n' 'ACP complex-density orientation archive: PASS'

