#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")" && pwd)
cd "$root"

sha256sum -c SHA256SUMS

for log in local/debug-ctest.log local/release-ctest.log linux/ctest.log; do
  grep -q '100% tests passed, 0 tests failed out of 2' "$log"
  grep -q 'tblite/exchange.*Passed' "$log"
  grep -q 'tblite/gxtb.*Passed' "$log"
  ! grep -qE '[[]FAILED[]]|tests failed, [1-9]' "$log"
done

test "$(tr -d '[:space:]' < linux/exit-code.txt)" = 0
grep -q '^Cpus_allowed_list:[[:space:]]*141$' linux/affinity-preexec.txt
grep -q '^Cpus_allowed_list:[[:space:]]*141$' linux/evidence-affinity.txt
grep -q '^computed_margin_kib=397339460$' linux/budget.txt
grep -q '^minimum_margin_kib=134217728$' linux/budget.txt
grep -q '^evidence complete$' linux/evidence-complete.txt

patch=source/0001-Qualify-exact-BvK-cache-identity.patch
for field in nao nsh maxsh nsh_id nao_sh ish_at iao_sh frscale omega lrscale \
  ondiag_scale hubbard_exp hubbard_exp_r0 gexp corr_exp hubbard onecxints \
  offdiag_scale rad kq corr_scale corr_rad; do
  grep -q "require_bvk_model_mismatch(exchange, cache, \"$field\"" "$patch"
done
grep -q 'BvK plan accepted a changed representative order' "$patch"
grep -q 'Signed-off-by: Thomas D. Kühne <tkuehne@cp2k.org>' "$patch"

grep -q '^4c4cc95530e6db966118e5abea5dff8870e84ce213cb6fbd1dd5cf5f43a539b1 ' \
  linux/final.sha256
grep -q '^d8805e2d29d2f4858e77c20d4b29099aaa5f3fb9771f610ff0a63c369d2f36e5 ' \
  linux/source-before-build.sha256

printf 'Exact BvK cache-identity archive verification passed\n'
