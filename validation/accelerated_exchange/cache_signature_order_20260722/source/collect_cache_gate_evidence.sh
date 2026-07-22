#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-acp-stream-20260722
source_root="$root/save_tblite-src"
gate="$root/cache_signature_order_20260722"
cpu=141

export OMP_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export BLIS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

printf 'pid=%s expected_cpu=%s\n' "$$" "$cpu" > "$gate/evidence-affinity.txt"
grep -E '^(Name|State|Cpus_allowed|Cpus_allowed_list):' /proc/$$/status \
  >> "$gate/evidence-affinity.txt"
if [[ $(awk '/^Cpus_allowed_list:/ {print $2}' /proc/$$/status) != "$cpu" ]]; then
  printf 'evidence process is not bound to singleton CPU %s\n' "$cpu" >&2
  exit 74
fi

grep -E '^(MemAvailable|MemTotal):' /proc/meminfo > "$gate/evidence-memory.txt"
ps -eLo user:20,pid,ppid,lwp,psr,stat,rss,vsz,nlwp,etime,comm,args \
  --sort=-rss > "$gate/evidence-all-live-rss.tsv"
git -C "$source_root" status --porcelain=v2 > "$gate/remote-source-status.txt"
git -C "$source_root" log -8 --date=iso-strict --format=fuller \
  --output="$gate/remote-source-history.txt"
git -C "$source_root" diff --binary --output="$gate/remote-source-working-tree.patch"
sha256sum \
  "$source_root/src/tblite/acp/type.f90" \
  "$source_root/src/tblite/cp2k_compat.f90" \
  "$source_root/test/unit/test_exchange.f90" \
  "$source_root/test/unit/test_gxtb.f90" \
  > "$gate/remote-source-files.sha256"
uname -a > "$gate/platform.txt"
printf 'evidence complete\n' > "$gate/evidence-complete.txt"
