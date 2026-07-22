#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-acp-stream-20260722
envroot=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
gate="$root/cache_signature_order_20260722"
build="$root/save_tblite-build-r3"
cpu=141
candidate_kib=$((16*1024*1024))
native_peak_kib=$((100*1024*1024))
minimum_margin_kib=$((128*1024*1024))

mkdir -p "$gate"
trap 'rc=$?; printf "%s\n" "$rc" > "$gate/exit-code.txt"' EXIT

export PATH="$envroot/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$envroot/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export BLIS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

printf 'pid=%s expected_cpu=%s\n' "$$" "$cpu" > "$gate/affinity-preexec.txt"
grep -E '^(Name|State|Cpus_allowed|Cpus_allowed_list):' /proc/$$/status \
  >> "$gate/affinity-preexec.txt"
if [[ $(awk '/^Cpus_allowed_list:/ {print $2}' /proc/$$/status) != "$cpu" ]]; then
  printf 'launcher is not bound to singleton CPU %s\n' "$cpu" >&2
  exit 74
fi

grep -E '^(MemAvailable|MemTotal):' /proc/meminfo > "$gate/memory.txt"
ps -eLo user:20,pid,ppid,lwp,psr,stat,rss,vsz,nlwp,etime,comm,args \
  --sort=-rss > "$gate/all-live-rss.tsv"

mem_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
remaining_kib=0
live_cp2k=0
while read -r pid; do
  [[ -r /proc/$pid/status ]] || continue
  state=$(awk '/^State:/ {print $2}' /proc/$pid/status)
  rss=$(awk '/^VmRSS:/ {print $2}' /proc/$pid/status)
  [[ $state != Z && ${rss:-0} -gt 0 ]] || continue
  allowance=$((native_peak_kib > rss ? native_peak_kib-rss : 0))
  remaining_kib=$((remaining_kib+allowance))
  live_cp2k=$((live_cp2k+1))
done < <(pgrep -x cp2k.psmp || true)
margin_kib=$((mem_kib-remaining_kib-candidate_kib))
{
  printf 'mem_available_kib=%s\n' "$mem_kib"
  printf 'live_cp2k=%s\n' "$live_cp2k"
  printf 'remaining_growth_allowance_kib=%s\n' "$remaining_kib"
  printf 'candidate_peak_kib=%s\n' "$candidate_kib"
  printf 'computed_margin_kib=%s\n' "$margin_kib"
  printf 'minimum_margin_kib=%s\n' "$minimum_margin_kib"
} > "$gate/budget.txt"
if ((margin_kib < minimum_margin_kib)); then
  printf 'launch margin is below 128 GiB\n' >&2
  exit 75
fi

git -C "$root/save_tblite-src" diff --check
sha256sum "$root/save_tblite-src/test/unit/test_exchange.f90" \
  > "$gate/source-before-build.sha256"
"$envroot/bin/cmake" --build "$build" --parallel 1 > "$gate/build.log" 2>&1
"$envroot/bin/ctest" --test-dir "$build" --output-on-failure -V \
  -R '^tblite/(exchange|gxtb)$' > "$gate/ctest.log" 2>&1
sha256sum "$build/test/unit/tblite-tester" \
  "$root/save_tblite-src/test/unit/test_exchange.f90" \
  > "$gate/final.sha256"
