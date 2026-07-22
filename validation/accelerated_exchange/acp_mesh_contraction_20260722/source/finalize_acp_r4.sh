#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-acp-stream-20260722
envroot=/home/kuehne88/work/gxtb-runtime-part-II-20260717/env
cpu=141
candidate_peak_kib=$((100*1024*1024))
native_peak_kib=$((100*1024*1024))
minimum_margin_kib=$((128*1024*1024))
cp2k_build="$root/cp2k-build-r3"
provider_install="$root/save_tblite-install-r3"
results="$root/results/final-r4"
provenance="$root/provenance/final-r4"
log="$root/logs/final-r4.log"
exit_file="$root/logs/final-r4.exit"

mkdir -p "$results" "$provenance"
trap 'rc=$?; printf "%s\n" "$rc" > "$exit_file"' EXIT

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

printf 'pid=%s expected_cpu=%s\n' "$$" "$cpu" > "$provenance/affinity-preexec.txt"
grep -E '^(Name|State|Cpus_allowed|Cpus_allowed_list):' /proc/$$/status \
  >> "$provenance/affinity-preexec.txt"
if [[ $(awk '/^Cpus_allowed_list:/ {print $2}' /proc/$$/status) != "$cpu" ]]; then
  printf 'final launcher is not bound to singleton CPU %s\n' "$cpu" >&2
  exit 74
fi

capture_budget() {
  local label=$1
  local mem_kib remaining_kib live_cp2k margin_kib pid state rss allowance

  grep -E '^(MemAvailable|MemTotal):' /proc/meminfo \
    > "$provenance/${label}-memory.txt"
  ps -e -o pid=,ppid=,sid=,psr=,rss=,vsz=,nlwp=,etimes=,stat=,comm=,args= --sort=-rss \
    > "$provenance/${label}-all-rss.tsv"
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
  margin_kib=$((mem_kib-remaining_kib-candidate_peak_kib))
  {
    printf 'mem_available_kib=%s\n' "$mem_kib"
    printf 'live_cp2k=%s\n' "$live_cp2k"
    printf 'remaining_growth_allowance_kib=%s\n' "$remaining_kib"
    printf 'candidate_peak_kib=%s\n' "$candidate_peak_kib"
    printf 'computed_margin_kib=%s\n' "$margin_kib"
    printf 'minimum_margin_kib=%s\n' "$minimum_margin_kib"
  } > "$provenance/${label}-budget.txt"
  if ((margin_kib < minimum_margin_kib)); then
    printf '%s launch margin is below 128 GiB\n' "$label" >&2
    exit 75
  fi
}

capture_budget rebuild
git -C "$root/save_tblite-src" diff --check
git -C "$root/cp2k" diff --check
"$envroot/bin/cmake" --build "$cp2k_build" --target cp2k.psmp --parallel 1 \
  > "$log" 2>&1

bin="$cp2k_build/bin/cp2k.psmp"
sha256sum "$bin" "$cp2k_build/src/libcp2k.so.2026.2" \
  "$provider_install/lib/libtblite.a" "$provenance/provider.patch" \
  "$provenance/cp2k.patch" > "$provenance/binary-and-source.sha256"
ldd "$bin" > "$provenance/cp2k-ldd.txt"
git -C "$root/save_tblite-src" status --short > "$provenance/provider-status.txt"
git -C "$root/cp2k" status --short > "$provenance/cp2k-status.txt"

run_case() {
  local name=$1 mode=$2 source_input=$3
  local dir="$results/$name"

  capture_budget "$name"
  mkdir -p "$dir"
  cp "$source_input" "$dir/input.inp"
  sha256sum "$dir/input.inp" > "$dir/input.sha256"
  (
    cd "$dir"
    env CP2K_GXTB_ACP_MESH_CONTRACTION="$mode" \
      "$bin" -i input.inp -o cp2k.out
  )
  grep -q 'PROGRAM ENDED AT' "$dir/cp2k.out"
  sha256sum "$dir/cp2k.out" > "$dir/output.sha256"
}

ch4="$root/cp2k/tests/xTB/regtest-tblite-gxtb/CH4_gxtb_kp_acceleration_production.inp"
ch4_3="$root/cp2k/tests/xTB/regtest-tblite-gxtb/CH4_gxtb_kp_symmetry_fused_production.inp"
si="$root/cp2k/tests/xTB/regtest-tblite-gxtb/Si_prim_gxtb_kp_shifted_full.inp"
o2="$root/cp2k/tests/xTB/regtest-tblite-gxtb/O2_gxtb_uks_kp_311.inp"

for mode in DENSE STREAMED QUALIFY; do
  mode_lc=${mode,,}
  run_case "ch4-2x2x2-$mode_lc" "$mode" "$ch4"
  run_case "ch4-3x3x3-$mode_lc" "$mode" "$ch4_3"
  run_case "si-shifted-$mode_lc" "$mode" "$si"
  run_case "o2-uks-$mode_lc" "$mode" "$o2"
done

capture_budget invalid-selector
invalid_dir="$results/invalid-selector"
mkdir -p "$invalid_dir"
cp "$ch4" "$invalid_dir/input.inp"
set +e
(
  cd "$invalid_dir"
  env CP2K_GXTB_ACP_MESH_CONTRACTION=BOGUS \
    "$bin" -i input.inp -o cp2k.out
)
invalid_rc=$?
set -e
printf '%s\n' "$invalid_rc" > "$invalid_dir/exit-code.txt"
if ((invalid_rc == 0)); then
  printf 'invalid ACP selector unexpectedly succeeded\n' >&2
  exit 76
fi
grep -q 'Unknown CP2K_GXTB_ACP_MESH_CONTRACTION value' "$invalid_dir/cp2k.out"

for case_name in ch4-2x2x2-streamed ch4-3x3x3-streamed si-shifted-streamed; do
  grep -q 'GXTB-ACP-MESH STREAMED nFull=' "$results/$case_name/cp2k.out"
  grep -q 'GXTB-ACP-MESH SPARSE-REVERSE projectorImages=' "$results/$case_name/cp2k.out"
done
for case_name in ch4-2x2x2-qualify ch4-3x3x3-qualify si-shifted-qualify; do
  grep -q 'ACP_SPARSE_RESPONSE_QUALIFY residual=' "$results/$case_name/cp2k.out"
  grep -q 'GXTB-QUALIFICATION_ONLY ACP-SPARSE-REVERSE' "$results/$case_name/cp2k.out"
done
for case_name in o2-uks-dense o2-uks-streamed o2-uks-qualify; do
  grep -q '# Total charge and spin         7.000000     5.000000' \
    "$results/$case_name/cp2k.out"
done

find "$results" -type f -print0 | sort -z | xargs -0 sha256sum \
  > "$provenance/results.sha256"
