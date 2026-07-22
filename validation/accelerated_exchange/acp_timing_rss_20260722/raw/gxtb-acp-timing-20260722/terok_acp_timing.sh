#!/usr/bin/env bash
set -euo pipefail

source_root=/home/kuehne88/work/gxtb-compact-partial-combination-20260722
run_root=/home/kuehne88/work/gxtb-acp-timing-20260722
cpu=90
candidate_peak_kib=$((100*1024*1024))
live_peak_kib=$((100*1024*1024))
minimum_margin_kib=$((128*1024*1024))
binary="$source_root/cp2k-build-r1/bin/cp2k.psmp"
input="$source_root/inputs/production.inp"
provider="$source_root/provider-install/lib/libtblite.a"

mkdir -p "$run_root/results" "$run_root/provenance"
trap 'rc=$?; printf "%s\n" "$rc" > "$run_root/campaign-exit-status.txt"' EXIT
exec 9>"/tmp/gxtb-cpu-${cpu}.lock"
flock -n 9 || {
  printf 'CPU %s reservation is already held\n' "$cpu" >&2
  exit 70
}
exec 8>"/tmp/gxtb-acp-timing.lock"
flock -n 8 || {
  printf 'ACP timing campaign is already active\n' >&2
  exit 71
}

export OMP_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export BLIS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export LD_LIBRARY_PATH="$source_root/provider-install/lib:$source_root/cp2k-build-r1/lib:${LD_LIBRARY_PATH:-}"

sha256sum "$binary" "$provider" "$input" > "$run_root/provenance/binary-provider-input.sha256"
ldd "$binary" > "$run_root/provenance/cp2k-ldd.txt"
{
  printf 'campaign_start=%s\n' "$(date -Is)"
  printf 'host=%s\n' "$(hostname)"
  printf 'cpu=%s\n' "$cpu"
  env | grep -E '^(OMP|OPENBLAS|MKL|BLIS|GOTO|VECLIB)_' | sort
} > "$run_root/provenance/campaign.txt"
ps -e -o pid=,ppid=,sid=,psr=,rss=,vsz=,nlwp=,etimes=,stat=,user=,comm=,args= --sort=-rss \
  > "$run_root/provenance/pre-campaign-all-rss.tsv"
grep -E '^(MemAvailable|MemTotal):' /proc/meminfo \
  > "$run_root/provenance/pre-campaign-memory.txt"

printf 'sequence\trepetition\tmode\tmeasured\twall_s\tpeak_rss_kib\tenergy_eh\texit_code\n' \
  > "$run_root/summary.tsv"

capture_budget() {
  local label=$1
  local mem_kib remaining_kib live_cp2k margin_kib pid state rss allowance
  local budget_file="$run_root/provenance/${label}-budget.txt"

  grep -E '^(MemAvailable|MemTotal):' /proc/meminfo \
    > "$run_root/provenance/${label}-memory.txt"
  ps -e -o pid=,ppid=,sid=,psr=,rss=,vsz=,nlwp=,etimes=,stat=,user=,comm=,args= --sort=-rss \
    > "$run_root/provenance/${label}-all-rss.tsv"
  mem_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  remaining_kib=0
  live_cp2k=0
  : > "$run_root/provenance/${label}-live-cp2k.tsv"
  while read -r pid; do
    [[ -r /proc/$pid/status ]] || continue
    state=$(awk '/^State:/ {print $2}' /proc/$pid/status)
    rss=$(awk '/^VmRSS:/ {print $2}' /proc/$pid/status)
    [[ $state != Z && ${rss:-0} -gt 0 ]] || continue
    allowance=$((live_peak_kib > rss ? live_peak_kib-rss : 0))
    remaining_kib=$((remaining_kib+allowance))
    live_cp2k=$((live_cp2k+1))
    printf '%s\t%s\t%s\t%s\n' "$pid" "$state" "$rss" "$allowance" \
      >> "$run_root/provenance/${label}-live-cp2k.tsv"
  done < <(pgrep -x cp2k.psmp || true)
  margin_kib=$((mem_kib-remaining_kib-candidate_peak_kib))
  {
    printf 'mem_available_kib=%s\n' "$mem_kib"
    printf 'live_cp2k=%s\n' "$live_cp2k"
    printf 'remaining_growth_allowance_kib=%s\n' "$remaining_kib"
    printf 'candidate_peak_kib=%s\n' "$candidate_peak_kib"
    printf 'computed_margin_kib=%s\n' "$margin_kib"
    printf 'minimum_margin_kib=%s\n' "$minimum_margin_kib"
  } > "$budget_file"
  if ((margin_kib < minimum_margin_kib)); then
    printf '%s launch margin is below 128 GiB\n' "$label" >&2
    exit 72
  fi
}

run_case() {
  local sequence=$1 repetition=$2 mode=$3 measured=$4
  local label dir start_ns end_ns wall_s peak_rss rss state rc energy
  label=$(printf '%02d-r%s-%s' "$sequence" "$repetition" "${mode,,}")
  dir="$run_root/results/$label"
  mkdir -p "$dir"
  cp "$input" "$dir/input.inp"
  capture_budget "$label"

  start_ns=$(date +%s%N)
  setsid taskset -c "$cpu" bash -c '
    set -euo pipefail
    cpu=$1
    dir=$2
    binary=$3
    mode=$4
    printf "pid=%s expected_cpu=%s mode=%s\n" "$$" "$cpu" "$mode" > "$dir/affinity-preexec.txt"
    grep -E "^(Name|State|Cpus_allowed|Cpus_allowed_list):" /proc/$$/status >> "$dir/affinity-preexec.txt"
    test "$(awk "/^Cpus_allowed_list:/ {print \$2}" /proc/$$/status)" = "$cpu"
    exec env CP2K_GXTB_ACP_MESH_CONTRACTION="$mode" \
      "$binary" -i "$dir/input.inp" -o "$dir/cp2k.out"
  ' _ "$cpu" "$dir" "$binary" "$mode" &
  local run_pid=$!

  peak_rss=0
  printf 'elapsed_ms\trss_kib\n' > "$dir/rss-series.tsv"
  while [[ -r /proc/$run_pid/status ]]; do
    state=$(awk '/^State:/ {print $2}' /proc/$run_pid/status 2>/dev/null || true)
    [[ $state != Z ]] || break
    rss=$(ps -o rss= --sid "$run_pid" 2>/dev/null | \
      awk '{sum+=$1} END {print sum+0}' || true)
    rss=${rss:-0}
    ((rss > peak_rss)) && peak_rss=$rss
    printf '%s\t%s\n' "$((($(date +%s%N)-start_ns)/1000000))" "$rss" \
      >> "$dir/rss-series.tsv"
    sleep 0.02
  done
  set +e
  wait "$run_pid"
  rc=$?
  set -e
  end_ns=$(date +%s%N)
  wall_s=$(awk -v ns="$((end_ns-start_ns))" 'BEGIN {printf "%.6f", ns/1.0e9}')
  printf '%s\n' "$rc" > "$dir/exit-status.txt"
  [[ $rc -eq 0 ]]
  grep -q 'PROGRAM ENDED AT' "$dir/cp2k.out"
  if [[ $mode == DENSE ]]; then
    ! grep -q 'GXTB-ACP-MESH STREAMED' "$dir/cp2k.out"
  else
    grep -q 'GXTB-ACP-MESH STREAMED nFull=' "$dir/cp2k.out"
    grep -q 'GXTB-ACP-MESH SPARSE-REVERSE projectorImages=' "$dir/cp2k.out"
  fi
  energy=$(awk '/ENERGY\| Total FORCE_EVAL/ {value=$NF} END {print value}' "$dir/cp2k.out")
  sha256sum "$dir/input.inp" "$dir/cp2k.out" "$dir/affinity-preexec.txt" \
    "$dir/rss-series.tsv" > "$dir/SHA256SUMS"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sequence" "$repetition" "$mode" "$measured" "$wall_s" "$peak_rss" "$energy" "$rc" \
    >> "$run_root/summary.tsv"
}

sequence=0
for mode in DENSE STREAMED; do
  sequence=$((sequence+1))
  run_case "$sequence" 0 "$mode" no
done
for repetition in 1 2 3 4 5; do
  if ((repetition % 2)); then
    modes=(DENSE STREAMED)
  else
    modes=(STREAMED DENSE)
  fi
  for mode in "${modes[@]}"; do
    sequence=$((sequence+1))
    run_case "$sequence" "$repetition" "$mode" yes
  done
done

python3 - "$run_root/summary.tsv" > "$run_root/statistics.tsv" <<'PY'
import csv
import statistics
import sys

with open(sys.argv[1], newline="") as handle:
    rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["measured"] == "yes"]

print("mode\tn\tmedian_wall_s\tmad_wall_s\tmedian_peak_rss_kib\tmad_peak_rss_kib\tenergy_eh")
for mode in ("DENSE", "STREAMED"):
    selected = [row for row in rows if row["mode"] == mode]
    walls = [float(row["wall_s"]) for row in selected]
    rss = [int(row["peak_rss_kib"]) for row in selected]
    median_wall = statistics.median(walls)
    median_rss = statistics.median(rss)
    mad_wall = statistics.median(abs(value-median_wall) for value in walls)
    mad_rss = statistics.median(abs(value-median_rss) for value in rss)
    energies = {row["energy_eh"] for row in selected}
    if len(energies) != 1:
        raise SystemExit(f"non-identical {mode} energies: {sorted(energies)}")
    print(f"{mode}\t{len(selected)}\t{median_wall:.6f}\t{mad_wall:.6f}\t{median_rss:.0f}\t{mad_rss:.0f}\t{next(iter(energies))}")
PY

find "$run_root" -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
  > "$run_root/SHA256SUMS"
printf 'campaign_end=%s\n' "$(date -Is)" >> "$run_root/provenance/campaign.txt"
