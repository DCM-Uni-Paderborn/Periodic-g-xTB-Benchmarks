#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  printf 'usage: %s JOB_ID CPU BINARY INPUT RESULT_DIR\n' "$0" >&2
  exit 64
fi

job_id=$1
cpu=$2
binary=$3
input=$4
result_dir=$5

if [[ ! $job_id =~ ^[A-Za-z0-9._-]+$ ]] || [[ ! $cpu =~ ^[0-9]+$ ]]; then
  printf 'invalid job identifier or CPU: %s %s\n' "$job_id" "$cpu" >&2
  exit 64
fi
if [[ ! -x $binary ]] || [[ ! -r $input ]]; then
  printf 'missing executable or input: %s %s\n' "$binary" "$input" >&2
  exit 66
fi

reservation_root="/tmp/gxtb-cpu-reservations-${USER}"
mkdir -p "$reservation_root" "$result_dir"
sha256sum "$binary" > "$result_dir/binary.sha256"
sha256sum "$input" > "$result_dir/input.sha256"
grep -E '^(MemTotal|MemAvailable):' /proc/meminfo > "$result_dir/memory_prelaunch.txt"
ps -eo pid,psr,rss,etimes,comm,args --sort=-rss > "$result_dir/processes_prelaunch.txt"

exec 9>"$reservation_root/.lock"
flock -x 9

cpu_in_list() {
  local needle=$1 item lo hi
  IFS=',' read -ra items
  for item in "${items[@]}"; do
    if [[ $item == *-* ]]; then
      lo=${item%-*}
      hi=${item#*-}
    else
      lo=$item
      hi=$item
    fi
    if (( needle >= lo && needle <= hi )); then
      return 0
    fi
  done
  return 1
}

for reservation_file in "$reservation_root"/*.reservation; do
  [[ -e $reservation_file ]] || continue
  read -r reserved_pid reserved_cpu < "$reservation_file" || true
  if [[ -z ${reserved_pid:-} || ! -d /proc/$reserved_pid ]]; then
    rm -f "$reservation_file"
    continue
  fi
  state=$(awk '/^State:/{print $2}' "/proc/$reserved_pid/status" 2>/dev/null || true)
  if [[ $state == Z ]]; then
    rm -f "$reservation_file"
    continue
  fi
  if [[ $reserved_cpu == "$cpu" ]]; then
    printf 'CPU %s is reserved by PID %s (%s)\n' "$cpu" "$reserved_pid" "$reservation_file" >&2
    exit 75
  fi
done

for proc in /proc/[0-9]*; do
  pid=${proc##*/}
  [[ $pid == $$ ]] && continue
  [[ $(stat -c %U "$proc" 2>/dev/null || true) == "$USER" ]] || continue
  comm=$(cat "$proc/comm" 2>/dev/null || true)
  case "$comm" in
    cp2k.psmp|cp2k.popt|mpirun|mpiexec|orterun) ;;
    *) continue ;;
  esac
  state=$(awk '/^State:/{print $2}' "$proc/status" 2>/dev/null || true)
  [[ $state == Z ]] && continue
  allowed=$(awk '/^Cpus_allowed_list:/{print $2}' "$proc/status" 2>/dev/null || true)
  if [[ -n $allowed ]] && items=$allowed cpu_in_list "$cpu"; then
    printf 'CPU %s overlaps live %s PID %s with affinity %s\n' "$cpu" "$comm" "$pid" "$allowed" >&2
    exit 75
  fi
done

reservation="$reservation_root/$job_id.reservation"
printf '%s %s\n' "$$" "$cpu" > "$reservation"
flock -u 9
trap 'rm -f "$reservation"' EXIT INT TERM

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

proof="$result_dir/affinity_preexec.txt"
output="$result_dir/cp2k.out"
profile="$result_dir/time_verbose.txt"
taskset -c "$cpu" bash -c '
  set -euo pipefail
  expected=$1
  proof=$2
  result_dir=$3
  profile=$4
  shift 4
  allowed=$(awk "/^Cpus_allowed_list:/{print \$2}" "/proc/$$/status")
  {
    printf "pid=%s expected_cpu=%s allowed=%s\\n" "$$" "$expected" "$allowed"
    grep -E "^(Name|State|Cpus_allowed|Cpus_allowed_list):" "/proc/$$/status"
    env | grep -E "^(OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|BLIS_NUM_THREADS|VECLIB_MAXIMUM_THREADS|NUMEXPR_NUM_THREADS)=" | sort
  } > "$proof"
  [[ $allowed == "$expected" ]] || exit 97
  cd "$result_dir"
  start_ns=$(date +%s%N)
  "$@" &
  calculation_pid=$!
  peak_rss_kib=0
  peak_hwm_kib=0
  printf "elapsed_s rss_kib hwm_kib\n" > rss_samples.tsv
  while kill -0 "$calculation_pid" 2>/dev/null; do
    now_ns=$(date +%s%N)
    rss_kib=$(awk "/^VmRSS:/{print \$2}" "/proc/$calculation_pid/status" 2>/dev/null || true)
    hwm_kib=$(awk "/^VmHWM:/{print \$2}" "/proc/$calculation_pid/status" 2>/dev/null || true)
    rss_kib=${rss_kib:-0}
    hwm_kib=${hwm_kib:-0}
    (( rss_kib > peak_rss_kib )) && peak_rss_kib=$rss_kib
    (( hwm_kib > peak_hwm_kib )) && peak_hwm_kib=$hwm_kib
    elapsed_ns=$((now_ns - start_ns))
    printf "%d.%09d %s %s\n" "$((elapsed_ns / 1000000000))" \
      "$((elapsed_ns % 1000000000))" "$rss_kib" "$hwm_kib" >> rss_samples.tsv
    sleep 0.1
  done
  set +e
  wait "$calculation_pid"
  calculation_status=$?
  set -e
  end_ns=$(date +%s%N)
  elapsed_ns=$((end_ns - start_ns))
  {
    printf "elapsed_seconds=%d.%09d\n" "$((elapsed_ns / 1000000000))" \
      "$((elapsed_ns % 1000000000))"
    printf "peak_sampled_rss_kib=%s\n" "$peak_rss_kib"
    printf "peak_sampled_hwm_kib=%s\n" "$peak_hwm_kib"
    printf "calculation_exit_status=%s\n" "$calculation_status"
  } > "$profile"
  exit "$calculation_status"
' bash "$cpu" "$proof" "$result_dir" "$profile" "$binary" -i "$input" -o "$output" &
child=$!
printf '%s %s\n' "$child" "$cpu" > "$reservation"
set +e
wait "$child"
status=$?
set -e
printf '%s\n' "$status" > "$result_dir/exit_status"
grep -E '^(MemTotal|MemAvailable):' /proc/meminfo > "$result_dir/memory_postrun.txt"
exit "$status"
