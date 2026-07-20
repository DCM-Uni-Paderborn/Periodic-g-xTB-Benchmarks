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

for reservation in "$reservation_root"/*.reservation; do
  [[ -e $reservation ]] || continue
  read -r reserved_pid reserved_cpu < "$reservation" || true
  if [[ -z ${reserved_pid:-} || ! -d /proc/$reserved_pid ]]; then
    rm -f "$reservation"
    continue
  fi
  state=$(awk '/^State:/{print $2}' "/proc/$reserved_pid/status" 2>/dev/null || true)
  if [[ $state == Z ]]; then
    rm -f "$reservation"
    continue
  fi
  if [[ $reserved_cpu == "$cpu" ]]; then
    printf 'CPU %s is reserved by PID %s (%s)\n' "$cpu" "$reserved_pid" "$reservation" >&2
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
taskset -c "$cpu" bash -c '
  set -euo pipefail
  expected=$1
  proof=$2
  result_dir=$3
  shift 3
  allowed=$(awk "/^Cpus_allowed_list:/{print \$2}" "/proc/$$/status")
  {
    printf "pid=%s expected_cpu=%s allowed=%s\\n" "$$" "$expected" "$allowed"
    grep -E "^(Name|State|Cpus_allowed|Cpus_allowed_list):" "/proc/$$/status"
  } > "$proof"
  [[ $allowed == "$expected" ]] || exit 97
  cd "$result_dir"
  exec "$@"
' bash "$cpu" "$proof" "$result_dir" "$binary" -i "$input" -o "$output" &
child=$!
printf '%s %s\n' "$child" "$cpu" > "$reservation"
set +e
wait "$child"
status=$?
set -e
printf '%s\n' "$status" > "$result_dir/exit_status"
exit "$status"
