#!/usr/bin/env bash
set -eu

cpu_base=$1
proof_dir=$2
shift 2

rank=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
cpu=$((cpu_base + rank))

exec taskset -c "$cpu" bash -c '
proof=$1
rank=$2
cpu=$3
shift 3
printf "pid=%s rank=%s expected_cpu=%s\n" "$$" "$rank" "$cpu" > "$proof"
grep -E "^(Name|State|Cpus_allowed|Cpus_allowed_list):" /proc/$$/status >> "$proof"
exec "$@"
' run "$proof_dir/affinity-rank${rank}.txt" "$rank" "$cpu" "$@"
