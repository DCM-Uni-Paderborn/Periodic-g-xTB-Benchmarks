#!/usr/bin/env bash
set -euo pipefail

rank=${OMPI_COMM_WORLD_RANK:?missing OMPI_COMM_WORLD_RANK}
size=${OMPI_COMM_WORLD_SIZE:?missing OMPI_COMM_WORLD_SIZE}
proof_dir=${GXTB_AFFINITY_DIR:?missing GXTB_AFFINITY_DIR}
expected_csv=${GXTB_EXPECTED_CORES:?missing GXTB_EXPECTED_CORES}
expected_ranks=${GXTB_EXPECTED_RANKS:?missing GXTB_EXPECTED_RANKS}

if [[ "$size" != "$expected_ranks" ]]; then
  printf 'rank-count mismatch: OMPI=%s expected=%s\n' "$size" "$expected_ranks" >&2
  exit 91
fi

IFS=',' read -r -a expected <<<"$expected_csv"
if [[ ${#expected[@]} -ne $size ]]; then
  printf 'expected core/rank cardinality mismatch\n' >&2
  exit 92
fi

allowed=$(awk '/^Cpus_allowed_list:/{print $2}' /proc/self/status)
processor=$(awk '{print $39}' "/proc/$$/stat")
expected_core=${expected[$rank]}
if [[ "$allowed" != "$expected_core" ]]; then
  printf 'rank %s has non-singleton/wrong mask %s, expected %s\n' \
    "$rank" "$allowed" "$expected_core" >&2
  exit 93
fi

proof="$proof_dir/preexec_rank_${rank}.tsv"
tmp="$proof.tmp.$$"
printf 'rank\tpid\tprocessor\tallowed\texpected\n%s\t%s\t%s\t%s\t%s\n' \
  "$rank" "$$" "$processor" "$allowed" "$expected_core" >"$tmp"
mv "$tmp" "$proof"

for _ in $(seq 1 500); do
  count=$(find "$proof_dir" -maxdepth 1 -type f -name 'preexec_rank_*.tsv' | wc -l)
  [[ "$count" -eq "$size" ]] && break
  sleep 0.01
done
count=$(find "$proof_dir" -maxdepth 1 -type f -name 'preexec_rank_*.tsv' | wc -l)
if [[ "$count" -ne "$size" ]]; then
  printf 'rank %s timed out waiting for affinity proofs: %s/%s\n' "$rank" "$count" "$size" >&2
  exit 94
fi

mapfile -t observed < <(
  find "$proof_dir" -maxdepth 1 -type f -name 'preexec_rank_*.tsv' -print0 \
    | sort -z \
    | xargs -0 -n1 tail -n1 \
    | awk -F '\t' '{print $4}' \
    | sort -n
)
mapfile -t wanted < <(printf '%s\n' "${expected[@]}" | sort -n)
if [[ ${#observed[@]} -ne "$size" ]] || [[ "${observed[*]}" != "${wanted[*]}" ]]; then
  printf 'rank %s sees duplicate/incomplete masks: observed=%s expected=%s\n' \
    "$rank" "${observed[*]}" "${wanted[*]}" >&2
  exit 95
fi

exec "$@"
