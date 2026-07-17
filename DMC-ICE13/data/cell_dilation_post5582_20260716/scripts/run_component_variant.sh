#!/usr/bin/env bash
set -uo pipefail

if [ "$#" -ne 3 ]; then
  echo "usage: $0 VARIANT CPU_START ROOT" >&2
  exit 2
fi

variant=$1
cpu_start=$2
root=$3
case "$variant" in
  no_exchange) parameter_name=gxtb_no_exchange.toml ;;
  frozen_qvszp) parameter_name=gxtb_frozen_qvszp.toml ;;
  no_anisotropic_multipole) parameter_name=gxtb_no_anisotropic_multipole.toml ;;
  no_acp) parameter_name=gxtb_no_acp.toml ;;
  no_exchange_no_acp) parameter_name=gxtb_no_exchange_no_acp.toml ;;
  *) echo "unsupported variant: $variant" >&2; exit 2 ;;
esac
cp2k=/home/kuehne88/work/codex-gxtb-post5582-clean-20260714/build/cp2k/bin/cp2k.psmp
mpirun=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env/bin/mpirun
campaign="$root/component_campaigns/$variant"
inputs="$campaign/inputs"
runs="$campaign/runs"
export OMP_NUM_THREADS=1
mkdir -p "$runs"
: > "$campaign/launched.tsv"
printf 'case\tpid\tcpus\tstarted_utc\n' >> "$campaign/launched.tsv"

pids=()
index=0
for input in "$inputs"/*.inp; do
  case_name=$(basename "$input" .inp)
  run="$runs/$case_name"
  mkdir -p "$run"
  cp "$input" "$run/input.inp"
  cp "$root/component_variants/$parameter_name" "$run/parameter.toml"
  first=$((cpu_start + 4 * index))
  last=$((first + 3))
  cpus="$first-$last"
  (
    cd "$run" || exit 98
    {
      sha256sum input.inp
      sha256sum "$cp2k"
      param=$(awk '$1 == "PARAM" {print $2; exit}' input.inp)
      sha256sum "$param"
    } > SHA256SUMS.initial
    date -u +%Y-%m-%dT%H:%M:%SZ > started_utc.txt
    taskset -c "$cpus" "$mpirun" --bind-to none -np 4 "$cp2k" -i input.inp -o cp2k.out > launcher.log 2>&1
    rc=$?
    printf '%s\n' "$rc" > returncode.txt
    date -u +%Y-%m-%dT%H:%M:%SZ > completed_utc.txt
    find . -maxdepth 1 -type f ! -name SHA256SUMS.final -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS.final
    exit "$rc"
  ) &
  pid=$!
  pids+=("$pid")
  printf '%s\t%s\t%s\t%s\n' "$case_name" "$pid" "$cpus" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$campaign/launched.tsv"
  index=$((index + 1))
done

failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    failed=$((failed + 1))
  fi
done

date -u +%Y-%m-%dT%H:%M:%SZ > "$campaign/completed_utc.txt"
printf 'completed=%s failed=%s\n' "${#pids[@]}" "$failed" > "$campaign/campaign.log"
(
  cd "$campaign" || exit 98
  find . -type f ! -name SHA256SUMS -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS
)
exit "$failed"
