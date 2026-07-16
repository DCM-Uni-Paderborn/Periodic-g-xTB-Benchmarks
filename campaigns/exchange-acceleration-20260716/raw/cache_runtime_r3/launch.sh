#!/usr/bin/env bash
set -u

root=/home/kuehne88/work/codex-gxtb-exchange-cache-cp2k-runtime-r3-20260716
binary=/home/kuehne88/work/codex-gxtb-exchange-cache-cp2k-20260716/build/bin/cp2k.psmp
data=/home/kuehne88/work/codex-gxtb-exchange-cache-cp2k-20260716/source/data
inputs="$root/inputs"
runs="$root/runs"
mkdir -p "$runs"

threads=8
index=0
for input in "$inputs"/*.inp; do
  name=${input##*/}
  name=${name%.inp}
  run="$runs/$name"
  first=$((index * threads))
  last=$((first + threads - 1))
  index=$((index + 1))
  mkdir -p "$run"

  if [[ -f "$run/returncode.txt" ]]; then
    continue
  fi
  if [[ -f "$run/launcher.pid" ]] && kill -0 "$(<"$run/launcher.pid")" 2>/dev/null; then
    continue
  fi

  qualifier=
  if [[ "$name" == ch4_spglib_oracle_iter2 ]]; then
    qualifier=2
  fi

  cp "$input" "$run/input.inp"
  sha256sum "$run/input.inp" "$binary" > "$run/SHA256SUMS.initial"
  nohup bash -c '
    run=$1
    binary=$2
    data=$3
    first=$4
    last=$5
    qualifier=$6
    cd "$run" || exit 97
    export OMP_NUM_THREADS=8
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export BLIS_NUM_THREADS=1
    export CP2K_DATA_DIR=$data
    if [[ -n "$qualifier" ]]; then
      export CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION=$qualifier
    fi
    start_ns=$(date +%s%N)
    taskset -c "$first-$last" "$binary" -i input.inp -o cp2k.out &
    cp2k_pid=$!
    max_rss_kb=0
    while kill -0 "$cp2k_pid" 2>/dev/null; do
      rss_kb=$(grep "^VmRSS:" "/proc/$cp2k_pid/status" 2>/dev/null | tr -s " " | cut -d " " -f 2)
      if [[ -n "$rss_kb" ]] && ((rss_kb > max_rss_kb)); then
        max_rss_kb=$rss_kb
      fi
      sleep 0.2
    done
    wait "$cp2k_pid"
    rc=$?
    end_ns=$(date +%s%N)
    printf "wall_ns %s\nmax_rss_kb %s\n" "$((end_ns - start_ns))" "$max_rss_kb" > time.txt
    printf "%s\n" "$rc" > returncode.txt
    sha256sum input.inp cp2k.out launcher.log returncode.txt time.txt > SHA256SUMS.final
    exit "$rc"
  ' runtime-job "$run" "$binary" "$data" "$first" "$last" "$qualifier" >"$run/launcher.log" 2>&1 </dev/null &
  printf '%s\n' "$!" > "$run/launcher.pid"
  printf '%s\t%s\t%s-%s\n' "$name" "$!" "$first" "$last" >> "$root/launched.tsv"
done
