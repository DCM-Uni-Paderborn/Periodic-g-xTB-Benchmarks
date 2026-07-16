#!/usr/bin/env bash
set -u

root=/home/kuehne88/work/codex-gxtb-symmetry-v2-runtime-20260716
binary="$root/build/cp2k/bin/cp2k.psmp"
inputs="$root/inputs"
runs="$root/runs"
data="$root/sources/cp2k/data"
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

  nohup bash -c '
    input=$1
    output=$2
    returncode=$3
    binary=$4
    data=$5
    first=$6
    last=$7
    qualifier=$8
    export OMP_NUM_THREADS=8
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export CP2K_DATA_DIR=$data
    if [[ -n "$qualifier" ]]; then
      export CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION=$qualifier
    fi
    taskset -c "$first-$last" "$binary" -i "$input" -o "$output"
    printf "%s\n" "$?" > "$returncode"
  ' runtime-job "$input" "$run/cp2k.out" "$run/returncode.txt" "$binary" "$data" \
    "$first" "$last" "$qualifier" >"$run/launcher.log" 2>&1 </dev/null &
  printf '%s\n' "$!" > "$run/launcher.pid"
done

