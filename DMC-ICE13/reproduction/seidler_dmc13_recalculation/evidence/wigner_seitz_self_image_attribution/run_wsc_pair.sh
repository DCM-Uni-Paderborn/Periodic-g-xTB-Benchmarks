#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 LABEL EXECUTABLE" >&2
  exit 2
fi

label=$1
executable=$2
root=/Users/tkuehne/Documents/g-xTB/wsc-index-attribution-20260720
package=/Users/tkuehne/.cache/gxtb-part-i-working/DMC-ICE13/reproduction/seidler_dmc13_recalculation

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

for phase in Ih VII; do
  run="$root/results/$label/$phase"
  structure="$package/raw/mstore_inorganic_cli/k222/$phase/POSCAR"
  mkdir -p "$run"
  command=("$executable" run --method gxtb --acc 0.1 --iterations 300 --no-restart --json result.json "$structure")
  printf '%q ' "${command[@]}" > "$run/command.txt"
  printf '\n' >> "$run/command.txt"
  shasum -a 256 "$executable" > "$run/binary.sha256"
  shasum -a 256 "$structure" > "$run/input.sha256"
  set +e
  (
    cd "$run"
    nice -n 10 "${command[@]}" > process.out 2> process.err
  )
  status=$?
  set -e
  printf '%s\n' "$status" > "$run/exit_status"
  [[ $status -eq 0 ]]
done
