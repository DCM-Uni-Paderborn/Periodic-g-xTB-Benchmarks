#!/usr/bin/env bash
set -euo pipefail

root=/Users/tkuehne/Documents/g-xTB/mstore-pbc-component-comparison-20260720
package=/Users/tkuehne/.cache/gxtb-part-i-working/DMC-ICE13/reproduction/seidler_dmc13_recalculation
mstore=/Users/tkuehne/Documents/g-xTB/save_tblite_mstore_inorganic_audit_20260719/build-release-v4/app/tblite
pbc=/Users/tkuehne/Documents/g-xTB/save_tblite_pbc_source_baseline/build-tests-packages/app/tblite
mstore_sha=8df9fcc990f15600f0b99316602d1d6adfad43f85a2b0203ae14aad44ad4b1aa
pbc_sha=81f1d9690ff040836c2f40cfe0eaf6aa33822681ec029479c5633785537d1aee

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

test "$(shasum -a 256 "$mstore" | awk '{print $1}')" = "$mstore_sha"
test "$(shasum -a 256 "$pbc" | awk '{print $1}')" = "$pbc_sha"

providers=(mstore pbc)
modes=(full no_exchange no_acp no_exchange_no_acp)
phases=(Ih VII)
for provider in "${providers[@]}"; do
  executable=${!provider}
  for mode in "${modes[@]}"; do
    for phase in "${phases[@]}"; do
      run="$root/results/$provider/$mode/$phase"
      mkdir -p "$run"
      structure="$package/raw/mstore_inorganic_cli/k222/$phase/POSCAR"
      if [[ -f "$run/exit_status" ]] && [[ "$(tr -d '[:space:]' < "$run/exit_status")" == 0 ]] && [[ -s "$run/result.json" ]]; then
        continue
      fi
      command=("$executable" run --method gxtb --acc 0.1 --iterations 300 --no-restart --json result.json)
      if [[ "$mode" != full ]]; then
        command+=(--param "$root/parameters/gxtb_${mode}.toml")
      fi
      command+=("$structure")
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
      if [[ "$status" != 0 ]]; then
        exit "$status"
      fi
    done
  done
done

python3 "$root/evaluate_component_matrix.py"
printf 'complete\n' > "$root/controller.status"
