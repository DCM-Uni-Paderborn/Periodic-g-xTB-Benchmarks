#!/usr/bin/env bash
set -eo pipefail

CP2K=${CP2K:-/home/kuehne88/bin/cp2k-current-tblite.psmp}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

phases=(Ih II III IV VI VII VIII IX XI XIII XIV XV XVII)
methods=(GFN1 GFN2)

cd "$(dirname "$0")"

for method in "${methods[@]}"; do
  for phase in "${phases[@]}"; do
    run_dir="runs/${method}/${phase}"
    input="ice_${phase}_${method}.inp"
    output="ice_${phase}_${method}.out"
    mkdir -p "${run_dir}"
    cp "inputs/${input}" "${run_dir}/"
    (
      cd "${run_dir}"
      if grep -q "ENERGY| Total FORCE_EVAL" "${output}" 2>/dev/null; then
        echo "SKIP ${method} ${phase}"
      else
        echo "RUN  ${method} ${phase}"
        rm -f "${output}"
        "${CP2K}" -i "${input}" -o "${output}"
      fi
    )
  done
done

python3 dmc_ice13_pipeline.py > results.out
