#!/usr/bin/env bash
set -euo pipefail

work=/home/kuehne88/work/provider-component-attribution-20260719
launcher=/home/kuehne88/work/gxtb-native-bvk-20260718/tools/launch_pinned_command.sh
structure=/home/kuehne88/work/gxtb-native-bvk-20260718/seidler-reproduction/structures/k222/VII/POSCAR
current=/home/kuehne88/work/gxtb-provenance-d933e083-15915c94-20260717-v2/save_tblite-build/app/tblite
pbc=/home/kuehne88/work/gxtb-pbc-author-exact-20260718/build-release/app/tblite
parameter_root=/home/kuehne88/work/gxtb-final-derivative-ablation-20260719/provenance

expected_current=f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a
expected_pbc=692b7da28bdef43cf8795fb33a2b60ccca49a116115715ffda2bda85fefe565d
test "$(sha256sum "$current" | awk '{print $1}')" = "$expected_current"
test "$(sha256sum "$pbc" | awk '{print $1}')" = "$expected_pbc"

mkdir -p "$work/results" "$work/validation"
providers=(current pbc)
modes=(full no_exchange no_acp no_exchange_no_acp)
cpus=(68 69 70 71 72 73 74 75)
pids=()
index=0

for provider in "${providers[@]}"; do
  executable=${!provider}
  for mode in "${modes[@]}"; do
    result="$work/results/$provider/$mode"
    mkdir -p "$result"
    command=("$executable" run --method gxtb --acc 0.1 --iterations 300 --no-restart --json result.json)
    if [[ $mode != full ]]; then
      command+=(--param "$parameter_root/gxtb_${mode}.toml")
    fi
    command+=("$structure")
    "$launcher" "provider-$provider-$mode" "${cpus[$index]}" "$result" -- "${command[@]}" &
    pids+=("$!")
    index=$((index + 1))
  done
done

status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=$?
done
printf '%s\n' "$status" > "$work/controller_exit_status"
test "$status" = 0

python3 - "$work" "$expected_current" "$expected_pbc" <<'PY'
import hashlib
import json
import math
import re
import sys
from pathlib import Path

work = Path(sys.argv[1])
expected = {"current": sys.argv[2], "pbc": sys.argv[3]}
modes = ("full", "no_exchange", "no_acp", "no_exchange_no_acp")
rows = []

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

for mode in modes:
    energies = {}
    records = {}
    for provider in ("current", "pbc"):
        result = work / "results" / provider / mode
        data = json.loads((result / "result.json").read_text())
        output = (result / "process.out").read_text(errors="replace")
        affinity = (result / "affinity_preexec.txt").read_text()
        allowed = re.search(r"allowed=(\d+)", affinity)
        energy = float(data["energy"])
        if not math.isfinite(energy):
            raise AssertionError((provider, mode, energy))
        if int((result / "exit_status").read_text().strip()) != 0:
            raise AssertionError((provider, mode, "exit"))
        energies[provider] = energy
        records[provider] = {
            "energy_hartree_supercell": energy,
            "result_sha256": digest(result / "result.json"),
            "output_sha256": digest(result / "process.out"),
            "singleton_cpu": int(allowed.group(1)) if allowed else None,
            "executable_sha256": expected[provider],
            "scc_converged": "SCC did not converge" not in output,
        }
    delta = energies["pbc"] - energies["current"]
    rows.append({
        "mode": mode,
        "current": records["current"],
        "pbc": records["pbc"],
        "pbc_minus_current_hartree_supercell": delta,
        "pbc_minus_current_hartree_primitive": delta / 8.0,
    })

payload = {
    "status": "PASS",
    "phase": "VII",
    "mesh": "2x2x2 explicit BvK supercell",
    "primitive_repetitions": 8,
    "rows": rows,
}
(work / "validation" / "comparison.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(payload, indent=2, sort_keys=True))
PY

(
  cd "$work"
  find results validation -type f -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
  sha256sum -c SHA256SUMS
)
