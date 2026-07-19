#!/usr/bin/env bash
set -euo pipefail

work=/home/kuehne88/work/dmc13-k222-viii-component-ablation-20260719
binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
launcher=/home/kuehne88/work/gxtb-native-bvk-20260718/launch_pinned_cp2k.sh
expected_binary=b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f
modes=(no_exchange no_acp no_exchange_no_acp)
routes=(full reduced)
cpus=(68 69 70 71 72 73)

test "$(sha256sum "$binary" | awk '{print $1}')" = "$expected_binary"
mkdir -p "$work/results" "$work/validation"

pids=()
index=0
for mode in "${modes[@]}"; do
  for route in "${routes[@]}"; do
    input=$work/inputs/${mode}_${route}.inp
    result=$work/results/$mode/$route
    mkdir -p "$result"
    "$launcher" "dmc-viii-$mode-$route" "${cpus[$index]}" "$binary" \
      "$input" "$result" > "$result/launcher.log" 2>&1 &
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

python3 - "$work" <<'PY'
import hashlib
import json
import re
import sys
from pathlib import Path

work = Path(sys.argv[1])
modes = ("no_exchange", "no_acp", "no_exchange_no_acp")
labels = (
    "Core Hamiltonian energy",
    "Repulsive potential energy",
    "Electrostatic energy",
    "Self-consistent dispersion energy",
    "Non-self consistent dispersion energy",
    "Correction for halogen bonding",
)

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def parse(mode, route):
    run = work / "results" / mode / route
    output = run / "cp2k.out"
    text = output.read_text(errors="replace")
    assert (run / "exit_status").read_text().strip() == "0"
    assert "PROGRAM ENDED AT" in text
    energy = re.search(r"^\s*ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$", text, re.M)
    assert energy
    components = {}
    for label in labels:
        match = re.search(rf"^\s*{re.escape(label)}:\s*([-+0-9.Ee]+)\s*$", text, re.M)
        assert match, label
        components[label] = float(match.group(1))
    return {
        "energy_hartree": float(energy.group(1)),
        "components_hartree": components,
        "input_sha256": digest(work / "inputs" / f"{mode}_{route}.inp"),
        "output_sha256": digest(output),
        "binary_sha256": (run / "binary.sha256").read_text().split()[0],
        "affinity_proof": (run / "affinity_preexec.txt").read_text().strip(),
    }

rows = []
for mode in modes:
    full = parse(mode, "full")
    reduced = parse(mode, "reduced")
    delta = full["energy_hartree"] - reduced["energy_hartree"]
    component_delta = {
        label: full["components_hartree"][label] - reduced["components_hartree"][label]
        for label in labels
    }
    rows.append({
        "mode": mode,
        "full": full,
        "reduced": reduced,
        "full_minus_reduced_energy_hartree": delta,
        "full_minus_reduced_components_hartree": component_delta,
        "parity_within_5e_12_hartree": abs(delta) <= 5.0e-12,
    })

payload = {
    "status": "PASS",
    "phase": "VIII",
    "mesh": "2x2x2",
    "eps_scf": 1.0e-12,
    "rows": rows,
}
target = work / "validation" / "comparison.json"
target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps(payload, indent=2, sort_keys=True))
PY

(
  cd "$work"
  find inputs results validation -type f -print0 | sort -z | \
    xargs -0 sha256sum > SHA256SUMS
  sha256sum -c SHA256SUMS
)
