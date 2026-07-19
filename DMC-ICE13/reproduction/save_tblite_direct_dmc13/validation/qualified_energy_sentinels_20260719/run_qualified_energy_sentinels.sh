#!/usr/bin/env bash
set -euo pipefail

root=/home/kuehne88/work/gxtb-native-bvk-20260718
qualified=/home/kuehne88/work/gxtb-qualified-reference-20260719
baseline_binary=/home/kuehne88/work/gxtb-final-clean-20260718/cp2k-build/bin/cp2k.psmp
qualified_binary="$qualified/cp2k-build/bin/cp2k.psmp"
launcher="$root/launch_pinned_cp2k.sh"
result_root="$root/runs/diagnostics/qualified-energy-sentinels"
log="$root/status/qualified-energy-sentinels.log"

mkdir -p "$root/status" "$result_root"

while [[ ! -f "$qualified/STATUS" ]] || \
      ! grep -q '^status=BUILD_PASS$' "$qualified/STATUS"; do
  sleep 30
done

while (( $(awk '/^MemAvailable:/{printf "%d\n", $2/1024/1024}' /proc/meminfo) < 80 )); do
  sleep 30
done

run_one() {
  local label=$1 cpu=$2 binary=$3 input=$4 result=$5
  if grep -q 'PROGRAM ENDED AT' "$result/cp2k.out" 2>/dev/null; then
    return
  fi
  if [[ -d "$result" ]]; then
    mv "$result" "$result-preserved-$(date +%Y%m%dT%H%M%S)"
  fi
  "$launcher" "$label" "$cpu" "$binary" "$input" "$result" \
    >"$result-launcher.log" 2>&1
  test "$(tr -d '\n' <"$result/exit_status")" = 0
  grep -q 'PROGRAM ENDED AT' "$result/cp2k.out"
}

run_one qualified-k222-Ih 97 "$qualified_binary" \
  "$root/inputs/k222/Ih/input.inp" "$result_root/qualified-k222-Ih"
run_one qualified-k222-VII 99 "$qualified_binary" \
  "$root/inputs/k222/VII/input.inp" "$result_root/qualified-k222-VII"
run_one baseline-k333-VII 97 "$baseline_binary" \
  "$root/inputs/k333-reduced/VII/input.inp" "$result_root/baseline-k333-VII"
run_one qualified-k333-VII 99 "$qualified_binary" \
  "$root/inputs/k333-reduced/VII/input.inp" "$result_root/qualified-k333-VII"

python3 - "$root" "$result_root" <<'PY' >"$result_root/comparison.json"
import json
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
result = Path(sys.argv[2])
pattern = re.compile(r"ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$")

def energy(path: Path) -> float:
    values = [float(match.group(1)) for line in path.read_text().splitlines()
              if (match := pattern.search(line))]
    if not values or "PROGRAM ENDED AT" not in path.read_text():
        raise SystemExit(f"incomplete output: {path}")
    return values[-1]

reference = {
    "k222_Ih": energy(root / "runs/final-clean-k222/Ih/cp2k.out"),
    "k222_VII": energy(root / "runs/final-clean-k222/VII/cp2k.out"),
    "k333_VII": energy(result / "baseline-k333-VII/cp2k.out"),
}
candidate = {
    "k222_Ih": energy(result / "qualified-k222-Ih/cp2k.out"),
    "k222_VII": energy(result / "qualified-k222-VII/cp2k.out"),
    "k333_VII": energy(result / "qualified-k333-VII/cp2k.out"),
}
deltas = {key: candidate[key] - reference[key] for key in reference}
maximum = max(map(abs, deltas.values()))
if maximum > 1.0e-12:
    raise SystemExit(f"qualified energy sentinel mismatch: {maximum:.16e} Ha")
print(json.dumps({
    "reference_hartree": reference,
    "qualified_hartree": candidate,
    "qualified_minus_reference_hartree": deltas,
    "maximum_absolute_delta_hartree": maximum,
    "status": "PASS",
}, indent=2, sort_keys=True))
PY

printf '%s sentinel comparison complete\n' "$(date --iso-8601=seconds)" >>"$log"
