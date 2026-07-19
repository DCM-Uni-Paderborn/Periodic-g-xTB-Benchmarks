#!/usr/bin/env bash
set -euo pipefail

root=${DMC_ROOT:-/home/kuehne88/work/gxtb-native-bvk-20260718}
required_binary=${REQUIRED_BINARY_SHA256:-b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f}
convergence_threshold=${CONVERGENCE_THRESHOLD:-0.10}
maximum_mesh=${MAXIMUM_MESH:-12}
selector_status="$root/status/strict-adaptive-completion.status"
selector_result="$root/status/strict-adaptive-completion/adaptive_endpoints.json"
output_root="$root/status/independent-adaptive-verification"
log="$root/status/independent-adaptive-verification.log"
status_file="$root/status/independent-adaptive-verification.status"

mkdir -p "$output_root" "$root/status"
printf '%s waiting for strict adaptive endpoint selection\n' \
  "$(date --iso-8601=seconds)" >>"$log"
while pgrep -u "$USER" -f 'run_strict_adaptive_completion\.sh' >/dev/null; do
  sleep 30
done

if [[ ! -f "$selector_status" ]] || ! grep -qx 'status=PASS' "$selector_status"; then
  {
    printf 'status=NOT_READY\n'
    printf 'reason=strict adaptive endpoint selection did not pass\n'
  } >"$status_file"
  printf '%s strict adaptive endpoint selection not ready\n' \
    "$(date --iso-8601=seconds)" >>"$log"
  exit 1
fi

python3 "$root/tools/verify_adaptive_dmc13.py" \
  "$root" "$selector_result" \
  "$root/tools/dmc_ice13_relative_energies.csv" \
  --meshes "$(seq -s, 4 "$maximum_mesh")" \
  --threshold "$convergence_threshold" \
  --require-binary-sha256 "$required_binary" \
  --output-json "$output_root/independent-verification.json" \
  >"$output_root/independent-verification.stdout.json" \
  2>"$output_root/independent-verification.stderr"

sha256sum \
  "$root/tools/verify_adaptive_dmc13.py" \
  "$selector_result" \
  "$output_root/independent-verification.json" \
  >"$output_root/SHA256SUMS"
{
  printf 'status=PASS\n'
  printf 'required_binary_sha256=%s\n' "$required_binary"
  printf 'threshold_kj_mol_per_water=%s\n' "$convergence_threshold"
  printf 'maximum_mesh=%s\n' "$maximum_mesh"
} >"$status_file"
printf '%s independent adaptive verification passed\n' \
  "$(date --iso-8601=seconds)" >>"$log"
