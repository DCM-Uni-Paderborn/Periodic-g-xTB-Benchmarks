#!/usr/bin/env bash
set -euo pipefail

campaign=${1:?usage: verify_after_runs.sh CAMPAIGN EXPECTED_BINARY_SHA256}
expected_binary=${2:?usage: verify_after_runs.sh CAMPAIGN EXPECTED_BINARY_SHA256}
controller_pid=$(cat "$campaign/controller.pid")

while kill -0 "$controller_pid" 2>/dev/null; do
  sleep 30
done

python3 "$campaign/verify_current_binary_requalification.py" \
  "$campaign" \
  --expected-binary-sha256 "$expected_binary" \
  --output "$campaign/verification.json" \
  >"$campaign/verification.stdout.json"
