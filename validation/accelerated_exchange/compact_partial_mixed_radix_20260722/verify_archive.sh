#!/usr/bin/env sh
set -eu

cd "$(dirname "$0")"
python3 verify.py
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c SHA256SUMS
else
  shasum -a 256 -c SHA256SUMS
fi
