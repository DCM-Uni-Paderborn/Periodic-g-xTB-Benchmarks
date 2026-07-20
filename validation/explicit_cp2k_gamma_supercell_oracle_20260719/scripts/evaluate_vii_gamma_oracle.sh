#!/usr/bin/env bash
set -euo pipefail

oracle=$(cd "$(dirname "$0")/.." && pwd)
root=$(cd "$oracle/../.." && pwd)
native="$root/DMC-ICE13/reproduction/seidler_dmc13_recalculation/raw/cp2k_native/k222-reduced/VII"
cli="$root/DMC-ICE13/reproduction/seidler_dmc13_recalculation/raw/current_pbc_cli/cli-k222/VII"
gamma="$oracle/results/VII/gamma_supercell_k222"

python3 "$oracle/scripts/compare_gamma_supercell_oracle.py" \
  "$native/cp2k.out" \
  "$gamma/cp2k.out" \
  "$cli/tblite.json" \
  --replicas 8 \
  --parity-tolerance 2e-7 \
  --require-binary-sha256 \
    b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f \
  --require-cli-binary-sha256 \
    f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a \
  --native-input "$native/input.inp" \
  --gamma-input "$oracle/inputs/VII/input.inp" \
  --cli-input "$cli/POSCAR" \
  --require-native-input-sha256 \
    a3ae9af9154ff0278f0f27482c15c5ca199656ef490a107b505a5ea360e0beac \
  --require-gamma-input-sha256 \
    6481b70e1a437e92247649e06ff90798358f0606dedf60c07b11e4ecab557a21 \
  --require-cli-input-sha256 \
    4de281e3ab3632f443b22f99162e80a3a327c0b8124489c015d4b68fa62d1d91 \
  --output "$oracle/verification-vii.json"
