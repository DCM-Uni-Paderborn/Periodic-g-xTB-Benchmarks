#!/usr/bin/env python3
"""Compare printed CP2K forces and analytical stress tensors."""

from __future__ import annotations

import argparse
import pathlib
import re

import numpy as np


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def parse(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, float]:
    lines = path.read_text().splitlines()
    energy = float([line.split()[-1] for line in lines if line.startswith(" ENERGY|")][-1])
    forces = []
    for line in lines:
        match = re.match(rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$", line)
        if match:
            forces.append(tuple(map(float, match.groups())))
    marker = max(i for i, line in enumerate(lines) if line == " STRESS| Analytical stress tensor [bar]")
    stress = []
    for line in lines[marker + 2 : marker + 5]:
        fields = line.split()
        stress.append(tuple(map(float, fields[-3:])))
    if not forces or len(stress) != 3:
        raise ValueError(f"missing observables in {path}")
    return np.array(forces), np.array(stress), energy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("restart", type=pathlib.Path)
    parser.add_argument("cold", type=pathlib.Path)
    args = parser.parse_args()
    force_restart, stress_restart, energy_restart = parse(args.restart)
    force_cold, stress_cold, energy_cold = parse(args.cold)
    print(f"restart={args.restart}")
    print(f"cold={args.cold}")
    print(f"energy_restart_hartree={energy_restart:.15f}")
    print(f"energy_cold_hartree={energy_cold:.15f}")
    print(f"energy_delta_hartree={energy_restart - energy_cold:.17e}")
    print(f"force_restart_hartree_per_bohr={force_restart.tolist()}")
    print(f"force_cold_hartree_per_bohr={force_cold.tolist()}")
    print(f"force_max_abs_delta_hartree_per_bohr={np.max(np.abs(force_restart - force_cold)):.17e}")
    print(f"stress_restart_bar={stress_restart.tolist()}")
    print(f"stress_cold_bar={stress_cold.tolist()}")
    print(f"stress_max_abs_delta_bar={np.max(np.abs(stress_restart - stress_cold)):.17e}")


if __name__ == "__main__":
    main()
