# Strict 2x2x2 CP2K-native/direct-CLI parity control

This control repeats ice Ih and ice VII with the exact qualified CP2K and
`save_tblite` executables while tightening both calculation routes from an
effective accuracy of 0.1 to 0.01. The CP2K inputs additionally tighten
`EPS_SCF` from `1.0e-9` to `1.0e-10`. All runs use the same structures as the
qualified 2x2x2 comparison, no restart, and a recorded singleton CPU-affinity
proof.

The direct CLI evaluates the explicit 2x2x2 Born--von Karman supercell; its
total energy is divided by eight before comparison with the CP2K-native Bloch
calculation. Relative energies are normalized by the 12 water molecules in the
primitive cell.

Run

```console
python3 verify_tight_parity.py
```

to reproduce `verification.json`. The verifier checks exit status, executable
and input hashes, printed convergence controls, finite energies, and both
absolute- and relative-energy parity.

## Result

| phase | native minus CLI at accuracy 0.1 (Eh/primitive) | native minus CLI at accuracy 0.01 (Eh/primitive) |
|---|---:|---:|
| Ih | -7.3349611e-9 | -7.2398052e-9 |
| VII | -1.0520830e-7 | -1.0519966e-7 |

For the VII--Ih relative energy, the strict native-minus-CLI difference is
`-2.14328e-5 kJ mol-1 H2O-1`. Tightening the controls therefore does not remove
the already tiny discrepancy: the CP2K-native energies are unchanged at
printed precision, while the direct-CLI supercell energies change by less than
`1e-9 Eh` before normalization. This rules out the SCF stopping threshold and
an energy-unit or supercell-normalization error as the origin of the remaining
numerical difference.
