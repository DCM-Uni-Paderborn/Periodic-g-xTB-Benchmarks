# DMC-ICE13 phase XVII derivative component ablation

This same-binary `2x2x2` native-k data set compares explicit full-grid and
SPGLIB-reduced energy, forces, and stress for phase XVII after separately
disabling nonlocal exchange, ACP, and both components. The purpose is to
localize possible symmetry-response defects without changing structures,
mesh, SCC settings, or the executable.

All six calculations terminate normally and carry exact singleton CPU-affinity,
binary, input, and output records. The independent verifier reports zero
difference at printed precision for total energy, every Cartesian force
component, and every stress-tensor component in all three ablations. Together
with the complete-model derivative gates in the adjacent final low-k package,
this excludes a general force/stress failure in the full-to-reduced route.

Run the independent verification with:

```console
python3 verify_derivative_component_ablation.py
```
