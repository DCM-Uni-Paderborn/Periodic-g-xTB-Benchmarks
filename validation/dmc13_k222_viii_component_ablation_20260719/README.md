# DMC-ICE13 phase VIII component ablation

This directory records a same-binary, tight-SCC comparison of explicit full-grid
and SPGLIB-reduced native `2x2x2` calculations for phase VIII. Three model
variants disable (i) nonlocal exchange, (ii) ACP, or (iii) both terms. Every run
uses the CP2K binary identified in `verification.json`, terminates normally, and
has a recorded disjoint singleton CPU affinity.

The initial controller file `validation/comparison.json` used `status: PASS` only
to denote successful execution; its deliberately very tight `5e-12` hartree
diagnostic is not a physical acceptance threshold. The independent verifier
`verify_component_ablation.py` applies both a per-cell and a per-water numerical
equivalence threshold and writes `verification.json`.

The full/reduced residual is present in every ablation and ranges from roughly
`0.08` to `0.25` nanohartree per primitive cell. Its largest value corresponds
to less than `0.000001 kJ mol-1` per water molecule. Consequently, the residual
cannot be assigned uniquely to nonlocal exchange, ACP, or their coupling. It is
a numerically negligible consequence of the different symmetry-dependent
summation order.

Run the independent verification with:

```console
python3 verify_component_ablation.py
```
