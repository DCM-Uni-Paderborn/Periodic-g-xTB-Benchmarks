# Part-II publication figures

This directory contains publication-ready figures that belong exclusively to
the Part-II implementation, qualification, memory, and scaling study.  It does
not contain the Part-I application benchmarks.

## Automatic exact backend policy

`auto_backend_policy.pdf` and `auto_backend_policy.svg` show the deterministic
selection performed by the CP2K `MODE AUTO` policy.  The physical complete
regular k-point mesh remains the starting point in every branch; the selected
backend changes only exact storage, contraction ordering, transforms, and MPI
ownership.  Implicit qualification is deliberately disabled, and the dense
oracle remains selectable in manual mode.  The one-point branch retains the
compact dense ACP contraction.  For periodic multi-point meshes with ACP
parameters, automatic mode selects bounded forward Bloch batches and the
sparse projector--image reverse; the bounded ACP cache remains active.

Regenerate both formats from the repository root with:

```bash
python3 scripts/plot_auto_backend_policy.py
```

The script fixes the SVG hash salt and output metadata so that identical input
and Matplotlib versions produce byte-identical files.  Verify the retained
outputs with:

```bash
sha256sum -c validation/accelerated_exchange/figures/SHA256SUMS
```
