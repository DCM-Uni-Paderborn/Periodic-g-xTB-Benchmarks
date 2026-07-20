# Part-I implementation audit

This directory provides a single machine-readable gate over the completed
periodic g-xTB reference-implementation tests. It aggregates, without changing
their individual tolerances, the following independent checks:

- direct save_tblite CLI versus CP2K-native absolute-energy and independently
  reconstructed relative-energy parity for all 52 points from `1 x 1 x 1`
  through `4 x 4 x 4`;
- in-process direct-CLI parity for energy, Cartesian gradients, and the tested
  periodic virial using the qualified CP2K and CLI binaries;
- strict same-mesh and validated regular-mesh restart equivalence against
  independent cold-start energies and SCF iteration counts;
- operational restart recovery with singleton affinity, immutable checkpoint
  preservation, explicit CP2K acceptance markers, and the quoted-file-name
  negative control;
- reciprocal source-patch attribution of the historical Wigner--Seitz
  self-image error and the remaining second-order minimum-image model change,
  including an independent ice-XVII cross-check;
- Git-aware completeness and hash verification of the independent Seidler
  recalculation package, including all raw text outputs;
- source-level and exact-arithmetic equivalence of every archived MacDonald
  grid to the corresponding Gamma-supercell BvK folding grid;
- build-level identity of the provider source, static provider archive, and
  selected CP2K provider revision;
- exact successor-commit diffs proving that the current CP2K and
  `save_tblite` branch heads do not alter the qualified MacDonald/MOPAC
  production-energy path;
- independent decimal reconstruction of BvK normalization, per-water
  referencing, unit conversion, DMC errors, and aggregate statistics;
- same-accuracy and tighter-SCC comparisons;
- native Bloch k points versus explicit CP2K Gamma-point BvK supercells for
  both ice XVII and the high-residual ice-VII case;
- full-grid, K290, SPGLIB, and time-reversal equivalence;
- analytic energy, force, virial, and stress paths, including 1D and 2D PBC;
- the complete ice-XVII native-$2\times2\times2$ full-grid/SPGLIB derivative
  identity, independent collective-displacement and homogeneous-strain finite
  differences, and direct 144-atom `save_tblite` CLI derivative comparison;
- frozen-production-binary repetition of the complete ice-XVII full-grid,
  symmetry-reduced, force-displacement, and homogeneous-strain sequence, with
  exact executable and input hashes and independent force/stress finite
  differences;
- exchange/ACP component ablations and model-provider attribution;
- reciprocal Wigner--Seitz one-patch builds that reproduce the author `pbc`
  energies and explain more than 95% of the historical sparse-mesh branch gap;
- source-level periodic Hamiltonian, exchange, ACP, Coulomb, dispersion,
  repulsion, and Wigner--Seitz tests;
- the complete 75-case Hamiltonian group after resolving the inherited
  nonperiodic CeCl3 finite-difference threshold false negative;
- exact geometry equivalence and internal consistency of the current adaptive
  DMC-ICE13 statistics;
- deterministic reconstruction of the qualified phase-wise progress values at
  mesh limits 6x6x6, 7x7x7, and 8x8x8 and identity of the latest snapshot with
  the current adaptive phase table;
- independent reconstruction of the numerical hierarchy separating the
  CLI/native residual, pbc source-state evolution, `mstore-inorganic` model
  evolution, and the still-provisional dense-mesh DMC statistic.

Run from the repository root with

```bash
python3 validation/implementation_audit_20260720/verify_implementation_audit.py
```

For a stronger clean-checkout qualification that actually reruns every
completed child verifier, regenerates the derived low-k, three-route, and
Gamma-supercell results, checks the selected SHA-256 manifests, and requires
the tracked worktree to remain unchanged, run

```bash
python3 validation/implementation_audit_20260720/requalify_completed_evidence.py
```

The resulting `requalification.json` records every invoked check and its
deterministic output hashes.  The driver refuses to start from a dirty tracked
checkout and fails if any regenerated evidence differs from the published
files.

`verification.json` reports the completed qualification gates separately from
the still-running science endpoints. The direct CLI/native diagnostic matrix
is complete through `4 x 4 x 4`, and no declared implementation diagnostic
remains pending. A passing audit therefore
does not turn the provisional adaptive DMC-ICE13 statistic into a final one.

The pending science endpoints are derived from the qualified adaptive table,
not from a predeclared mesh schedule. The audit fails closed unless the table
contains exactly the twelve DMC-ICE13 polymorphs, every convergence flag agrees
with the one-step threshold of `0.10 kJ mol-1` per water, and the declared final
state agrees with the number of converged phases. A missing same-build ice-Ih
reference is inserted ahead of the corresponding phase endpoint automatically.
The present Part-I schedule is capped at `8 x 8 x 8`; an unresolved phase at
that cap is reported separately instead of generating a `9 x 9 x 9` job.
