# `pbc` versus current-provider component attribution

This archive compares the final upstream `pbc` command-line executable with the
current `cp2k-integration` command-line executable for the identical DMC-ICE13
phase-VII `2 x 2 x 2` explicit Born--von Karman supercell.  All eight calculations
completed their SCC cycles on distinct, pre-verified singleton CPUs.

The four rows are separate self-consistent calculations: the complete model,
exchange disabled, ACP disabled, and both exchange and ACP disabled.  Therefore
the rows are diagnostic ablations, not an additive fixed-density decomposition.

The complete-model difference (`pbc - current`) is
`-3.1710954436e-3` hartree per supercell.  Disabling exchange reduces its absolute
magnitude to `5.3545428273e-5` hartree, a 98.31% reduction.  This identifies the
exchange-containing self-consistent path as the dominant source of the small
provider difference.  Disabling ACP instead increases the difference, so ACP
partly compensates it for this structure.  The non-additive double ablation must
not be interpreted by subtracting independent component energies.

Run `python3 verify_attribution.py` to recheck all energies, digests, convergence
markers, executable identities, and disjoint CPU affinities and to regenerate
`verification.json`.  `SHA256SUMS` preserves the raw archive produced on Terok.
