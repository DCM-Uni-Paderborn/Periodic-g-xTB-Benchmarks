# g-xTB PBC V2: symmetry, memory, and scaling roadmap

Status: deferred design and manuscript note. Implementation starts only after
the V1 benchmark campaign has been completed and frozen. This file belongs to
the new g-xTB work and must not be copied into or used to modify the existing
GFN2-xTB Overleaf project.

## Scientific baseline

V1 is the correctness oracle. CP2K diagonalizes only the irreducible k points,
expands the density and overlap matrices to the complete coupled mesh, calls
save_tblite on that mesh, and folds the response back to the irreducible
representation. V2 may change storage, evaluation order, communication, and
transform algorithms, but not the physical sum or the converged result.

The distinction from ordinary semilocal DFT should be made explicit in the
paper. LDA, GGA, and meta-GGA first form symmetry-completed real-space fields
such as

\[
  n(\mathbf r)=\sum_{\mathbf k\in\mathrm{IBZ}}
  w_{\mathbf k} n_{\mathbf k}(\mathbf r),
\]

after which the Hartree and semilocal exchange-correlation potentials no
longer retain an individual k-point label. A g-xTB Fock-like contribution has
the schematic coupled form

\[
  K(\mathbf k)=\sum_{\mathbf k'}
  L(\mathbf k,\mathbf k')P(\mathbf k').
\]

A scalar star weight therefore cannot replace the rotated density matrices:
the atom and basis permutation, the basis-space transformation \(U_g\), Bloch
phases, and time-reversal conjugation must all be retained. Hartree terms can
be spatially nonlocal without having this k-pair structure. Hybrid/exact
exchange in DFT is the closer algorithmic analogue to the g-xTB term.

V2 must preserve every physical star contribution. "No persistent full mesh"
means that symmetry-related matrices are generated and contracted transiently;
it does not mean that their contributions are omitted or replaced by scalar
weights.

## V2 work packages

Each work package must have an independent runtime/build switch until its
individual and combined validation gates have passed.

### 1. On-the-fly star contraction

- Start from irreducible \(P_i\) and \(S_i\) representatives and their recorded
  star maps.
- Generate a required star member transiently using \(U_g\), atom/basis
  permutations, Bloch phases, and time reversal.
- Contract it immediately into the appropriate exchange/Fock accumulator and
  release the transient matrix instead of materializing complete \(P\), \(S\),
  and response meshes.
- Fold the accumulated response back to the irreducible representation with
  the same normalization and star multiplicities as V1.
- Only after this streaming path matches V1 should pair-orbit or
  \(\mathbf q=\mathbf k-\mathbf k'\) symmetry be used to avoid algebraically
  equivalent repeated contractions.

This first restores memory savings and some transformation savings. The
physical k-pair sum remains; further reduction requires a proven symmetry of
the two-index kernel rather than an assumed star weight.

### 2. Immutable-data caches

Cache only data whose lifetime and invalidation key are explicit:

- full/irreducible k maps, star membership, symmetry operations, time-reversal
  partners, atom and basis permutations;
- invariant Bloch phases and \(U_g\) transformations;
- k/R index maps, transform normalization, and FFT plans;
- ACP image/support tables and other geometry-dependent neighbour metadata.

The cache fingerprint must include at least cell, coordinates and atom order,
k mesh and shift, symmetry operations and tolerance, basis layout/fingerprint,
model parameters, ACP support, and every q-vSZP state that can change the
basis. Geometry, cell, basis, charge-dependent basis, mesh, or symmetry changes
must invalidate the affected entries before reuse. Density/Fock matrices,
charges, residuals, and mixer history are dynamic SCC state and must not be
hidden in an immutable-geometry cache.

Tests must cover both reuse and deliberate invalidation during ionic and cell
steps. A cache hit is never allowed to change a result relative to a forced
cache miss.

### 3. FFT k-to-R and R-to-k transforms

Replace dense transforms only for complete regular meshes whose indexing,
shift, normalization, and phase convention have been proven compatible.
Shifted even meshes require their explicit phase factor; they are not treated
as unshifted FFT grids. Keep the dense transform as a selectable oracle and as
a fallback for unsupported meshes.

Unit gates include forward/backward round trips, Parseval-like norm checks,
Hermiticity, time reversal, odd unshifted and even shifted meshes, and direct
matrix comparison with the dense transform. Integration gates cover converged
energies, Fock responses, forces, and stress.

### 4. MPI k groups

Distribute k points, star members, or proven k-pair/q orbits without requiring
every rank to retain the complete mesh. Define ownership and reductions for
energy, Fock, forces, stress, and SCC norms explicitly. A one-group execution
must reproduce the serial V1 oracle; multi-group reductions must remain inside
a frozen round-off envelope and should offer a reproducible reduction mode for
debugging.

Test rank counts that do and do not divide the number of full and irreducible k
points, empty local partitions, different rank mappings, and repeated runs.
Measure communication and imbalance separately from kernel time.

### 5. Validated 5x5x5 to 6x6x6 density restart

The two k meshes are not nested. Density matrices must never be copied by k
index or relabelled. A valid restart has to pass through a well-defined common
representation, for example:

1. transform the converged \(5^3\) density to a localized real-space/BvK
   representation;
2. apply a documented truncation/interpolation rule;
3. evaluate the resulting initial density on the shifted \(6^3\) mesh;
4. restore symmetry, Hermiticity, electron number, and the required trace
   normalization before the first SCC step.

A lower-risk first stage may restart only mesh-invariant atomic charges and
multipoles while rebuilding the density matrix. DIIS/Broyden history must not
be transferred blindly when vector dimensions or mesh conventions change.

Warm and cold starts must converge to the same physical state and final
energy, charges, density, forces, and stress. Iteration or wall-time savings
are reported only after that equality gate passes. Production benchmark data
must remain cold-start data until the restart path is independently accepted.

## Validation contract

The serial V1 expanded-full-mesh path and its frozen executable/results remain
available as the oracle. V2 acceptance is staged:

1. **Algebraic kernel tests:** fixed complex matrices; expansion/foldback,
   streaming contraction, dense/FFT transforms, cache hit/miss, and serial/MPI
   reductions compared before SCC.
2. **Single-step integration:** identical density input; compare total and
   decomposed energies, Fock matrices, populations/multipoles, electron count,
   Hermiticity, forces, and stress.
3. **Converged SCC:** compare final observables, residuals and state identity
   for the save_tblite Fock-DIIS default and the CP2K density/Fock mixer
   alternatives. Iteration trajectories may differ only within a documented
   round-off explanation; endpoints must agree.
4. **Derivative gates:** analytic forces and stress against both V1 and
   independent central finite differences, including cell changes that force
   cache invalidation.
5. **Restart gates:** cold versus warm \(5^3\!\to6^3\), including failed or
   deliberately poor restart input and clean fallback.
6. **Regression gates:** GFN1, GFN2, implicit/explicit Gamma, and unaffected DFT
   paths must remain unchanged.

The matrix must contain at least:

- P1 and high-symmetry structures;
- inversion and time-reversal cases;
- a nonsymmorphic operation with fractional translation;
- complex Bloch phases and atom/basis permutations;
- odd unshifted meshes and shifted even meshes;
- an anisotropic cell that exercises the full ACP image/support range;
- DMC-ICE13, X23b, and LC12 representatives, followed by full benchmark-level
  spot checks.

Before implementation results are inspected, freeze separate tolerances for
algebraic double-precision tests and converged SCF tests from repeated V1
baseline noise. Do not relax them after seeing a V2 discrepancy. As initial
targets, kernel-level comparisons should be near double-precision round-off;
end-to-end energy, force, and stress envelopes should be no looser than the
existing full-grid-versus-SPGLIB and finite-difference validation permits.

All feature combinations are required, not only each feature in isolation:
streaming plus cache, streaming plus FFT, streaming plus MPI, all three
together, and each combination with and without restart.

## Performance and memory evidence

Correctness gates precede performance claims. Use identical commits, compiler,
libraries, hardware allocation, convergence thresholds, and input hashes.
Record at least:

- \(N_k^{\rm full}\), \(N_k^{\rm irr}\), star and pair/q-orbit counts;
- peak resident memory and allocated matrix bytes per rank and in total;
- time split into symmetry transforms, k/R transforms, ACP preparation,
  exchange contraction, communication, diagonalization, and other SCC work;
- cold-cache and warm-cache timings;
- strong scaling across MPI k groups and load imbalance;
- cold-start versus restarted SCC iterations and wall time.

Benchmark small through dense meshes so setup overhead and asymptotic behavior
are both visible. Report cases where V2 is slower, and retain an automatic V1
or dense-transform fallback when that is the safer/faster regime.

## Paper record

The new g-xTB paper should preserve the following claims and distinctions:

- V1 is the explicit full-mesh correctness implementation behind CP2K's
  SPGLIB-reduced interface.
- V2 changes the algorithmic realization, not the g-xTB energy expression or
  the physical star sum.
- Semilocal DFT can reduce earlier because it forms symmetry-completed scalar
  or local real-space fields; exact/hybrid exchange is the appropriate DFT
  analogy for the coupled g-xTB exchange.
- Separate savings from irreducible diagonalization, avoided full-mesh
  materialization, pair/q symmetry, FFT transforms, caching, MPI distribution,
  and restart.
- Present numerical equivalence before speedup and scaling plots.
- Keep k-point convergence error distinct from algorithmic V1/V2 equivalence.
- Archive feature switches, build and source identities, input/output hashes,
  cache/restart metadata, MPI layout, and raw timing/memory counters with every
  reported result.

When V2 becomes executable, create a dated V2 campaign and a new immutable
build manifest. Never modify the hash-bound V1 manifest retroactively.
