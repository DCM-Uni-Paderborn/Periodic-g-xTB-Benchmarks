# FFT, symmetry-phase-cache, and cross-mesh restart evidence bundle

This index groups three independently verifiable acceleration archives without
combining their claims:

| Archive | Qualification | Explicit boundary |
|---|---|---|
| `../mixed_radix_fft_20260717/` | PASS: separable and mixed-radix FFT regular-mesh exchange transforms agree with the dense oracle for energy, force, and stress in the archived RKS/UKS and 1D/2D/3D matrix | Correctness archive; no speedup claim from the short serial runs |
| `../symmetry_phase_cache_20260717/` | PASS: cached and uncached symmetry-star phase paths are numerically equivalent in the archived matrix | Correctness archive; no speedup claim from the short serial runs |
| `../crossmesh_restart/` | PASS: opt-in metadata-validated BvK density transfer, including official regression, malformed-input fallback, RKS/UKS, force/stress and MPI-2 | Initial guess only; not guaranteed N-representable, and not enabled by default |

Each target directory has its own README, raw data, source patch, verification
helpers, provenance, and SHA-256 manifest.  The archives deliberately retain
the dense transform and cold restart as numerical oracles.  No manuscript,
SI, or application-benchmark content is part of this bundle.
