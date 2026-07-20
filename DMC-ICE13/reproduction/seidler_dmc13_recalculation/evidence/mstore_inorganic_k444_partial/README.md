# Historical mstore-inorganic 4x4x4 partial matrix

This directory preserves the direct save_tblite CLI calculations made with
the recorded historical `mstore-inorganic` executable.  Ih and eleven of the
twelve DMC-ICE13 benchmark phases completed.  Phase XIII was terminated by the
operating system with exit status -9 before the first SCC result because the
explicit 4x4x4 BvK supercell exceeded the available local memory.

Run:

```text
python3 evaluate_partial_k444.py
```

A passing verification proves the hashes, command line, convergence
thresholds, successful SCC terminations, and the exact failed XIII record.  It
also compares the same eleven completed phases with the qualified current-pbc
CP2K-native values at 2x2x2, 3x3x3, and 4x4x4.  It deliberately sets
`full_matrix_complete` and `usable_for_full_benchmark_statistics` to false.

The same-eleven-phase diagnostics show a strongly shrinking historical-source
gap as the mesh is refined.  The historical 3x3x3-to-4x4x4 changes are still
large, so the smaller sparse-mesh DMC error is finite-size error cancellation;
it is not a complete converged benchmark and does not contradict the direct
current-pbc CLI/CP2K-native parity proof.
