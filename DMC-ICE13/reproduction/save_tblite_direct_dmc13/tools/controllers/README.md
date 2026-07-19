# Qualified adaptive DMC-ICE13 completion

These controllers implement the production ordering and qualification policy
for the response-corrected Part-I DMC-ICE13 series.

`run_strict_adaptive_completion.sh` first qualifies ice VII at `8 x 8 x 8`,
then ice Ih at the same mesh, and only then calls the explicit CP2K
Gamma-supercell oracle.  It accepts a result only after normal termination and
an exact match to the required CP2K executable SHA-256.  The remaining series
uses the first adjacent pair satisfying
`|Delta R| <= 0.10 kJ mol^-1 per H2O`, retains the denser endpoint, and advances
only unresolved phases.  Mesh 8 and above are serialized to limit memory use.

Concurrency below that boundary is bounded by archived peak-memory data: all
twelve non-reference phases run together at `4^3` and `5^3`, at most six run
together at `6^3`, and at most two at `7^3`.  Every science process still
receives a disjoint singleton CPU through the reservation launcher.

`restart_parallel_completion_after_current.sh` safely replaces an older
suspended controller only after its current hash-qualified calculation has
completed and no live CP2K or MPI process remains.

`run_gamma_supercell_oracle.sh` refuses to run until the qualified VII/Ih
`8 x 8 x 8` pair exists.  It then evaluates the explicit Gamma-only
Born--von--Karman supercell with its independently frozen input hash.

`run_independent_adaptive_verification.sh` reparses the selected endpoints with
an independent verifier after the production controller passes.

All calculations are delegated to the archived pinned launcher.  That
launcher is responsible for disjoint singleton CPU reservations,
`OMP/BLAS=1`, the pre-execution `/proc` affinity proof, and cross-process
reservation locking.

`launch_pinned_command.sh` applies the same fail-closed policy to direct
save_tblite commands.  It records the executable and final input SHA-256
before acquiring the reservation and refuses overlapping live CP2K, MPI, or
save_tblite processes.  This keeps the direct CLI/native energy comparison
reproducible without treating it as manuscript benchmark data.

`monitor_qualified_mixed_mae.sh` can additionally receive an explicitly
labelled earlier-run comparator through `EARLIER_RUN_MAE`.  The comparator is
reported separately from the same-mesh paper value; it never participates in
endpoint selection or qualification.

Run the synthetic fail-closed controller test with:

```bash
python3 test_adaptive_controller.py -v
```

The test verifies the VII--Ih--oracle priority order, the `0.10` threshold,
phase-local pruning, exact binary qualification, and successful final status
without launching CP2K.
