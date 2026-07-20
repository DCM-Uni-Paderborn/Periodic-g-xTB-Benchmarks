# CP2K k-point restart operational safety

This gate validates the restart procedure used for the remaining large
Part-I DMC-ICE13 endpoints.  It is deliberately separate from the scientific
energy comparison: a restart is accepted only when CP2K explicitly reports a
validated k-point restart, not merely because the calculation eventually
converges to the same energy.

The test exposed two operational details that are easy to miss:

- CP2K must be started in the result directory so that newly written `.kp`
  checkpoints stay with the corresponding input, output, hashes, and affinity
  proof.
- In this input path, quoting `WFN_RESTART_FILE_NAME` makes CP2K treat the
  otherwise existing file as missing.  The archived one-line ablation therefore
  falls back to a cold guess and needs twelve SCF steps.  The unquoted control
  accepts the same checkpoint and needs one step.

`prepare_same_mesh_restart.py` emits only short, unquoted restart names and
rejects whitespace or quotes.  `resume_pinned_cp2k.sh` preserves an immutable
source checkpoint, archives the preceding attempt, launches through the
singleton-affinity wrapper, and fails unless CP2K prints either the strict
same-mesh or validated BvK-transfer acceptance marker.  The latter remains an
exact restart route and is included because the formal selector test requests
that mode explicitly.

The already-running ice-Ih `9 x 9 x 9` endpoint predates checkpoint activation
and is therefore not restartable.  It is left untouched.  The queued ice-VII,
ice-XIV, ice-XI, and ice-XVII inputs write a checkpoint after every completed
SCF cycle; once the first such cycle has completed, this validated procedure
can resume them without redoing the completed cycles.

Run the archived gate from the repository root with

```bash
python3 validation/restart_operational_safety_20260720/verify_restart_operational_safety.py
```
