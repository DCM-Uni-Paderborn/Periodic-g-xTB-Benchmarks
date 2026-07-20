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

At the time of this archived gate, the ice-Ih `9 x 9 x 9` endpoint predated
checkpoint activation and was therefore not restartable.  That calculation was
subsequently stopped and `9 x 9 x 9` removed from the active Part-I schedule.
The ice-VII `9 x 9 x 9` input remains here only as immutable restart-safety
evidence; it is not queued.  The active sub-cap ice-XIV and ice-XI inputs write
a checkpoint after every completed SCF cycle, so the validated procedure can
resume them after their first completed cycle without repeating it.

Run the archived gate from the repository root with

```bash
python3 validation/restart_operational_safety_20260720/verify_restart_operational_safety.py
```
