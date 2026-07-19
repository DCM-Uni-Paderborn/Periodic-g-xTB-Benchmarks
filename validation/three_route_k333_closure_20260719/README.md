# Complete 3 x 3 x 3 three-route closure

This gate combines two independently verified archives into one direct
all-phase comparison:

- the final author `pbc` `save_tblite` CLI;
- the current CP2K-integration `save_tblite` CLI;
- the CP2K-native symmetry-reduced Bloch implementation.

`generate_and_verify.py` first reruns both prerequisite provenance verifiers.
It then reparses every CP2K output, checks complete 13-structure coverage,
normal termination, primitive-cell water counts, absolute-energy parity, and
all Ih-referenced DMC-ICE13 relative energies.  The generated tables make the
two logically distinct differences explicit: the tiny current-CLI/native
roundoff residual and the small author-`pbc`/current-provider model revision.

Run:

```text
python3 generate_and_verify.py
```
