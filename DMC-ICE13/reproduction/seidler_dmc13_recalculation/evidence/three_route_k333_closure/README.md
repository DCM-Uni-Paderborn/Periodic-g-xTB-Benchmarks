# Complete 3 x 3 x 3 three-route closure

This gate combines the package-resident, hash-qualified raw data into one
direct all-phase comparison:

- the final author `pbc` `save_tblite` CLI;
- the current CP2K-integration `save_tblite` CLI;
- the CP2K-native symmetry-reduced Bloch implementation.

`generate_and_verify.py` reparses every current-CLI JSON and CP2K output,
checks the recorded executable and input hashes, complete 13-structure
coverage, normal termination, primitive-cell water counts, absolute-energy
parity, and all Ih-referenced DMC-ICE13 relative energies.  It has no external
archive dependency.  The generated tables make the
two logically distinct differences explicit: the tiny current-CLI/native
roundoff residual and the small author-`pbc`/current-provider model revision.

Run:

```text
python3 generate_and_verify.py
```

After regenerating the tables, refresh `SHA256SUMS` from this directory and
verify it with `shasum -a 256 -c SHA256SUMS`.
