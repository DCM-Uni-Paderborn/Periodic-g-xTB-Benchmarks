# Historical mstore-inorganic accuracy equivalence

This gate tests whether the direct-CLI DMC-ICE13 difference between the
historical `mstore-inorganic` source state and the current `pbc`-derived source
could be caused by the CLI accuracy setting.

All thirteen structures, including Ih, were evaluated at `3 x 3 x 3` with the
same independently rebuilt historical executable and identical POSCARs.  The
two matrices differ only in the requested CLI accuracy, `0.1` versus `0.01`.
Every admitted row carries a successful exit state plus executable and input
hashes.

Run

```text
python3 verify_mstore_accuracy_equivalence.py
```

to regenerate `verification.json`.  The maximum energy response to the tighter
threshold is `6.5e-11` hartree per evaluated supercell, corresponding after
same-mesh Ih referencing to only `2.26e-10 kJ mol-1 H2O-1`.  The direct
`result.json` energy is already in hartree and is not converted a second time.
This rules out the SCC stopping threshold as the origin of the substantial
`mstore-inorganic`/`pbc` DMC difference.
