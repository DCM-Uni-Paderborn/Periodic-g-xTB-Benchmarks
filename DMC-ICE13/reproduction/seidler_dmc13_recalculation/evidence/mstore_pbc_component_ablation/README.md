# `mstore-inorganic` versus `pbc` component ablation

This directory isolates the origin of the unusually large difference between
the historical `mstore-inorganic` and later `pbc` DMC-ICE13 energies.  The
same explicit `2 x 2 x 2` BvK cells for ice Ih and ice VII were evaluated with
the exact recorded source states and CLI accuracy `0.1`.  Four independently
self-consistent parameterizations were used for each provider:

- the complete g-xTB model;
- exchange disabled;
- ACP disabled;
- exchange and ACP disabled.

The full same-mesh ice-VII relative-energy difference is
`-148.119354589090 kJ mol-1 H2O-1` (`pbc` minus `mstore-inorganic`).  Disabling
exchange reduces its magnitude by `98.5692%`, to
`2.119286496877 kJ mol-1 H2O-1`.  Disabling ACP alone reduces the full gap by
only `4.2775%`.  With both exchange and ACP disabled, the gap is
`1.559210900355 kJ mol-1 H2O-1`, a `98.9473%` reduction.

This establishes that the different sparse-mesh behavior of the two source
states is controlled almost entirely by their different exchange paths.  It
does not indicate an error in the CP2K interface: the independently archived
same-provider direct-CLI/CP2K matrix agrees through `4 x 4 x 4`.  Because every
ablation was reconverged self-consistently, the numbers are a coupled response
test rather than an additive fixed-density energy decomposition.

`evaluate_component_matrix.py` rechecks all sixteen calculations, executable
and input hashes, printed convergence thresholds, SCC termination, archived
full-model values, and the derived relative-energy gaps.  It regenerates
`component_relative_energies.csv` and `verification.json`.
