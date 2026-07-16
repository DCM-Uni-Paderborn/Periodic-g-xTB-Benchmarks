# Focused provider qualification run

`focused_exchange.stdout` and `focused_exchange.stderr` are the unedited
streams from the command in `command.txt`.  GNU Fortran emitted known
array-temporary diagnostics on stderr; the focused test itself ended with
`bvk_exchange_supercell [PASSED]` and return code zero.

This run is the qualification basis for the provider cache/planner and the
matrix-lean forward stream only.  Its storage assertion proves only that the
stream does not retain full-k-space density/overlap input arrays; it does not
measure the three retained BvK-image tensors, the two dense phase tables, or
total process memory.  True bounded-memory R/image batching, the
reduced-memory reverse stream and CP2K consumer integration are outside its
passed scope.
