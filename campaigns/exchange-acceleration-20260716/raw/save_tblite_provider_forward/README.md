# Focused provider qualification run

`focused_exchange.stdout` and `focused_exchange.stderr` are the unedited
streams from the command in `command.txt`.  GNU Fortran emitted known
array-temporary diagnostics on stderr; the focused test itself ended with
`bvk_exchange_supercell [PASSED]` and return code zero.

This run is the qualification basis for the provider cache/planner and the
bounded-memory forward stream only.  The reduced-memory reverse stream and
the CP2K consumer integration are outside its passed scope.
