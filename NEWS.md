## SecAct (development version)

* `SecAct.activity.inference()` and its ST/scRNA-seq variants gain
  `rng_method = "legacy_r"`, which reproduces the permutation stream of
  the SecAct 1.0.0 pure-R function `SecAct.inference.r`
  (`set.seed(i); sample.int(n)` per permutation). Results match that
  function to floating-point precision (~1e-15), not bit-identically.
  It does not reproduce the 1.0.0 default entry points, which permuted
  with C `srand(0)/rand()`. Available on the pure-R backend only;
  `backend = "auto"` falls back to `"cpu-pure"`.

## SecAct 1.1.0

*2026-07-28*

* Add a new tutorial for analyzing Visium HD data.

## SecAct 1.0.0

*2024-08-03*

* Release SecAct package.
