## SecAct (development version)

* `SecAct.activity.inference()` and its ST/scRNA-seq variants gain
  `rng_method = "legacy_r"`, which reproduces the SecAct 1.0.0 pure-R
  permutation stream (`set.seed(i); sample.int(n)` per permutation) and
  therefore the published SecAct z-scores exactly. Available on the
  pure-R backend only; `backend = "auto"` falls back to `"cpu-pure"`.

## SecAct 1.1.0

*2026-07-28*

* Add a new tutorial for analyzing Visium HD data.

## SecAct 1.0.0

*2024-08-03*

* Release SecAct package.
