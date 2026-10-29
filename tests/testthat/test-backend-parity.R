test_that("cpu-pure and cpu-fast backends produce identical z-scores in mt19937 mode", {
  skip_if_not_installed("RidgeFast")

  set.seed(7)
  n <- 60; p <- 10; m <- 4
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("sig", 1:p)
  Y <- matrix(rnorm(n * m), n, m); colnames(Y) <- paste0("s", 1:m)
  rownames(X) <- rownames(Y) <- paste0("g", 1:n)
  X <- scale(X); Y <- scale(Y)

  r_pure <- SecAct:::.ridge_dispatch(X, Y, lambda = 100, nrand = 200,
                                     ncores = 1L, rng_method = "mt19937",
                                     backend = "cpu-pure")
  r_fast <- SecAct:::.ridge_dispatch(X, Y, lambda = 100, nrand = 200,
                                     ncores = 1L, rng_method = "mt19937",
                                     backend = "cpu-fast")

  # bit-identical up to machine epsilon (double precision)
  expect_lt(max(abs(r_pure$beta   - r_fast$beta)),   1e-12)
  expect_lt(max(abs(r_pure$se     - r_fast$se)),     1e-12)
  expect_lt(max(abs(r_pure$zscore - r_fast$zscore)), 1e-10)
  expect_identical(r_pure$pvalue, r_fast$pvalue)  # integer counts; exact match
})

test_that("pure-R backend is deterministic across runs", {
  set.seed(1)
  n <- 40; p <- 6; m <- 3
  X <- scale(matrix(rnorm(n * p), n, p))
  Y <- scale(matrix(rnorm(n * m), n, m))
  colnames(X) <- paste0("sig", 1:p); colnames(Y) <- paste0("s", 1:m)

  r1 <- SecAct:::.ridge_pureR(X, Y, lambda = 50, nrand = 100)
  r2 <- SecAct:::.ridge_pureR(X, Y, lambda = 50, nrand = 100)

  expect_identical(r1, r2)
})

test_that("dispatcher rejects unsupported rng_method on cpu-pure", {
  X <- scale(matrix(rnorm(30), 10, 3))
  Y <- scale(matrix(rnorm(20), 10, 2))
  colnames(X) <- paste0("s", 1:3); colnames(Y) <- paste0("y", 1:2)
  expect_error(
    SecAct:::.ridge_dispatch(X, Y, lambda = 1, nrand = 10,
                             ncores = 1L, rng_method = "srand",
                             backend = "cpu-pure"),
    "pure-R supports only 'mt19937'"
  )
})

# SecAct 1.0.0 SecAct.inference.r permutation loop (source of the published z-scores)
legacy_ref <- function(X, Y, lambda, nrand) {
  n <- nrow(Y)
  R <- chol(crossprod(X) + lambda * diag(ncol(X)))
  beta <- backsolve(R, forwardsolve(t(R), crossprod(X, Y)))
  aver <- aver_sq <- cnt <- 0
  for (i in 1:nrand) {
    set.seed(i)
    br <- backsolve(R, forwardsolve(t(R), crossprod(X, Y[sample.int(n), , drop = FALSE])))
    aver <- aver + br; aver_sq <- aver_sq + br^2; cnt <- cnt + (abs(br) >= abs(beta))
  }
  aver <- aver / nrand
  se <- sqrt(aver_sq / nrand - aver^2)
  list(beta = beta, se = se, zscore = (beta - aver) / se, pvalue = (cnt + 1) / (nrand + 1))
}

legacy_data <- function() {
  set.seed(11)
  n <- 50; p <- 8; m <- 3
  X <- scale(matrix(rnorm(n * p), n, p)); Y <- scale(matrix(rnorm(n * m), n, m))
  colnames(X) <- paste0("sig", 1:p); colnames(Y) <- paste0("s", 1:m)
  list(X = X, Y = Y)
}

test_that("legacy_r reproduces the SecAct 1.0.0 pure-R permutation stream", {
  d <- legacy_data()
  ref <- legacy_ref(d$X, d$Y, lambda = 100, nrand = 200)
  res <- SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 200,
                                  rng_method = "legacy_r", backend = "cpu-pure")
  expect_lt(max(abs(res$beta   - ref$beta)),   1e-12)
  expect_lt(max(abs(res$se     - ref$se)),     1e-12)
  expect_lt(max(abs(res$zscore - ref$zscore)), 1e-10)
  expect_equal(unname(res$pvalue), unname(ref$pvalue))
})

test_that("legacy_r leaves the caller's RNG state and kind untouched", {
  d <- legacy_data()
  RNGkind("L'Ecuyer-CMRG"); set.seed(99)
  on.exit(RNGkind("default", "default", "default"), add = TRUE)
  kind <- RNGkind(); seed <- .Random.seed
  res <- SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                  rng_method = "legacy_r", backend = "cpu-pure")
  expect_identical(RNGkind(), kind)
  expect_identical(.Random.seed, seed)
  # stream is independent of the caller's RNG kind
  RNGkind("default", "default", "default")
  res2 <- SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                   rng_method = "legacy_r", backend = "cpu-pure")
  expect_identical(res, res2)
})

test_that("legacy_r routes auto to cpu-pure and rejects explicit accelerators", {
  d <- legacy_data()
  r_auto <- suppressMessages(SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                                      rng_method = "legacy_r", backend = "auto"))
  r_pure <- SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                     rng_method = "legacy_r", backend = "cpu-pure")
  expect_identical(r_auto, r_pure)
  expect_error(SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                        rng_method = "legacy_r", backend = "cpu-fast"),
               "only available with backend='cpu-pure'")
})

test_that("legacy_r batched path matches the single-call path", {
  d <- legacy_data()
  r_one <- SecAct:::.ridge_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                    rng_method = "legacy_r", backend = "cpu-pure")
  r_bat <- SecAct:::.ridge_batch_dispatch(d$X, d$Y, lambda = 100, nrand = 50,
                                          rng_method = "legacy_r", backend = "cpu-pure",
                                          batch_size = 2L)
  for (nm in c("beta", "se", "zscore", "pvalue")) expect_equal(r_bat[[nm]], r_one[[nm]], tolerance = 1e-12)
})
