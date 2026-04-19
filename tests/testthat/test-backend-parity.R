test_that("cpu-pure and cpu-fast backends produce identical z-scores in mt19937 mode", {
  skip_if_not_installed("RidgeRegFast")

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
