make_fixture <- function(n = 250, p = 10, m = 50, seed = 42L) {
  set.seed(seed)
  X <- scale(matrix(rnorm(n * p), n, p)); colnames(X) <- paste0("sig", seq_len(p))
  Y <- scale(matrix(rnorm(n * m), n, m)); colnames(Y) <- paste0("s", seq_len(m))
  rownames(X) <- rownames(Y) <- paste0("g", seq_len(n))
  list(X = X, Y = Y)
}

test_that("pure-R batched matches single-call at ncores=1", {
  d <- make_fixture()
  r_single <- SecAct:::.ridge_pureR(d$X, d$Y, lambda = 1, nrand = 80L)
  r_batch  <- SecAct:::.ridge_pureR_batch(d$X, d$Y, lambda = 1, nrand = 80L,
                                          batch_size = 11L)
  for (nm in c("beta", "se", "zscore", "pvalue")) {
    expect_equal(r_batch[[nm]], r_single[[nm]], info = paste("field:", nm))
  }
})

test_that("dispatcher routes to cpu-pure when only pure-R available", {
  d <- make_fixture(m = 20)
  r <- SecAct:::.ridge_batch_dispatch(d$X, d$Y, lambda = 1, nrand = 40L,
                                      backend = "cpu-pure",
                                      batch_size = 7L)
  r_ref <- SecAct:::.ridge_pureR(d$X, d$Y, lambda = 1, nrand = 40L)
  expect_equal(r$beta, r_ref$beta)
})

test_that("dispatcher rejects srand on cpu-pure", {
  d <- make_fixture(m = 10)
  expect_error(
    SecAct:::.ridge_batch_dispatch(d$X, d$Y, lambda = 1, nrand = 10L,
                                   backend = "cpu-pure",
                                   rng_method = "srand",
                                   batch_size = 5L),
    "pure-R supports only"
  )
})

test_that("pure-R batch writes HDF5 round-trip matches in-memory", {
  skip_if_not_installed("rhdf5")
  d <- make_fixture(m = 16)
  out <- tempfile(fileext = ".h5"); on.exit(unlink(out))
  meta <- SecAct:::.ridge_pureR_batch(d$X, d$Y, lambda = 1, nrand = 30L,
                                      batch_size = 5L, output_h5 = out)
  expect_equal(meta$m, ncol(d$Y))
  r_mem <- SecAct:::.ridge_pureR(d$X, d$Y, lambda = 1, nrand = 30L)
  for (nm in c("beta", "se", "zscore", "pvalue")) {
    v <- rhdf5::h5read(out, nm)
    ref <- r_mem[[nm]]; dimnames(ref) <- NULL
    expect_equal(v, ref, info = paste("field:", nm))
  }
})
