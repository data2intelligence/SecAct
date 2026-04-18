# =========================================================
# Pure-R ridge regression with permutation testing
# and backend dispatch to optional accelerators.
# =========================================================

#' Pure-R ridge + permutation kernel
#'
#' Solves \eqn{\beta = (X'X + \lambda I)^{-1} X'Y} via Cholesky, then
#' estimates null via Y-row permutations.
#'
#' @param X Numeric matrix, n x p, column-scaled.
#' @param Y Numeric matrix, n x m, column-scaled.
#' @param lambda Ridge penalty.
#' @param nrand Number of permutations.
#' @param mode \code{"canonical"} uses GSL-compatible MT19937 (seed 0)
#'   so permutations match RidgeRegFast / RidgeRegCuda bitwise.
#'   \code{"fast"} uses R's builtin RNG per-permutation (faster setup,
#'   not cross-backend reproducible).
#' @return list(beta, se, zscore, pvalue) — each a p x m matrix.
#' @keywords internal
.ridge_pureR <- function(X, Y, lambda, nrand, mode = c("canonical", "fast")) {
  mode <- match.arg(mode)
  n <- nrow(Y); p <- ncol(X); m <- ncol(Y)

  A <- crossprod(X) + lambda * diag(p)
  R <- chol(A)
  Rt <- t(R)
  beta <- backsolve(R, forwardsolve(Rt, crossprod(X, Y)))

  if (mode == "canonical") {
    perm_table <- .gsl_mt19937_perm_table(n, nrand)
    get_perm <- function(i) perm_table[i, ] + 1L
  } else {
    get_perm <- function(i) { set.seed(i); sample.int(n) }
  }

  aver <- matrix(0, p, m)
  aver_sq <- matrix(0, p, m)
  pvalue_count <- matrix(0, p, m)

  for (i in seq_len(nrand)) {
    perm <- get_perm(i)
    beta_rand <- backsolve(R, forwardsolve(Rt, crossprod(X, Y[perm, , drop = FALSE])))
    aver <- aver + beta_rand
    aver_sq <- aver_sq + beta_rand^2
    pvalue_count <- pvalue_count + (abs(beta_rand) >= abs(beta))
  }

  aver <- aver / nrand
  aver_sq <- aver_sq / nrand
  se <- sqrt(aver_sq - aver^2)
  zscore <- (beta - aver) / se
  pvalue <- (pvalue_count + 1) / (nrand + 1)

  dn <- list(colnames(X), colnames(Y))
  dimnames(beta)   <- dn
  dimnames(se)     <- dn
  dimnames(zscore) <- dn
  dimnames(pvalue) <- dn

  list(beta = beta, se = se, zscore = zscore, pvalue = pvalue)
}

#' Resolve backend based on installed packages and GPU availability
#' @keywords internal
.resolve_backend <- function(backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  backend <- match.arg(backend)
  if (backend != "auto") return(backend)

  if (requireNamespace("RidgeRegCuda", quietly = TRUE)) {
    gpu_ok <- tryCatch(
      RidgeRegCuda::check_cuda_available(),
      error = function(e) FALSE
    )
    if (isTRUE(gpu_ok)) return("gpu")
  }
  if (requireNamespace("RidgeRegFast", quietly = TRUE)) return("cpu-fast")
  "cpu-pure"
}

#' Dispatch ridge+permutation call to installed backend
#'
#' Picks GPU (RidgeRegCuda) > CPU-fast (RidgeRegFast) > CPU-pure (this package)
#' based on \code{backend}. All backends return the same list shape.
#'
#' @param X,Y Scaled numeric matrices (n x p, n x m).
#' @param lambda,nrand Ridge + permutation parameters.
#' @param ncores Number of threads (ignored for cpu-pure, used by accelerators).
#' @param rng_method Only \code{"gsl_r"} (MT19937 seed 0) is supported for
#'   canonical mode; accelerators may accept others in \code{mode="fast"}.
#' @param mode \code{"canonical"} (default): bit-identical across backends.
#'   \code{"fast"}: native RNG per backend, faster, not reproducible.
#' @param backend \code{"auto"} (default) or force a specific backend.
#' @return list(beta, se, zscore, pvalue)
#' @keywords internal
.ridge_dispatch <- function(X, Y, lambda, nrand,
                            ncores = 1L,
                            rng_method = "gsl_r",
                            mode = c("canonical", "fast"),
                            backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  mode <- match.arg(mode)
  chosen <- .resolve_backend(backend)
  message("[SecAct] ridge backend: ", chosen, " (mode=", mode, ")")

  if (chosen == "gpu") {
    RidgeRegCuda::Ridge.Reg.Fast(X = X, Y = Y, lambda = lambda, nrand = nrand,
                                 ncores = ncores, rng_method = rng_method, mode = mode)
  } else if (chosen == "cpu-fast") {
    RidgeRegFast::Ridge.Reg.Fast(X = X, Y = Y, lambda = lambda, nrand = nrand,
                                 ncores = ncores, rng_method = rng_method, mode = mode)
  } else {
    .ridge_pureR(X, Y, lambda, nrand, mode = mode)
  }
}
