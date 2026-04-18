# =========================================================
# Pure-R ridge regression with permutation testing
# and backend dispatch to optional accelerators.
#
# API matches Beibei's RidgeRegFast spec:
#   Ridge.Reg.Fast(X, Y, lambda, nrand, ncores, rng_method)
#
# rng_method values:
#   "mt19937"  — canonical GSL-compatible MT19937 (seed 0) with
#                deterministic Fisher-Yates. Bit-identical across
#                backends when ncores=1. Default.
#   "srand"    — platform C stdlib rand() via accelerator backends
#                (not cross-platform reproducible). Pure-R does not
#                support this.
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
#' @param rng_method \code{"mt19937"} only in pure-R. See file header.
#' @return list(beta, se, zscore, pvalue) — each a p x m matrix.
#' @keywords internal
.ridge_pureR <- function(X, Y, lambda, nrand, rng_method = "mt19937") {
  if (rng_method != "mt19937") {
    stop("Pure-R backend only supports rng_method='mt19937'. ",
         "For rng_method='", rng_method, "' install RidgeRegFast or RidgeRegCuda.")
  }
  n <- nrow(Y); p <- ncol(X); m <- ncol(Y)

  A <- crossprod(X) + lambda * diag(p)
  R <- chol(A)
  Rt <- t(R)
  beta <- backsolve(R, forwardsolve(Rt, crossprod(X, Y)))

  perm_table <- .gsl_mt19937_perm_table(n, nrand)

  aver <- matrix(0, p, m)
  aver_sq <- matrix(0, p, m)
  pvalue_count <- matrix(0, p, m)

  for (i in seq_len(nrand)) {
    perm <- perm_table[i, ] + 1L
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
#' Picks GPU (RidgeRegCuda) > CPU-fast (RidgeRegFast) > CPU-pure (this
#' package) based on \code{backend}. All three accelerators share the
#' API \code{Ridge.Reg.Fast(X, Y, lambda, nrand, ncores, rng_method)}.
#' With \code{rng_method="mt19937"} and \code{ncores=1}, results are
#' bit-identical across backends.
#'
#' @keywords internal
.ridge_dispatch <- function(X, Y, lambda, nrand,
                            ncores = 1L,
                            rng_method = "mt19937",
                            backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  chosen <- .resolve_backend(backend)
  message("[SecAct] ridge backend: ", chosen, " (rng_method=", rng_method,
          ", ncores=", ncores, ")")

  if (chosen == "gpu") {
    RidgeRegCuda::Ridge.Reg.Fast(X = X, Y = Y, lambda = lambda, nrand = nrand,
                                 ncores = ncores, rng_method = rng_method)
  } else if (chosen == "cpu-fast") {
    RidgeRegFast::Ridge.Reg.Fast(X = X, Y = Y, lambda = lambda, nrand = nrand,
                                 ncores = ncores, rng_method = rng_method)
  } else {
    .ridge_pureR(X, Y, lambda, nrand, rng_method = rng_method)
  }
}
