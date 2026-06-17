#' Pure-R ridge + permutation kernel
#'
#' Solves \eqn{\beta = (X'X + \lambda I)^{-1} X'Y} via Cholesky, then
#' estimates null via Y-row permutations. Precomputes the projection
#' \eqn{T = (X'X + \lambda I)^{-1} X'} once so each permutation is a
#' single p-by-m dgemm instead of two triangular solves plus a crossprod.
#'
#' @param X Numeric matrix, n x p, column-scaled.
#' @param Y Numeric matrix, n x m, column-scaled.
#' @param lambda Ridge penalty.
#' @param nrand Number of permutations.
#' @return list(beta, se, zscore, pvalue) — each a p x m matrix.
#' @keywords internal
.ridge_pureR <- function(X, Y, lambda, nrand) {
  n <- nrow(Y); p <- ncol(X); m <- ncol(Y)

  A <- crossprod(X) + lambda * diag(p)
  R <- chol(A)
  Rt <- t(R)
  T_proj <- backsolve(R, forwardsolve(Rt, t(X)))   # p x n
  beta <- T_proj %*% Y

  perm_table <- .gsl_mt19937_perm_table(n, nrand)

  aver <- matrix(0, p, m)
  aver_sq <- matrix(0, p, m)
  pvalue_count <- matrix(0, p, m)

  for (i in seq_len(nrand)) {
    perm <- perm_table[i, ] + 1L
    beta_rand <- T_proj %*% Y[perm, , drop = FALSE]
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

#' Resolve backend based on installed FlashReg and GPU availability
#'
#' Legacy backend names ("gpu", "cpu-fast", "cpu-pure") are still
#' accepted and remap onto the corresponding FlashReg backends so
#' downstream code keeps working without changes.
#' @keywords internal
.resolve_backend <- function(backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  backend <- match.arg(backend)
  if (backend != "auto") return(backend)

  if (requireNamespace("FlashReg", quietly = TRUE)) {
    gpu_ok <- tryCatch(FlashReg::cuda_available(), error = function(e) FALSE)
    if (isTRUE(gpu_ok)) return("gpu")
    return("cpu-fast")
  }
  "cpu-pure"
}

# Map SecAct's legacy backend name onto FlashReg's backend name.
.secact_to_flashreg_backend <- function(chosen) {
  switch(chosen,
         gpu        = "cuda_native",
         `cpu-fast` = "omp",
         `cpu-pure` = "pure_r",
         stop("internal: unknown SecAct backend label '", chosen, "'"))
}

#' Dispatch ridge+permutation call to FlashReg or the in-house pure-R loop
#'
#' Picks GPU (\code{FlashReg::ridge(backend="cuda_native")}) > CPU-fast
#' (\code{FlashReg::ridge(backend="omp")}) > CPU-pure (this package's
#' \code{.ridge_pureR}). FlashReg replaces the historical pair of
#' optional packages (RidgeFast for CPU, RidgeCuda for GPU); legacy
#' backend names are kept as aliases for backward compatibility.
#'
#' With \code{rng_method="mt19937"} and \code{ncores=1}, FlashReg's
#' backends are bit-identical to the pure-R loop.
#'
#' @keywords internal
.ridge_dispatch <- function(X, Y, lambda, nrand,
                            ncores = 1L,
                            rng_method = "mt19937",
                            backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  rng_method <- match.arg(tolower(rng_method), c("mt19937", "srand"))
  chosen <- .resolve_backend(backend)
  if (isTRUE(getOption("SecAct.verbose", FALSE))) {
    message("[SecAct] ridge backend: ", chosen, " (rng_method=", rng_method,
            ", ncores=", ncores, ")")
  }

  if (chosen %in% c("gpu", "cpu-fast")) {
    if (!requireNamespace("FlashReg", quietly = TRUE)) {
      stop("backend='", chosen, "' requires the FlashReg package. ",
           "Install it from https://github.com/data2intelligence/FlashReg ",
           "or set backend='cpu-pure' to use the in-house pure-R loop.")
    }
    FlashReg::ridge(X = X, Y = Y, lambda = lambda, nrand = nrand,
                    backend = .secact_to_flashreg_backend(chosen),
                    ncores = ncores, rng_method = rng_method)
  } else {
    if (rng_method != "mt19937") {
      stop("rng_method='", rng_method, "' requires FlashReg; ",
           "pure-R supports only 'mt19937'.")
    }
    .ridge_pureR(X, Y, lambda, nrand)
  }
}








#' Dispatch ridge+permutation over Y-column batches
#'
#' Memory-efficient version of \code{\link{.ridge_dispatch}}. Processes
#' \code{Y} in column-batches, calling the resolved backend on each
#' batch. With FlashReg installed, the per-batch kernel is the
#' \code{FlashReg::ridge()} dispatcher (GPU or C+OpenMP CPU);
#' otherwise the pure-R loop is used. The same column-batching code
#' handles HDF5 Y input and streaming HDF5 output for all backends.
#'
#' @keywords internal
.ridge_batch_dispatch <- function(X, Y, lambda, nrand,
                                  ncores = 1L,
                                  rng_method = "mt19937",
                                  backend = c("auto", "gpu", "cpu-fast", "cpu-pure"),
                                  batch_size = 5000L,
                                  reader = NULL, n_samples = NULL,
                                  output_h5 = NULL) {
  rng_method <- match.arg(tolower(rng_method), c("mt19937", "srand"))
  chosen <- .resolve_backend(backend)
  if (isTRUE(getOption("SecAct.verbose", FALSE))) {
    message("[SecAct] ridge_batch backend: ", chosen,
            " (batch_size=", batch_size, ", rng_method=", rng_method, ")")
  }

  if (chosen %in% c("gpu", "cpu-fast")) {
    if (!requireNamespace("FlashReg", quietly = TRUE)) {
      stop("backend='", chosen, "' requires the FlashReg package. ",
           "Install it from https://github.com/data2intelligence/FlashReg ",
           "or set backend='cpu-pure' to use the in-house pure-R loop.")
    }
    flashreg_backend <- .secact_to_flashreg_backend(chosen)
    kernel_fn <- function(X, Y_batch, lambda, nrand) {
      FlashReg::ridge(X = X, Y = Y_batch, lambda = lambda, nrand = nrand,
                      backend = flashreg_backend,
                      ncores = ncores, rng_method = rng_method)
    }
  } else {
    if (rng_method != "mt19937") {
      stop("rng_method='", rng_method, "' requires FlashReg; ",
           "pure-R supports only 'mt19937'.")
    }
    kernel_fn <- .ridge_pureR
  }

  .ridge_pureR_batch(X, Y, lambda, nrand,
                     batch_size = batch_size,
                     reader = reader, n_samples = n_samples,
                     output_h5 = output_h5,
                     kernel_fn = kernel_fn)
}

#' Pure-R batched ridge+permutation
#'
#' Column-batched wrapper around \code{.ridge_pureR}. Supports in-memory,
#' HDF5-path, or user-reader \code{Y} sources and optional HDF5 output.
#' Results are bit-identical to a single \code{.ridge_pureR} call since
#' statistics are per-column and \code{T} is re-derived deterministically
#' per batch from the same \code{X}.
#'
#' @keywords internal
.ridge_pureR_batch <- function(X, Y, lambda, nrand,
                               batch_size = 5000L,
                               reader = NULL, n_samples = NULL,
                               output_h5 = NULL,
                               kernel_fn = NULL) {
  if (is.null(kernel_fn)) kernel_fn <- .ridge_pureR
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  sig_names <- colnames(X)
  if (is.null(sig_names)) sig_names <- paste0("sig", seq_len(p))

  batch_size <- as.integer(batch_size)
  if (is.na(batch_size) || batch_size < 1L) stop("batch_size must be >= 1.")

  src <- .secact_resolve_y_source(Y, reader, n_samples, n)
  m <- src$m
  samp_names <- src$samp_names
  reader_fn <- src$reader_fn

  h5_out <- !is.null(output_h5)
  if (h5_out) {
    if (!requireNamespace("rhdf5", quietly = TRUE))
      stop("'output_h5' requires the 'rhdf5' package.")
    if (file.exists(output_h5)) file.remove(output_h5)
    rhdf5::h5createFile(output_h5)
    chunk_m <- min(batch_size, m)
    for (nm in c("beta", "se", "zscore", "pvalue")) {
      rhdf5::h5createDataset(output_h5, nm, dims = c(p, m),
                             storage.mode = "double",
                             chunk = c(p, chunk_m))
    }
    rhdf5::h5write(sig_names, output_h5, "signature_names")
    rhdf5::h5write(samp_names, output_h5, "sample_names")
  } else {
    out_beta   <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
    out_se     <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
    out_zscore <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
    out_pvalue <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
  }

  num_batches <- as.integer(ceiling(m / batch_size))
  for (b in seq_len(num_batches)) {
    s <- (b - 1L) * batch_size + 1L
    e <- min(as.integer(b * batch_size), m)
    Y_batch <- reader_fn(s, e)
    if (!is.matrix(Y_batch)) Y_batch <- as.matrix(Y_batch)
    storage.mode(Y_batch) <- "double"
    if (is.null(colnames(Y_batch))) colnames(Y_batch) <- samp_names[s:e]

    res <- kernel_fn(X, Y_batch, lambda, nrand)

    if (h5_out) {
      for (nm in c("beta", "se", "zscore", "pvalue")) {
        rhdf5::h5write(res[[nm]], output_h5, nm, index = list(NULL, s:e))
      }
    } else {
      out_beta[, s:e]   <- res$beta
      out_se[, s:e]     <- res$se
      out_zscore[, s:e] <- res$zscore
      out_pvalue[, s:e] <- res$pvalue
    }
  }

  if (h5_out) {
    rhdf5::h5closeAll()
    invisible(list(path = output_h5, p = p, m = m, num_batches = num_batches))
  } else {
    list(beta = out_beta, se = out_se, zscore = out_zscore, pvalue = out_pvalue)
  }
}

.secact_resolve_y_source <- function(Y, reader, n_samples, n) {
  if (is.character(Y) && length(Y) == 1L) {
    if (!requireNamespace("rhdf5", quietly = TRUE))
      stop("Reading Y from an HDF5 file requires the 'rhdf5' package.")
    info <- rhdf5::h5ls(Y)
    y_row <- info[info$name == "Y", ]
    if (nrow(y_row) == 0L) stop("HDF5 file '", Y, "' must contain a dataset named 'Y'.")
    dims <- as.integer(strsplit(y_row$dim[1], " x ", fixed = TRUE)[[1]])
    if (length(dims) != 2L) stop("'Y' dataset must be 2-dimensional.")
    if (dims[1] != n) stop(sprintf("HDF5 Y rows (%d) do not match nrow(X) (%d).", dims[1], n))
    m <- dims[2]
    samp_names <- if ("sample_names" %in% info$name) {
      as.character(rhdf5::h5read(Y, "sample_names"))
    } else paste0("s", seq_len(m))
    path <- Y
    list(m = m, samp_names = samp_names,
         reader_fn = function(s, e) {
           rhdf5::h5read(path, "Y", index = list(NULL, s:e))
         })
  } else if (!is.null(reader)) {
    if (is.null(n_samples)) stop("'reader' requires 'n_samples'.")
    m <- as.integer(n_samples)
    list(m = m, samp_names = paste0("s", seq_len(m)), reader_fn = reader)
  } else {
    if (!is.matrix(Y)) Y <- as.matrix(Y)
    if (nrow(Y) != n) stop(sprintf("nrow(Y) (%d) != nrow(X) (%d).", nrow(Y), n))
    m <- ncol(Y)
    samp_names <- if (is.null(colnames(Y))) paste0("s", seq_len(m)) else colnames(Y)
    list(m = m, samp_names = samp_names,
         reader_fn = function(s, e) Y[, s:e, drop = FALSE])
  }
}
