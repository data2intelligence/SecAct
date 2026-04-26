#' Dispatch ridge+permutation over Y-column batches
#'
#' Memory-efficient version of \code{\link{.ridge_dispatch}}. Processes
#' \code{Y} in column-batches, forwarding to the chosen backend's
#' \code{ridge_batch()} (RidgeCuda > RidgeFast) or to the in-house
#' pure-R batched path. Supports HDF5 Y input and streaming HDF5 output.
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

  if (chosen == "gpu") {
    RidgeCuda::ridge_batch(X = X, Y = Y, lambda = lambda, nrand = nrand,
                           ncores = ncores, rng_method = rng_method,
                           batch_size = batch_size, reader = reader,
                           n_samples = n_samples, output_h5 = output_h5)
  } else if (chosen == "cpu-fast") {
    RidgeFast::ridge_batch(X = X, Y = Y, lambda = lambda, nrand = nrand,
                           ncores = ncores, rng_method = rng_method,
                           batch_size = batch_size, reader = reader,
                           n_samples = n_samples, output_h5 = output_h5)
  } else {
    if (rng_method != "mt19937") {
      stop("rng_method='", rng_method, "' requires RidgeFast or RidgeCuda; ",
           "pure-R supports only 'mt19937'.")
    }
    .ridge_pureR_batch(X, Y, lambda, nrand,
                       batch_size = batch_size,
                       reader = reader, n_samples = n_samples,
                       output_h5 = output_h5)
  }
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
                                output_h5 = NULL) {
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

    res <- .ridge_pureR(X, Y_batch, lambda, nrand)

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
