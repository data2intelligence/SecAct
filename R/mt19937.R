# =========================================================
# Pure R implementation of GSL's MT19937 (Mersenne Twister)
#
# CANONICAL SOURCE for this file. A copy lives in
# RidgeRegCuda/R/mt19937.R — keep them in sync.
#
# Produces output identical to gsl_rng_mt19937 with seed 0.
# Algorithm: Matsumoto & Nishimura (1998) with 2002 init.
# =========================================================

#' @keywords internal
.xor32 <- function(a, b) {
  bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
}

#' @keywords internal
.and32 <- function(a, b) {
  bitwAnd(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwAnd(as.integer(a %% 65536), as.integer(b %% 65536))
}

#' @keywords internal
.or32 <- function(a, b) {
  bitwOr(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwOr(as.integer(a %% 65536), as.integer(b %% 65536))
}

#' @keywords internal
.shr32 <- function(a, k) a %/% (2^k)

#' @keywords internal
.shl32 <- function(a, k) (a * (2^k)) %% 4294967296

#' @keywords internal
.mul32 <- function(a, b) {
  a_lo <- a %% 65536;  a_hi <- a %/% 65536
  b_lo <- b %% 65536;  b_hi <- b %/% 65536
  (a_lo * b_lo + (a_lo * b_hi + a_hi * b_lo) * 65536) %% 4294967296
}

#' @keywords internal
.mt19937_new <- function(seed = 0) {
  N <- 624L
  if (seed == 0) seed <- 4357  # GSL default
  mt <- double(N)
  mt[1] <- seed %% 4294967296
  for (i in 2:N) {
    prev <- mt[i - 1L]
    mt[i] <- (.mul32(1812433253, .xor32(prev, .shr32(prev, 30))) + (i - 1L)) %% 4294967296
  }
  list(mt = mt, mti = N)
}

#' @keywords internal
.mt19937_twist <- function(mt) {
  N <- 624L; M <- 397L
  UPPER <- 2147483648
  LOWER <- 2147483647
  MAGIC <- 2567483615
  mt_old <- mt

  y <- .or32(.and32(mt_old[1:(N - 1L)], UPPER), .and32(mt_old[2:N], LOWER))
  sy <- .xor32(.shr32(y, 1), ifelse(y %% 2 == 1, MAGIC, 0))

  idx1 <- 1L:(N - M)
  mt[idx1] <- .xor32(mt_old[idx1 + M], sy[idx1])

  idx2a <- (N - M + 1L):(2L * (N - M))
  mt[idx2a] <- .xor32(mt[idx2a + M - N], sy[idx2a])

  idx2b <- (2L * (N - M) + 1L):(N - 1L)
  mt[idx2b] <- .xor32(mt[idx2b + M - N], sy[idx2b])

  y <- .or32(.and32(mt_old[N], UPPER), .and32(mt[1L], LOWER))
  mt[N] <- .xor32(mt[M], .xor32(.shr32(y, 1), if (y %% 2 == 1) MAGIC else 0))
  mt
}

#' @keywords internal
.mt19937_generate <- function(state, count) {
  values <- double(count)
  pos <- 1L
  while (pos <= count) {
    if (state$mti >= 624L) {
      state$mt <- .mt19937_twist(state$mt)
      state$mti <- 0L
    }
    avail <- 624L - state$mti
    take  <- min(avail, count - pos + 1L)
    raw   <- state$mt[(state$mti + 1L):(state$mti + take)]

    k <- raw
    k <- .xor32(k, .shr32(k, 11))
    k <- .xor32(k, .and32(.shl32(k, 7), 2636928640))
    k <- .xor32(k, .and32(.shl32(k, 15), 4022730752))
    k <- .xor32(k, .shr32(k, 18))

    values[pos:(pos + take - 1L)] <- k
    state$mti <- state$mti + take
    pos <- pos + take
  }
  list(values = values, state = state)
}

#' Pure R permutation table using GSL-compatible MT19937
#'
#' Generates the same permutation table as GSL's MT19937 (seed 0) with
#' Fisher-Yates shuffle. Used as the canonical RNG across SecAct,
#' RidgeRegFast (CPU), and RidgeRegCuda (GPU) so permutation sequences
#' match bitwise across backends.
#'
#' @param n Number of samples (array length to shuffle).
#' @param nrand Number of permutations.
#' @return Integer matrix (nrand x n), 0-indexed.
#' @keywords internal
.gsl_mt19937_perm_table <- function(n, nrand) {
  MT_MAX <- 4294967295

  total <- as.integer(n - 1L) * nrand
  state <- .mt19937_new(seed = 0)
  gen   <- .mt19937_generate(state, total)
  vals  <- gen$values

  divisors <- MT_MAX %/% (n:2L) + 1
  vals_mat <- matrix(vals, nrow = n - 1L, ncol = nrand)
  j_mat    <- vals_mat %/% divisors + seq_len(n - 1L)
  storage.mode(j_mat) <- "integer"

  arr   <- seq_len(n) - 1L
  table <- matrix(0L, nrow = nrand, ncol = n)

  for (perm in seq_len(nrand)) {
    j_vec <- j_mat[, perm]
    for (i in seq_len(n - 1L)) {
      j <- j_vec[i]
      tmp <- arr[j]; arr[j] <- arr[i]; arr[i] <- tmp
    }
    table[perm, ] <- arr
  }
  table
}
