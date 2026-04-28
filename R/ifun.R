sweep_sparse <- function(m, margin, stats, fun)
{
  f <- match.fun(fun)

  if(margin==1)
  {
    idx <- m@i + 1
  }else{
    if(inherits(m, "dgCMatrix"))
    {
      idx <- rep(1:m@Dim[2], diff(m@p))
    }else{
      idx <- m@j + 1
    }
  }

  m@x <- f(m@x, stats[idx])
  m
}

transferSymbol <- function(x)
{
  alias2symbol <- read.csv(system.file("extdata", 'NCBI_20251008_gene_result_alias2symbol.csv.gz', package = 'SecAct'),as.is=T)
  alias2symbol[is.na(alias2symbol[,"Alias"]),"Alias"] <- "NA"

  x[x%in%alias2symbol[,1]] <- alias2symbol[
    match(
      x[x%in%alias2symbol[,1]],
      alias2symbol[,1]
    ), 2]

  x
}

rm_duplicates <- function(mat)
{
  gene_count <- table(rownames(mat))
  gene_dupl <- names(gene_count)[gene_count>1]

  if(length(gene_dupl) > 0){
    gene_unique <- names(gene_count)[gene_count==1]
    gene_unique_index <- which(rownames(mat)%in%gene_unique)

    gene_dupl_index <- c()
    for(gene in gene_dupl)
    {
      gene_dupl_index_gene <- which(rownames(mat)%in%gene)
      mat_dupl_gene <- mat[gene_dupl_index_gene,]
      dupl_sum <- Matrix::rowSums(mat_dupl_gene)
      max_flag <- which(dupl_sum==max(dupl_sum))
      gene_dupl_index <- c(gene_dupl_index,gene_dupl_index_gene[max_flag[1]])
    }

    mat <- mat[sort(c(gene_unique_index,gene_dupl_index)),]
  }

  return(mat)
}

scalar1 <- function(x)
{
  x / sqrt(sum(x^2))
}

find_neighbors <- function(coordinate_mat, radius, k = 100)
{
  nn_result <- RANN::nn2(
    coordinate_mat[,c("coordinate_x_um","coordinate_y_um")],
    k=k, searchtype="radius", radius=radius
  )

  neighbor_indices <- nn_result$nn.idx
  neighbor_distances <- nn_result$nn.dists

  i <- rep(1:nrow(neighbor_indices), each=ncol(neighbor_indices))
  j <- as.vector(t(neighbor_indices))
  x <- as.vector(t(neighbor_distances))

  valid <- x <= radius & x > 0
  list(i=i[valid], j=j[valid], x=x[valid], n=nrow(neighbor_indices))
}

calWeights <- function(spotCoordinates, radius, sigma=100, diagAsZero=TRUE)
{
  nn_result <- RANN::nn2(spotCoordinates, searchtype="radius", radius=radius, k=nrow(spotCoordinates))

  neighbor_indices <- nn_result$nn.idx
  neighbor_distances <- nn_result$nn.dists

  i <- rep(1:nrow(neighbor_indices), each=ncol(neighbor_indices)) # row indices (cell index)
  j <- as.vector(t(neighbor_indices))
  x <- as.vector(t(neighbor_distances))

  valid <- x<=radius & x>0
  i <- i[valid]          # Keep only valid indices
  j <- j[valid]          # Valid neighbor indices
  x <- x[valid]          # Valid distances

  # transform distance to weight
  x <- exp( -x^2 / (2*sigma^2) )

  # Create the sparse matrix using the 'i', 'j', and 'x' vectors
  W <- Matrix::sparseMatrix(i=i, j=j, x=x, dims=c(nrow(neighbor_indices), nrow(neighbor_indices)), repr="T")
  rownames(W) <- rownames(spotCoordinates)
  colnames(W) <- rownames(spotCoordinates)

  if(diagAsZero==FALSE) diag(W) <- 1

  W
}

CoxPH_best_separation = function(X, Y, margin)
{
  # part 1: continuous regression
  errflag = F

  coxph.fit = tryCatch(
    survival::coxph(Y~., data=X),
    error = function(e) errflag <<- T,
    warning = function(w) errflag <<- T)

  if(errflag) return (NA)

  n_r = nrow(X)
  n_c = ncol(X)

  arr_result = summary(coxph.fit)$coef[n_c,]

  if(is.na(arr_result["z"])) return (NA)

  # no need to find optimal threshold
  if(is.null(margin)) return (arr_result)

  # part 2: find the optimal threshold
  arr = X[, n_c]
  vthres = sort(arr)

  # these are missing values, not NULL not existing values
  zscore_opt = thres_opt = NA

  for(i in (margin+1):(n_r-margin))
  {
    X[, n_c] = as.numeric(arr >= vthres[i])

    errflag = F
    coxph.fit = tryCatch(
      survival::coxph(Y~., data=X),
      error = function(e) errflag <<- T,
      warning = function(w) errflag <<- T)

    if(errflag) next

    z = summary(coxph.fit)$coef[n_c, "z"]
    if(is.na(z)) next

    if (is.na(zscore_opt)){
      zscore_opt = z
      thres_opt= vthres[i]

    }else if(arr_result['z'] > 0){
      if(z > zscore_opt){
        zscore_opt = z
        thres_opt = vthres[i]
      }

    }else{ # arr_result['z'] <= 0
      if(z < zscore_opt){
        zscore_opt = z
        thres_opt = vthres[i]
      }
    }
  }

  arr_result['thres.opt'] = thres_opt
  arr_result['z.opt'] = zscore_opt

  return (arr_result)
}

load_sig_matrix <- function(sigMatrix, lambda = NULL)
{
  if(sigMatrix == "SecAct")
  {
    Xfile <- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile, sep="\t", check.names=FALSE)
    if(is.null(lambda)) lambda <- 5e+05

  }else if(grepl("SecAct-", sigMatrix, fixed=TRUE)){
    Xfile <- paste0("https://hpc.nih.gov/~Jiang_Lab/SecAct_Package/", sigMatrix, "_filterByPan_ds3_vst.tsv")
    X <- read.table(Xfile, sep="\t", check.names=FALSE)
    if(is.null(lambda)) lambda <- 5e+05

  }else if(sigMatrix == "CytoSig"){
    Xfile <- "https://raw.githubusercontent.com/data2intelligence/CytoSig/refs/heads/master/CytoSig/signature.centroid"
    X <- read.table(Xfile, sep="\t", check.names=FALSE)
    if(is.null(lambda)) lambda <- 10000

  }else{
    X <- read.table(sigMatrix, sep="\t", check.names=FALSE)
  }

  list(X = X, lambda = lambda)
}

cluster_signatures <- function(X,is.group.cor=0.9)
{
  dis <- as.dist(1-cor(X,method="pearson"))
  hc <- hclust(dis,method="complete")

  group_labels <- cutree(hc, h = 1 - is.group.cor)
  newsig <- data.frame()

  for(j in unique(group_labels))
  {
    geneGroups <- names(group_labels)[group_labels==j]

    newsig[rownames(X),paste0(geneGroups,collapse="|")] <- rowMeans(X[,geneGroups,drop=F])
  }

  newsig
}

expand_rows <- function(mat)
{
  new_rows <- lapply(1:nrow(mat), function(i) {
    names <- strsplit(rownames(mat)[i], "\\|")[[1]]
    do.call(rbind, replicate(length(names), mat[i, , drop = FALSE], simplify = FALSE)) |>
      `rownames<-`(names)
  })
  do.call(rbind, new_rows)
}

normalize_log_sparse <- function(expr, scale.factor)
{
  stats <- Matrix::colSums(expr)
  expr <- sweep_sparse(expr, 2, stats, "/")
  expr@x <- expr@x * scale.factor
  expr@x <- log2(expr@x + 1)
  expr
}

extract_ccc_data <- function(data)
{
  if(inherits(data, "SpaCET"))
  {
    ccc <- data@results$SecAct_output$SecretedProteinCCC
  }else if(inherits(data, "Seurat")){
    ccc <- data@misc$SecAct_output$SecretedProteinCCC
  }else{
    stop("Please input a SpaCET or Seurat object.")
  }
  ccc <- cbind(ccc, communication=1)
  ccc <- cbind(ccc, senderReceiver=paste0(ccc[,"sender"], "-", ccc[,"receiver"]))
  ccc
}

build_ccc_matrix <- function(ccc)
{
  mat <- reshape2::acast(ccc[,c("sender","receiver","communication")], sender~receiver, length, value.var="communication")
  cellTypes <- sort(unique(c(rownames(mat), colnames(mat))))
  for(cellType in cellTypes)
  {
    if(cellType %in% rownames(mat) & cellType %in% colnames(mat))
    {
      mat[cellType, cellType] <- NA
    }
  }
  mat
}

extract_seurat_counts <- function(obj)
{
  if(inherits(obj@assays$RNA, "Assay5"))
  {
    counts <- obj@assays$RNA@layers$counts
    colnames(counts) <- rownames(obj@assays$RNA@cells)
    rownames(counts) <- rownames(obj@assays$RNA@features)
  }else{
    counts <- obj@assays$RNA@counts
  }
  counts
}

compute_spatial_correlation <- function(act_new, exp_new, exp_new_aggr)
{
  genes <- rownames(act_new)
  n <- length(genes)
  rs <- numeric(n)
  ps <- numeric(n)
  names(rs) <- names(ps) <- genes

  for(i in seq_len(n))
  {
    gene <- genes[i]
    if(gene %in% rownames(exp_new))
    {
      cor_res <- cor.test(act_new[gene,], exp_new_aggr[gene,], method="spearman")
      rs[i] <- cor_res$estimate
      ps[i] <- cor_res$p.value
    }else{
      rs[i] <- NA
      ps[i] <- NA
    }
  }
  data.frame(r=rs, p=ps, padj=p.adjust(ps, method="BH"))
}


# =========================================================
# Pure R implementation of GSL's MT19937 (Mersenne Twister)
#
# CANONICAL SOURCE for this file. A copy lives in
# RidgeCuda/R/mt19937.R — keep them in sync.
#
# Drift canary: tests/testthat/test-backend-parity.R calls
# .ridge_dispatch with backend="cpu-pure" and backend="cpu-fast"
# and asserts bit-identical outputs. If this file drifts from
# the GSL MT19937 (seed 0) contract the cpu-pure output diverges
# from RidgeFast's C-side build_perm_table and the test fails.
# The RidgeCuda copy is not exercised by SecAct's tests — run
# RidgeCuda/tests/test_gpu_parity.sbatch after editing it.
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
#' RidgeFast (CPU), and RidgeCuda (GPU) so permutation sequences
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
