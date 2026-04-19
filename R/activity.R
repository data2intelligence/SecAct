#' @title Secreted protein activity inference (pure R)
#' @description Infer the activity of over 1000 secreted proteins from gene
#'   expression profiles using a pure-R ridge + permutation kernel. No GSL
#'   or OpenMP required. For large datasets, install \code{RidgeFast}
#'   (CPU accelerator) or \code{RidgeCuda} (GPU accelerator) and use
#'   \code{SecAct.activity.inference} with \code{backend="auto"}.
#' @param Y Gene expression matrix with gene symbol (row) x sample (column).
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomizations in the permutation test.
#' @return list(beta, se, zscore, pvalue).
#' @rdname SecAct.inference.r
#' @export
SecAct.inference.r <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000)
{
  sig <- load_sig_matrix(SigMat, lambda)
  X <- sig$X
  lambda <- sig$lambda

  olp <- intersect(rownames(Y), rownames(X))
  X <- as.matrix(X[olp, , drop = FALSE])
  Y <- as.matrix(Y[olp, , drop = FALSE])

  X <- scale(X)
  Y <- scale(Y)

  res <- .ridge_pureR(X, Y, lambda, nrand)

  for (nm in names(res)) res[[nm]] <- expand_rows(res[[nm]])
  idx <- sort(rownames(res$beta))
  for (nm in names(res)) res[[nm]] <- res[[nm]][idx, , drop = FALSE]

  res
}


#' @title Secreted protein activity inference
#' @description Infer the signaling activity of over 1000 secreted proteins from gene expression profiles.
#' @param inputProfile Gene expression matrix with gene symbol (row) x sample (column).
#' @param inputProfile_control Gene expression matrix with gene symbol (row) x sample (column).
#' @param is.differential A logical flag indicating whether inputProfile has been differential profiles against the control (Default: FALSE).
#' @param is.paired A logical flag indicating whether you want a paired operation of differential profiles between inputProfile and inputProfile_control if samples in inputProfile and inputProfile_control are paired (Default: FALSE).
#' @param is.singleSampleLevel A logical flag indicating whether to calculate activity change for each single sample between inputProfile and inputProfile_control (Default: FALSE). If FALSE, calculate the overall activity change between two phenotypes .
#' @param sigMatrix Secreted protein signature matrix. Could be "SecAct", "CytoSig", "SecAct-Breast", "SecAct-Colorectal", "SecAct-Glioblastoma", "SecAct-Kidney", "SecAct-Liver", "SecAct-Lung-Adeno", "SecAct-Ovarian", "SecAct-Pancreatic", "SecAct-Prostate". SecAct signatures were derived from all cancer ST samples; SecAct-XXX signatures were derived from XXX cancer ST samples.
#' @param is.filter.sig A logical flag indicating whether to filter the secreted protein signatures based on the genes from inputProfile (Default: FALSE). Because some sequencing platforms (e.g., CosMx) cover only a subset of secreted proteins, setting this option to TRUE restricts the activity inference on those proteins.
#' @param is.group.sig A logical flag indicating whether to group similar signatures (Default: TRUE). Many secreted proteins, such as cytokines with similar cell surface receptors and downstream pathways, have cellular effects that appear redundant within a cellular context. When enabled, this option clusters secreted proteins based on Pearson correlations among their composite signatures. The output still reports activity estimates for all secreted proteins prior to clustering. Secreted proteins assigned to the same non-redundant cluster share the same inferred activity.
#' @param is.group.cor A numeric value specifying the correlation cutoff used to define similar signatures (Default: 0.90). When r > 0.90, 1,170 secreted protein signatures are grouped into 657 non-redundant signature groups.
#' @param lambda Penalty factor in the ridge regression. If NULL, lambda will be assigned as 5e+05 or 10000 when sigMatrix = "SecAct" or "CytoSig", respectively.
#' @param nrand Number of randomization in the permutation test, with a default value of 1000.
#' @param ncores Number of threads for accelerator backends (ignored by pure R). Default 1.
#' @param backend One of \code{"auto"}, \code{"gpu"}, \code{"cpu-fast"}, \code{"cpu-pure"}.
#'   \code{"auto"} picks GPU (RidgeCuda) > CPU-fast (RidgeFast) > CPU-pure
#'   depending on what is installed. Default \code{"auto"}.
#' @param rng_method RNG for permutations. \code{"mt19937"} (default) is
#'   GSL-compatible MT19937 seed 0 — bit-identical across backends when
#'   \code{ncores=1}. Accelerators may support \code{"srand"} for faster,
#'   non-reproducible runs. Pure-R supports only \code{"mt19937"}.
#' @param batch_size Optional positive integer. When supplied,
#'   permutation inference is performed over column-batches of \code{Y}
#'   via the backend's \code{ridge_batch} (memory-efficient path for
#'   large sample counts). Default \code{NULL} (no batching).
#' @param output_h5 Optional path to an HDF5 file. When supplied
#'   (requires \code{batch_size} and the \pkg{rhdf5} package), the four
#'   result matrices are streamed to HDF5 datasets
#'   \code{beta/se/zscore/pvalue} and the function returns metadata in
#'   place of the matrices. Use when even the result matrices do not fit
#'   in RAM. Not compatible with \code{is.group.sig=TRUE}.
#' @return
#'
#' A list with four items, each is a matrix.
#' beta: regression coefficients
#' se: standard errors of coefficients
#' zscore: beta/se
#' pvalue: statistical significance
#'
#' @rdname SecAct.activity.inference
#' @export
#'
SecAct.activity.inference <- function(
  inputProfile,
  inputProfile_control=NULL,
  is.differential=FALSE,
  is.paired=FALSE,
  is.singleSampleLevel=FALSE,
  sigMatrix="SecAct",
  is.filter.sig=FALSE,
  is.group.sig=TRUE,
  is.group.cor=0.9,
  lambda=5e+05,
  nrand=1000,
  ncores=1L,
  backend="auto",
  rng_method="mt19937",
  batch_size=NULL,
  output_h5=NULL
)
{
  if(inherits(inputProfile, "SpaCET"))
  {
    stop("Please use 'SecAct.activity.inference.ST'.")
  }
  if(inherits(inputProfile, "Seurat"))
  {
    stop("Please use 'SecAct.activity.inference.scRNAseq'.")
  }

  if(is.differential)
  {
    Y <- inputProfile
    if(ncol(Y)==1) colnames(Y) <- "Change"

  }else{
    if(is.null(inputProfile_control))
    {
      if(ncol(inputProfile)==1) stop("Please include at least two samples in 'inputProfile' or set 'inputProfile_control'.")
      Y <- inputProfile - rowMeans(inputProfile)

    }else{
      if(is.paired)
      {
        Y <- inputProfile - inputProfile_control[,colnames(inputProfile),drop=FALSE]
      }else{
        Y <- inputProfile - rowMeans(inputProfile_control)
      }

      if(is.singleSampleLevel==FALSE)
      {
        Y <- matrix(rowMeans(Y), ncol=1, dimnames=list(rownames(Y), "Change"))
      }

    }
  }

  sig <- load_sig_matrix(sigMatrix, lambda)
  X <- sig$X
  lambda <- sig$lambda

  if(is.filter.sig==TRUE)
  {
    X <- X[,colnames(X)%in%rownames(Y)]
  }

  if(is.group.sig==TRUE)
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

    X <- newsig
  }

  olp <- intersect(rownames(Y),rownames(X))

  if(length(olp)<2) stop("The overlapped genes between your expression matrix and our signature matrix are too few!")

  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

  X <- scale(X)
  Y <- scale(Y)

  if (!is.null(batch_size)) {
    if (!is.null(output_h5) && isTRUE(is.group.sig)) {
      stop("output_h5 cannot be combined with is.group.sig=TRUE (group ",
           "expansion requires matrices in memory).")
    }
    res <- .ridge_batch_dispatch(X, Y, lambda, nrand,
                                 ncores = ncores, rng_method = rng_method,
                                 backend = backend,
                                 batch_size = as.integer(batch_size),
                                 output_h5 = output_h5)
    if (!is.null(output_h5)) return(invisible(res))
  } else {
    res <- .ridge_dispatch(X, Y, lambda, nrand,
                           ncores = ncores, rng_method = rng_method, backend = backend)
  }

  if(is.group.sig==TRUE)
  {
    for(nm in names(res)) res[[nm]] <- expand_rows(res[[nm]])
    idx <- sort(rownames(res$beta))
    for(nm in names(res)) res[[nm]] <- res[[nm]][idx,,drop=F]
  }

  res
}


#' @title Spot activity inference from spatial data
#' @description Calculate secreted protein signaling activity of spots from spatial transcriptomocs data.
#' @param inputProfile A SpaCET object.
#' @param inputProfile_control A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param sigMatrix Secreted protein signature matrix. Could be "SecAct", "CytoSig", "SecAct-Breast", "SecAct-Colorectal", "SecAct-Glioblastoma", "SecAct-Kidney", "SecAct-Liver", "SecAct-Lung-Adeno", "SecAct-Ovarian", "SecAct-Pancreatic", "SecAct-Prostate". SecAct signatures were derived from all cancer ST samples; SecAct-XXX signatures were derived from XXX cancer ST samples.
#' @param is.filter.sig A logical flag indicating whether to filter the secreted protein signatures based on the genes from inputProfile (Default: FALSE). Because some sequencing platforms (e.g., CosMx) cover only a subset of secreted proteins, setting this option to TRUE restricts the activity inference on those proteins.
#' @param is.group.sig A logical flag indicating whether to group similar signatures (Default: TRUE). Many secreted proteins, such as cytokines with similar cell surface receptors and downstream pathways, have cellular effects that appear redundant within a cellular context. When enabled, this option clusters secreted proteins based on Pearson correlations among their composite signatures. The output still reports activity estimates for all secreted proteins prior to clustering. Secreted proteins assigned to the same non-redundant cluster share the same inferred activity.
#' @param is.group.cor A numeric value specifying the correlation cutoff used to define similar signatures (Default: 0.90). When r > 0.90, 1,170 secreted protein signatures are grouped into 657 non-redundant signature groups.
#' @param lambda Penalty factor in the ridge regression. If NULL, lambda will be assigned as 5e+05 or 10000 when sigMatrix = "SecAct" or "CytoSig", respectively.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @return A SpaCET object.
#' @rdname SecAct.activity.inference.ST
#' @export
#'
SecAct.activity.inference.ST <- function(
    inputProfile,
    inputProfile_control = NULL,
    scale.factor = 1e+05,
    sigMatrix="SecAct",
    is.filter.sig=FALSE,
    is.group.sig=TRUE,
    is.group.cor=0.9,
    lambda=5e+05,
    nrand=1000,
    ncores=1L,
    backend="auto",
    rng_method="mt19937"
)
{
  if(!inherits(inputProfile, "SpaCET"))
  {
    stop("Please input a SpaCET object.")
  }

  # extract count matrix
  expr <- inputProfile@input$counts
  expr <- expr[Matrix::rowSums(expr)>0,]
  rownames(expr) <- transferSymbol(rownames(expr))
  expr <- rm_duplicates(expr)

  expr <- normalize_log_sparse(expr, scale.factor)

  if(is.null(inputProfile_control))
  {
    expr.diff <- expr - Matrix::rowMeans(expr)

  }else{
    expr_control <- inputProfile_control@input$counts
    expr_control <- expr_control[Matrix::rowSums(expr_control)>0,]
    rownames(expr_control) <- transferSymbol(rownames(expr_control))
    expr_control <- rm_duplicates(expr_control)

    expr_control <- normalize_log_sparse(expr_control, scale.factor)

    olp <- intersect(rownames(expr), rownames(expr_control))
    expr.diff <- expr[olp,] - Matrix::rowMeans(expr_control[olp,])
  }

  res <- SecAct.activity.inference(
    inputProfile = expr.diff,
    is.differential = TRUE,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    backend = backend,
    rng_method = rng_method
  )

  inputProfile @results $SecAct_output $SecretedProteinActivity <- res

  inputProfile
}


#' @title Cell state activity inference from single cell data
#' @description Calculate secreted protein signaling activity of cell states from single cell RNA-Sequencing data.
#' @param inputProfile A Seurat object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param is.singleCellLevel A logical flag indicating whether to calculate for each single cell (Default: FALSE).
#' @param sigMatrix Secreted protein signature matrix. Could be "SecAct", "CytoSig", "SecAct-Breast", "SecAct-Colorectal", "SecAct-Glioblastoma", "SecAct-Kidney", "SecAct-Liver", "SecAct-Lung-Adeno", "SecAct-Ovarian", "SecAct-Pancreatic", "SecAct-Prostate". SecAct signatures were derived from all cancer ST samples; SecAct-XXX signatures were derived from XXX cancer ST samples.
#' @param is.filter.sig A logical flag indicating whether to filter the secreted protein signatures based on the genes from inputProfile (Default: FALSE). Because some sequencing platforms (e.g., CosMx) cover only a subset of secreted proteins, setting this option to TRUE restricts the activity inference on those proteins.
#' @param is.group.sig A logical flag indicating whether to group similar signatures (Default: TRUE). Many secreted proteins, such as cytokines with similar cell surface receptors and downstream pathways, have cellular effects that appear redundant within a cellular context. When enabled, this option clusters secreted proteins based on Pearson correlations among their composite signatures. The output still reports activity estimates for all secreted proteins prior to clustering. Secreted proteins assigned to the same non-redundant cluster share the same inferred activity.
#' @param is.group.cor A numeric value specifying the correlation cutoff used to define similar signatures (Default: 0.90). When r > 0.90, 1,170 secreted protein signatures are grouped into 657 non-redundant signature groups.
#' @param lambda Penalty factor in the ridge regression. If NULL, lambda will be assigned as 5e+05 or 10000 when sigMatrix = "SecAct" or "CytoSig", respectively.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @return A Seurat object.
#' @rdname SecAct.activity.inference.scRNAseq
#' @export
#'
SecAct.activity.inference.scRNAseq <- function(
    inputProfile,
    cellType_meta,
    is.singleCellLevel=FALSE,
    sigMatrix="SecAct",
    is.filter.sig=FALSE,
    is.group.sig=TRUE,
    is.group.cor=0.9,
    lambda=5e+05,
    nrand=1000,
    ncores=1L,
    backend="auto",
    rng_method="mt19937"
)
{
  if(!inherits(inputProfile, "Seurat"))
  {
    stop("Please input a Seurat object.")
  }

  counts <- extract_seurat_counts(inputProfile)

  rownames(counts) <- transferSymbol(rownames(counts))
  counts <- rm_duplicates(counts)

  if(is.singleCellLevel==FALSE)
  {
    cellType_vec <- inputProfile@meta.data[,cellType_meta]

    # generate pseudo bulk
    expr <- data.frame()
    for(cellType in sort(unique(cellType_vec)))
    {
      expr[rownames(counts),cellType] <- Matrix::rowSums(counts[,cellType_vec==cellType,drop=F])
    }

    # normalize to TPM
    expr <- sweep(expr, 2, Matrix::colSums(expr), "/") *1e6

  }else{
    expr <- sweep(counts, 2, Matrix::colSums(counts), "/") *1e5

  }
  rm(counts);gc()

  # transform to log space
  expr <- log2(expr + 1)

  # normalized with the control samples
  expr.diff <- expr - Matrix::rowMeans(expr)

  rm(expr);gc()

  inputProfile @misc $SecAct_output $SecretedProteinActivity <-
    SecAct.activity.inference(
    inputProfile = expr.diff,
    is.differential = TRUE,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    backend = backend,
    rng_method = rng_method
  )

  inputProfile
}
