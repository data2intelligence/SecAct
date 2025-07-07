#' @title Secreted protein activity inference
#' @description Infer the activity of over 1000 secreted proteins from tumor gene expression profiles.
#' @param Y Gene expression matrix with gene symbol (row) x sample (column).
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomizations in the permutation test, with a default value 1000.
#' @return
#'
#' A list with four items, each is a matrix.
#' beta: regression coefficients
#' se: standard errors of coefficients
#' zscore: beta/se
#' pvalue: statistical significance
#'
#' @rdname SecAct.inference
#' @export
#'
SecAct.inference <- function(Y, SigMat="SecAct", lambda=5e+5, nrand=1000)
{
    if(SigMat=="SecAct")
    {
      Xfile<- file.path(system.file(package = "SecAct"), "extdata/vst_condition_logUMI_cellType.tsv.gz")
      X <- read.table(Xfile,sep="\t",check.names=F)
    }else{
      X <- read.table(SigMat,sep="\t",check.names=F)
    }

    olp <- intersect(row.names(Y),row.names(X))
    X <- as.matrix(X[olp,,drop=F])
    Y <- as.matrix(Y[olp,,drop=F])

    X <- scale(X)
    Y <- scale(Y)

    n <- length(olp)
    p <- ncol(X)
    m <- ncol(Y)

    res <- .C("ridgeReg",
      X=as.double(t(X)),
      Y=as.double(t(Y)),
      as.integer(n),
      as.integer(p),
      as.integer(m),
      as.double(lambda),
      as.double(nrand),
      beta=double(p*m),
      se=double(p*m),
      zscore=double(p*m),
      pvalue=double(p*m)
    )

    beta <- matrix(res$beta,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
    se <- matrix(res$se,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
    zscore <- matrix(res$zscore,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
    pvalue <- matrix(res$pvalue,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))

    beta <- expand_rows(beta)
    se <- expand_rows(se)
    zscore <- expand_rows(zscore)
    pvalue <- expand_rows(pvalue)

    res <- list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)

    res
}


#' @title Secreted protein activity inference
#' @description Infer the signaling activity of 1248 secreted proteins from gene expression profiles.
#' @param inputProfile Gene expression matrix with gene symbol (row) x sample (column).
#' @param inputProfile_control Gene expression matrix with gene symbol (row) x sample (column).
#' @param is.differential A logical indicating whether inputProfile has been differential profiles against to control.
#' @param is.paired A logical indicating whether you want a paired operation of differential profiles between inputProfile and inputProfile_control if samples in inputProfile and inputProfile_control are paired.
#' @param is.singleSampleLevel A logical indicating whether to calculate activity change for each single sample between inputProfile and inputProfile_control. If FALSE, calculate the overall activity change between two phenotypes.
#' @param sigMatrix Secreted protein signature matrix.
#' @param is.group.sig A logical indicating whether to group similar signatures.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @param sigFilter A logical indicating whether filter the secreted protein signatures with the genes from inputProfile.
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
  is.group.sig=TRUE,
  lambda=5e+5,
  nrand=1000,
  sigFilter=FALSE
)
{
  if(class(inputProfile)[1]=="SpaCET")
  {
    stop("Please use 'SecAct.activity.inference.ST'.")
  }
  if(class(inputProfile)[1]=="Seurat")
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

  if(sigMatrix=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/vst_condition_logUMI_cellType.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(sigMatrix,sep="\t",check.names=F)
  }

  if(is.group.sig==TRUE)
  {

    dis <- as.dist(1-cor(X,method="pearson"))
    hc <- hclust(dis,method="complete")

    group_labels <- cutree(hc, h = 1 - 0.9)
    newsig <- data.frame()

    for(j in unique(group_labels))
    {
      geneGroups <- names(group_labels)[group_labels==j]

      newsig[rownames(X),paste0(geneGroups,collapse="|")] <- rowMeans(X[,geneGroups,drop=F])
    }

    X <- newsig
  }


  if(sigFilter==TRUE)
  {
    #X <- X[,colnames(X)%in%row.names(Y)]

    mat <- X
    vec2 <- row.names(Y)

    # Process column names
    keep_idx <- sapply(colnames(mat), function(x) {
      any(strsplit(x, "\\|")[[1]] %in% vec2)
    })

    new_colnames <- sapply(colnames(mat)[keep_idx], function(x) {
      parts <- strsplit(x, "\\|")[[1]]
      matched <- intersect(parts, vec2)
      if (length(matched) > 0) {
        paste(matched, collapse = "|")
      } else {
        NA  # shouldn't happen due to filter
      }
    })

    # Final matrix
    mat_clean <- mat[, keep_idx, drop = FALSE]
    colnames(mat_clean) <- new_colnames

    X <- mat_clean
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

  X <- scale(X)
  Y <- scale(Y)

  n <- length(olp)
  p <- ncol(X)
  m <- ncol(Y)

  res <- .C("ridgeReg",
    X=as.double(t(X)),
    Y=as.double(t(Y)),
    as.integer(n),
    as.integer(p),
    as.integer(m),
    as.double(lambda),
    as.double(nrand),
    beta=double(p*m),
    se=double(p*m),
    zscore=double(p*m),
    pvalue=double(p*m)
  )

  beta <- matrix(res$beta,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  se <- matrix(res$se,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  zscore <- matrix(res$zscore,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  pvalue <- matrix(res$pvalue,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))

  beta <- expand_rows(beta)
  se <- expand_rows(se)
  zscore <- expand_rows(zscore)
  pvalue <- expand_rows(pvalue)

  res <- list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)

  res
}


#' @title Spot activity inference from spatial data
#' @description Calculate secreted protein signaling activity of spots from spatial transcriptomocs data.
#' @param inputProfile A SpaCET object.
#' @param inputProfile_control A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param is.group.sig A logical indicating whether to group similar signatures.
#' @param sigMatrix Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @param sigFilter A logical indicating whether filter the secreted protein signatures with the genes from inputProfile.
#' @return A SpaCET object.
#' @rdname SecAct.activity.inference.ST
#' @export
#'
SecAct.activity.inference.ST <- function(
    inputProfile,
    inputProfile_control = NULL,
    scale.factor = 1e+05,
    sigMatrix="SecAct",
    is.group.sig=TRUE,
    lambda=5e+5,
    nrand=1000,
    sigFilter=FALSE
)
{
  if(!class(inputProfile)[1]=="SpaCET")
  {
    stop("Please input a SpaCET object.")
  }

  # extract count matrix
  expr <- inputProfile@input$counts
  expr <- expr[Matrix::rowSums(expr)>0,]
  rownames(expr) <- transferSymbol(rownames(expr))
  expr <- rm_duplicates_sparse(expr)

  # normalize to TPM
  stats <- Matrix::colSums(expr)
  expr <- sweep_sparse(expr,2,stats,"/")
  expr@x <- expr@x * scale.factor

  # transform to log space
  expr@x <- log2(expr@x + 1)

  if(is.null(inputProfile_control))
  {
    # normalized with the control samples
    expr.diff <- expr - Matrix::rowMeans(expr)

  }else{
    # extract count matrix
    expr_control <- inputProfile_control@input$counts
    expr_control <- expr_control[Matrix::rowSums(expr_control)>0,]
    rownames(expr_control) <- transferSymbol(rownames(expr_control))
    expr_control <- rm_duplicates(expr_control)

    # normalize to TPM
    stats <- Matrix::colSums(expr_control)
    expr_control <- sweep_sparse(expr_control,2,stats,"/")
    expr_control@x <- expr_control@x * scale.factor

    # transform to log space
    expr_control@x <- log2(expr_control@x + 1)

    olp <- intersect(rownames(expr), rownames(expr_control))
    expr.diff <- expr[olp,] - Matrix::rowMeans(expr_control[olp,])
  }

  res <- SecAct.activity.inference(
    inputProfile = expr.diff,
    is.differential = TRUE,
    sigMatrix = sigMatrix,
    is.group.sig = is.group.sig,
    lambda = lambda,
    nrand = nrand,
    sigFilter = sigFilter
  )

  inputProfile @results $SecAct_output $SecretedProteinActivity <- res

  inputProfile
}


#' @title Cell state activity inference from single cell data
#' @description Calculate secreted protein signaling activity of cell states from single cell RNA-Sequencing data.
#' @param inputProfile A Seurat object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param sigMatrix Secreted protein signature matrix.
#' @param is.group.sig A logical indicating whether to group similar signatures.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @param sigFilter A logical indicating whether filter the secreted protein signatures with the genes from inputProfile.
#' @return A Seurat object.
#' @rdname SecAct.activity.inference.scRNAseq
#' @export
#'
SecAct.activity.inference.scRNAseq <- function(
    inputProfile,
    cellType_meta,
    sigMatrix="SecAct",
    is.group.sig=TRUE,
    lambda=5e+5,
    nrand=1000,
    sigFilter=FALSE
)
{
  if(!class(inputProfile)[1]=="Seurat")
  {
    stop("Please input a Seurat object.")
  }

  # extract count matrix
  expr <-  inputProfile@assays$RNA@counts
  cellType_vec <- inputProfile@meta.data[,cellType_meta]

  # generate pseudo bulk
  bulk <- data.frame()
  for(cellType in sort(unique(cellType_vec)))
  {
    bulk[rownames(expr),cellType] <- Matrix::rowSums(expr[,cellType_vec==cellType,drop=F])
  }
  expr <- bulk

  # normalize to TPM
  expr <- sweep(expr, 2, Matrix::colSums(expr), "/") *1e6

  # transform to log space
  expr <- log2(expr + 1)

  # normalized with the control samples
  expr.diff <- expr - rowMeans(expr)

  res <- SecAct.activity.inference(
    inputProfile = expr.diff,
    is.differential = TRUE,
    sigMatrix = sigMatrix,
    is.group.sig = is.group.sig,
    lambda = lambda,
    nrand = nrand,
    sigFilter = sigFilter
  )

  inputProfile @misc $SecAct_output $SecretedProteinActivity <- res

  inputProfile
}
