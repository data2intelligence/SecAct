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
      Xfile<- file.path(system.file(package = "SecAct"), "extdata/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv.gz")
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
    colnames(Y) <- "Change"
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
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(sigMatrix,sep="\t",check.names=F)
  }

  if(sigFilter==TRUE)
  {
    X <- X[,colnames(X)%in%row.names(Y)]
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

  list(
    beta = matrix(res$beta,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y))),
    se = matrix(res$se,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y))),
    zscore = matrix(res$zscore,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y))),
    pvalue = matrix(res$pvalue,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  )

}


#' @title Spot activity inference from spatial data
#' @description Calculate secreted protein signaling activity of spots from spatial transcriptomocs data.
#' @param inputProfile A SpaCET object.
#' @param inputProfile_control A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
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
    lambda = lambda,
    nrand = nrand,
    sigFilter = sigFilter
  )

  inputProfile @misc $SecAct_output $SecretedProteinActivity <- res

  inputProfile
}


#' @title Cell type activity inference between two conditions from single cell data
#' @description Calculate the changes in secreted protein signaling activity for each cell type between two conditions from single cell RNA-Sequencing data.
#' @param inputProfile A Seurat object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param condition_meta Column name in meta data that includes condition information.
#' @param conditionCase Case condition.
#' @param conditionControl Control condition.
#' @param sigMatrix Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @param sigFilter A logical indicating whether filter the secreted protein signatures with the genes from inputProfile.
#' @return A Seurat object.
#' @rdname SecAct.activity.inference.scRNAseq.2
#' @export
#'
SecAct.activity.inference.scRNAseq.2 <- function(
    inputProfile,
    cellType_meta,
    condition_meta,
    conditionCase,
    conditionControl,
    sigMatrix="SecAct",
    lambda=5e+5,
    nrand=1000,
    sigFilter=FALSE
)
{
  if(!class(inputProfile)[1]=="Seurat")
  {
    stop("Please input a Seurat object.")
  }

  counts <-  inputProfile@assays$RNA@counts
  rownames(counts) <- transferSymbol(rownames(counts))
  counts <- rm_duplicates_sparse(counts)

  meta <- inputProfile@meta.data

  cellTypes <- intersect(
    meta[meta[,condition_meta]==conditionCase,cellType_meta],
    meta[meta[,condition_meta]==conditionControl,cellType_meta]
  )

  print("Step 1: calculating changes in secreted protein activity.")

  bulk.diff <- data.frame()
  for(cellType in cellTypes)
  {
      expr <- counts[,meta[,condition_meta]==conditionCase&meta[,cellType_meta]==cellType]
      expr <- Matrix::rowSums(expr)
      expr <- matrix(expr,ncol=1,dimnames = list(names(expr),"bulk"))

      # normalize to TPM
      expr <- t(t(expr)*1e6/colSums(expr))

      # transform to log space
      expr <- log2(expr + 1)

      expr_case <- expr


      expr <- counts[,meta[,condition_meta]==conditionControl&meta[,cellType_meta]==cellType]
      expr <- Matrix::rowSums(expr)
      expr <- matrix(expr,ncol=1,dimnames = list(names(expr),"bulk"))

      # normalize to TPM
      expr <- t(t(expr)*1e6/colSums(expr))

      # transform to log space
      expr <- log2(expr + 1)

      expr_control <- expr


      # normalized with the control samples
      bulk.diff[rownames(expr_case),cellType] <- expr_case - expr_control
  }

  inputProfile @misc $SecAct_output $SecretedProteinActivity <- SecAct.activity.inference(bulk.diff, is.differential = TRUE, sigMatrix = sigMatrix)

  inputProfile
}
