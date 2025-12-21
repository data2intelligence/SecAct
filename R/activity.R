#' @title Secreted protein activity inference
#' @description Infer the activity of over 1000 secreted proteins from tumor gene expression profiles.
#' @param Y Gene expression matrix with gene symbol (row) x sample (column).
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @return
#'
#' A list with four items, each is a matrix.
#' beta: regression coefficients
#' se: standard errors of coefficients
#' zscore: beta/se
#' pvalue: statistical significance
#'
#' @rdname SecAct.inference.gsl
#' @export
#'
SecAct.inference.gsl <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
{
  if(SigMat=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
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
#' @rdname SecAct.inference.r
#' @export
#'
SecAct.inference.r <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
{
  if(SigMat=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

  X <- scale(X)
  Y <- scale(Y)

  n <- nrow(Y)
  p <- ncol(X)
  m <- ncol(Y)

  A <- crossprod(X) + lambda * diag(p)  # SPD
  R <- chol(A)                          # A = R'R
  beta <- backsolve(R, forwardsolve(t(R), crossprod(X, Y)))

  for(i in 1:nrand)
  {
    set.seed(i)
    beta_rand <- backsolve(R, forwardsolve(t(R), crossprod(X, Y[sample.int(n),])))

    if(i==1)
    {
      aver <- beta_rand
      aver_sq <- beta_rand^2
      pvalue <- abs(beta_rand)>=abs(beta)
    }else{
      aver <- aver + beta_rand
      aver_sq <- aver_sq + beta_rand^2
      pvalue <- pvalue + (abs(beta_rand)>=abs(beta))
    }
  }

  aver <- aver/nrand
  aver_sq <- aver_sq/nrand
  aver_sq <- sqrt(aver_sq - aver*aver)

  zscore <- (beta-aver)/aver_sq

  pvalue <- (pvalue+1)/(nrand+1)

  rownames(beta) <- colnames(X)
  colnames(beta) <- colnames(Y)

  rownames(aver_sq) <- colnames(X)
  colnames(aver_sq) <- colnames(Y)

  rownames(zscore) <- colnames(X)
  colnames(zscore) <- colnames(Y)

  rownames(pvalue) <- colnames(X)
  colnames(pvalue) <- colnames(Y)

  beta <- expand_rows(beta)
  aver_sq <- expand_rows(aver_sq)
  zscore <- expand_rows(zscore)
  pvalue <- expand_rows(pvalue)

  beta <- beta[sort(rownames(beta)),,drop=F]
  aver_sq <- aver_sq[sort(rownames(beta)),,drop=F]
  zscore <- zscore[sort(rownames(beta)),,drop=F]
  pvalue <- pvalue[sort(rownames(beta)),,drop=F]

  res <- list(beta=beta, se=aver_sq, zscore=zscore, pvalue=pvalue)

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
#' @param is.group.sig A logical indicating whether group similar signatures.
#' @param is.group.cor Correlation cutoff of similar signatures.
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
  is.group.cor=0.9,
  lambda=5e+05,
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
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(sigMatrix,sep="\t",check.names=F)
  }

  if(sigFilter==TRUE)
  {
    X <- X[,colnames(X)%in%row.names(Y)]
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

  olp <- intersect(row.names(Y),row.names(X))

  if(length(olp)<2) stop("The overlapped genes between your expression matrix and our signature matrix are too few!")

  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

  X <- scale(X)
  Y <- scale(Y)

  n <- nrow(Y)
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

  if(is.group.sig==TRUE)
  {
    beta <- expand_rows(beta)
    se <- expand_rows(se)
    zscore <- expand_rows(zscore)
    pvalue <- expand_rows(pvalue)

    beta <- beta[sort(rownames(beta)),,drop=F]
    se <- se[sort(rownames(beta)),,drop=F]
    zscore <- zscore[sort(rownames(beta)),,drop=F]
    pvalue <- pvalue[sort(rownames(beta)),,drop=F]
  }

  res <- list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)

  res
}


#' @title Spot activity inference from spatial data
#' @description Calculate secreted protein signaling activity of spots from spatial transcriptomocs data.
#' @param inputProfile A SpaCET object.
#' @param inputProfile_control A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param sigMatrix Secreted protein signature matrix.
#' @param is.group.sig A logical indicating whether to group similar signatures.
#' @param is.group.cor Correlation cutoff of similar signatures.
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
    is.group.cor=0.9,
    lambda=5e+05,
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
  expr <- rm_duplicates(expr)

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
    is.group.cor = is.group.cor,
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
#' @param is.singleCellLevel A logical indicating whether to calculate for each single cell.
#' @param is.group.sig A logical indicating whether to group similar signatures.
#' @param is.group.cor Correlation cutoff of similar signatures.
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
    is.singleCellLevel=FALSE,
    is.group.sig=TRUE,
    is.group.cor=0.9,
    lambda=5e+05,
    nrand=1000,
    sigFilter=FALSE
)
{
  if(!class(inputProfile)[1]=="Seurat")
  {
    stop("Please input a Seurat object.")
  }

  if(class(inputProfile@assays$RNA)=="Assay5")
  {
    counts <- inputProfile@assays$RNA@layers$counts
    colnames(counts) <- rownames(inputProfile@assays$RNA@cells)
    rownames(counts) <- rownames(inputProfile@assays$RNA@features)
  }else{
    counts <-  inputProfile@assays$RNA@counts
  }

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
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    sigFilter = sigFilter
  )

  inputProfile
}
