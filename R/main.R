#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Y PARAM_DESCRIPTION
#' @param SigMat PARAM_DESCRIPTION
#' @param lambda PARAM_DESCRIPTION
#' @param nrand PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#'
#' Yfile <- file.path(system.file(package = "SecAct"), "extdata/IFNG_GSE100093.diff")
#' Y <- read.table(Yfile,sep="\t",check.names=F)
#' res <- SecAct.inference(Y, SigMat=NULL, lambda=10000, nrand=1000)
#' names(res)
#' head(res$zscore)
#'
#' @rdname SecAct.inference
#' @export
#'
SecAct.inference <- function(Y, SigMat=NULL, lambda=10000, nrand=1000)
{
  Y_type <- "matrix"
  if(class(Y)=="SpaCET")
  {
    Y_type <- "SpaCET"
    SpaCET_obj <- Y

    Y <- SpaCET_obj@input$counts

    Y <- Matrix::t(Matrix::t(Y)*1e5/Matrix::colSums(Y))
    Y@x <- log2(Y@x + 1)
    Y <- Y - Matrix::rowMeans(Y)
  }

  if(is.null(SigMat))
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/signature.centroid")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,])
  Y <- as.matrix(Y[olp,])

  #write.table(X,"/Users/rub2/Workspace/GitHub/SecAct/inst/extdata/X",quote=F,sep="\t")
  #write.table(Y,"/Users/rub2/Workspace/GitHub/SecAct/inst/extdata/Y",quote=F,sep="\t")

  X <- scale(X)
  Y <- scale(Y)

  #write.table(X,"/Users/rub2/Workspace/GitHub/SecAct/inst/extdata/X_z",quote=F,sep="\t")
  #write.table(Y,"/Users/rub2/Workspace/GitHub/SecAct/inst/extdata/Y_z",quote=F,sep="\t")

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

  if(Y_type=="SpaCET")
  {
    SpaCET_obj@results$SecAct_res <- res
    SpaCET_obj
  }else{
    res
  }

}


calWeights <- function(SpotIDs, r=3, diag0=TRUE)
{
  d <- matrix(Inf,ncol=length(SpotIDs),nrow=length(SpotIDs))
  colnames(d) <- SpotIDs
  rownames(d) <- SpotIDs

  for(i in 1:ncol(d))
  {
    xy <- rownames(d)[i]
    x <- as.numeric(unlist(strsplit(xy,"x")))[1]
    y <- as.numeric(unlist(strsplit(xy,"x")))[2]
    xm <- r-1
    ym <- 2*(r-1)+1
    x_y <- expand.grid((x-xm):(x+xm),(y-ym):(y+ym))
    x_y_d <- cbind(x_y,d=sqrt( (0.5*sqrt(3)*(x_y[,1]-x))^2 + (0.5*(x_y[,2]-y))^2) )
    rownames(x_y_d) <- paste0(x_y[,1],"x",x_y[,2])
    x_y_d <- x_y_d[rownames(x_y_d)%in%rownames(d),]

    d[xy,rownames(x_y_d)] <- x_y_d[,3]
    d[rownames(x_y_d),xy] <- x_y_d[,3]
  }

  W <- exp(-d^2/2)

  if(diag0==TRUE) diag(W) <- 0

  W <- W[,colSums(W)>0] # remove spot island
  W <- W[rowSums(W)>0,] # remove spot island

  W
}
