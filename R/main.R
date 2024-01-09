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
#' library(SecAct)
#' Yfile<- file.path(system.file(package = "SecAct"), "extdata/IFNG_GSE100093.diff")
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
  if(is.null(SigMat))
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/signature.centroid")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(X),row.names(Y))
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
    X_olp=as.double(t(X)),
    Y_olp=as.double(t(Y)),
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

  list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)
}
