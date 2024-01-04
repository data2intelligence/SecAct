#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname hello
#' @export
hello <- function(n) {
  .C("hello", as.integer(n))
  #print(paste0("***",n,"***"))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' Yfile<- file.path(system.file(package = "SecAct"), "extdata/IFNG_GSE100093.diff")
#' Y <- read.table(Yfile,sep="\t",check.names=F)
#' res <- SecAct.inference(Y, lambda=10000, nrand=1000)
#' names(res)
#' head(res$zscore)
#' @rdname SecAct.inference
#' @export
SecAct.inference <- function(Y,lambda=10000,nrand=1000)
{
  Xfile<- file.path(system.file(package = "SecAct"), "extdata/signature.centroid")
  X <- read.table(Xfile,sep="\t",check.names=F)

  olp <- intersect(row.names(X),row.names(Y))
  X_olp <- as.matrix(X[olp,])
  Y_olp <- as.matrix(Y[olp,])

  n <- length(olp)
  p <- ncol(X_olp)
  m <- ncol(Y_olp)

  res2 <- .C("ridgeReg2", as.double(t(X_olp)), as.double(t(Y_olp)), as.integer(n), as.integer(p), as.integer(m), as.double(lambda), as.double(nrand), beta=double(p*m), se=double(p*m), zscore=double(p*m), pvalue=double(p*m))

  beta <- matrix(res2[[8]],byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  se <- matrix(res2[[9]],byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  zscore <- matrix(res2[[10]],byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  pvalue <- matrix(res2[[11]],byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))

  list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)
}


