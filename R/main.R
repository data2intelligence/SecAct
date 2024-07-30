#' @title Secreted protein activity inference
#' @description Infer the activity of over 1000 secreted proteins from tumor gene expression profiles.
#' @param expr Gene expression matrix with gene symbol (row) x sample (column).
#' @param sigMatrix Secreted protein signature matrix.
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
#' @examples
#'
#' exprPath <- file.path(system.file(package = "SecAct"), "extdata/IFNG_GSE100093.diff")
#' expr <- read.table(expr, sep="\t", check.names=F)
#' res <- SecAct.inference(expr, lambda=10000, nrand=1000)
#' head(res$zscore)
#'
#' @rdname SecAct.inference
#' @export
#'
SecAct.inference <- function(expr, SigMat="SecAct", lambda=5e+5, nrand=1000)
{
  if(SigMat=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  Y <- expr

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

#' @title Secreted protein signaling velocity
#' @description Calculate the signaling velocity of secreted proteins based on their activities.
#' @param SpaCET_obj An SpaCET object.
#' @param gene Gene symbol coding a secreted protein.
#' @param signalMode Mode of signaling velocity, i.e., "receiving", "sending", and "both".
#' @return A ggplot2 object.
#' @details
#' The velocity direction starts from the source cell producing a secreted protein and moves to sink cells receiving the secreted protein signal. The velocity magnitude represents the product between the secreted protein-coding gene expression at source cells and signaling activities at sink cells.
#'
#' @examples
#' SecAct.spatial.velocity(SpaCET_obj, gene="TGFB1", signalMode="receiving")
#' SecAct.spatial.velocity(SpaCET_obj, gene="TGFB1", signalMode="sending")
#' SecAct.spatial.velocity(SpaCET_obj, gene="TGFB1", signalMode="both")
#'
#' @rdname SecAct.spatial.velocity
#' @export
#'
SecAct.spatial.velocity <- function(SpaCET_obj, gene, signalMode="both")
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }

  counts <- SpaCET_obj@input$counts
  exp <- Matrix::t(Matrix::t(counts)*1e5/Matrix::colSums(counts))
  exp <- log2(exp+1)

  act <- SpaCET_obj@results$activity$zscore
  act[act<0] <- 0

  weights <- calWeights(colnames(exp), r=3, diag0=TRUE)

  act_new <- act[,colnames(weights)] # remove spot island
  exp_new <- exp[,colnames(weights)] # remove spot island

  if(!gene%in%rownames(exp_new))
  {
    weights_new <- weights * 0
  }else{
    weights_new <- weights * exp_new[gene,]
  }
  weights_new <- t(t(weights_new) * act_new[gene,])

  coordi <- data.frame(
    SpaCET_obj@input$spotCoordinates[colnames(weights_new),1],
    SpaCET_obj@input$spotCoordinates[colnames(weights_new),2])

  image <- SpaCET_obj@input$image
  xDiml <- dim(image$grob$raster)[1] # dim pixel
  yDiml <- dim(image$grob$raster)[2] # dim pixel

  coordi[,1] <- xDiml-coordi[,1]

  colnames(weights_new) <- paste0(coordi[,1],"x",coordi[,2])
  rownames(weights_new) <- paste0(coordi[,1],"x",coordi[,2])


  # sending
  startend <- data.frame()
  for(i in 1:nrow(weights_new))
  {
    vector_len <- sum(weights_new[i,])
    if(vector_len==0) next

    spot <- rownames(weights_new)[i]
    neighbors <- colnames(weights_new)[weights_new[i,]>0]
    neighbors_value <- weights_new[i,weights_new[i,]>0]

    spot_2col <- t(matrix(as.numeric(unlist(strsplit(spot,"x"))),nrow=2))
    neighbors_2col <- t(matrix(as.numeric(unlist(strsplit(neighbors,"x"))),nrow=2))
    neighbors_2col[,1] <- neighbors_2col[,1] - spot_2col[,1]
    neighbors_2col[,2] <- neighbors_2col[,2] - spot_2col[,2]

    neighbors_2col <- t(apply(neighbors_2col,1,scalar1))
    rownames(neighbors_2col) <- neighbors

    neighbors_2col <- neighbors_2col*neighbors_value
    neighbors_2col <- colMeans(neighbors_2col)
    neighbors_2col <- scalar1(neighbors_2col)
    neighbors_2col <- neighbors_2col * vector_len

    startend[spot,"x_start"] <- spot_2col[,1]
    startend[spot,"y_start"] <- spot_2col[,2]
    startend[spot,"x_change"] <- neighbors_2col[1]
    startend[spot,"y_change"] <- neighbors_2col[2]
    startend[spot,"x_end"] <- spot_2col[,1] + neighbors_2col[1]
    startend[spot,"y_end"] <- spot_2col[,2] + neighbors_2col[2]
    startend[spot,"vec_len"] <- sqrt(neighbors_2col[1]^2 + neighbors_2col[2]^2)
  }

  if(nrow(startend)!=0)
  {
    startend[,3] <- startend[,3]*10/max(abs(startend[,3]))
    startend[,4] <- startend[,4]*10/max(abs(startend[,4]))

    startend[,5] <- startend[,1] + startend[,3]
    startend[,6] <- startend[,2] + startend[,4]

    startend[startend[,"vec_len"]<0.1,"vec_len"] <- 0.01
    startend[startend[,"vec_len"]>=0.1,"vec_len"] <- 0.08
  }



  # receiving
  startend2 <- data.frame()
  for(i in 1:ncol(weights_new))
  {
    vector_len <- sum(weights_new[,i])
    if(vector_len==0) next

    spot <- colnames(weights_new)[i]
    neighbors <- rownames(weights_new)[weights_new[,i]>0]
    neighbors_value <- weights_new[weights_new[,i]>0,i]

    spot_2col <- t(matrix(as.numeric(unlist(strsplit(spot,"x"))),nrow=2))
    neighbors_2col <- t(matrix(as.numeric(unlist(strsplit(neighbors,"x"))),nrow=2))
    neighbors_2col[,1] <- -(neighbors_2col[,1] - spot_2col[,1])
    neighbors_2col[,2] <- -(neighbors_2col[,2] - spot_2col[,2])

    neighbors_2col <- t(apply(neighbors_2col,1,scalar1))
    rownames(neighbors_2col) <- neighbors

    neighbors_2col <- neighbors_2col*neighbors_value
    neighbors_2col <- colMeans(neighbors_2col)
    neighbors_2col <- scalar1(neighbors_2col)
    neighbors_2col <- neighbors_2col * vector_len

    startend2[spot,"x_start"] <- spot_2col[,1] - neighbors_2col[1]
    startend2[spot,"y_start"] <- spot_2col[,2] - neighbors_2col[2]
    startend2[spot,"x_change"] <- neighbors_2col[1]
    startend2[spot,"y_change"] <- neighbors_2col[2]
    startend2[spot,"x_end"] <- spot_2col[,1]
    startend2[spot,"y_end"] <- spot_2col[,2]
    startend2[spot,"vec_len"] <- sqrt(neighbors_2col[1]^2 + neighbors_2col[2]^2)
  }

  if(nrow(startend2)!=0)
  {
    startend2[,3] <- startend2[,3]*10/max(abs(startend2[,3]))
    startend2[,4] <- startend2[,4]*10/max(abs(startend2[,4]))

    startend2[,1] <- startend2[,5] - startend2[,3]
    startend2[,2] <- startend2[,6] - startend2[,4]

    startend2[startend2[,"vec_len"]<0.1,"vec_len"] <- 0.01
    startend2[startend2[,"vec_len"]>=0.1,"vec_len"] <- 0.08

    #startend2[,"vec_len"] <- startend2[,"vec_len"]/max(startend2[,"vec_len"])
    #startend2[,"vec_len"] <- startend2[,"vec_len"]/10
  }



  library(ggplot2)
  library(patchwork)

  fig.df <- data.frame(
    x=coordi[,1],
    y=coordi[,2],
    Exp=exp_new[gene,]
  )

  fig.df[fig.df[,3]>5,3] <- 5

  p1 <- ggplot(fig.df,aes(x=x,y=y))+
    annotation_custom(image$grob)+
    geom_point(aes(colour=Exp))+
    scale_color_gradient(low="blue",high="red")+
    scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0)) +
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), data=startend, arrow = arrow(length = unit(startend$vec_len, "cm")))+
    ggtitle(paste0(gene," (sending)"))+
    theme_void()+
    theme(
      plot.title = element_text(hjust = 0.5)
    )+coord_flip()


  fig.df2 <- data.frame(
    x=coordi[,1],
    y=coordi[,2],
    Act=act_new[gene,]
  )

  p2 <- ggplot(fig.df2,aes(x=x,y=y))+
    annotation_custom(image$grob)+
    geom_point(aes(colour=Act))+
    #scale_color_gradientn(colors=c("#a5a6ff","#ff72a1","brown"))+
    scale_color_gradientn(colors=c("#b8e186","#de77ae","#c51b7d"))+
    scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0)) +
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), data=startend2, arrow = arrow(length = unit(startend2$vec_len, "cm")))+
    ggtitle(paste0(gene," (receiving)"))+
    theme_void()+
    theme(
      plot.title = element_text(hjust = 0.5)
    )+coord_flip()
  p1+p2

  if(signalMode=="sending")
  {
    p1
  }else if(signalMode=="receiving"){
    p2
  }else{
    p1+p2
  }

}

scalar1 <- function(x)
{
  x / sqrt(sum(x^2))
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
