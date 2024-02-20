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

  if(class(Y)[1]=="SpaCET")
  {
    Y_type <- "SpaCET"
    SpaCET_obj <- Y

    Y <- SpaCET_obj@input$counts
    Y <- Matrix::t(Matrix::t(Y)*1e5/Matrix::colSums(Y))
    Y@x <- log2(Y@x + 1)
    Y <- Y - Matrix::rowMeans(Y)
  }else if (class(Y)[1]=="Seurat"){
    Y_type <- "Seurat"
    Seurat_obj <- Y

    Y <- Seurat_obj@assays$RNA@data
    Y <- Y - Matrix::rowMeans(Y)
  }else{
    Y_type <- "matrix"
  }

  if(is.null(SigMat))
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/signature.centroid")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

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
  }else if(Y_type=="Seurat"){
    res_z <- t(res$zscore)
    colnames(res_z) <- paste0("SecAct_",colnames(res_z))

    Seurat_obj@meta.data <- cbind(Seurat_obj@meta.data,res_z)
    Seurat_obj
  }else{
    res
  }

}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param SpaCET_obj PARAM_DESCRIPTION
#' @param gene PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#'
#' @rdname SecAct.signalling.direction
#' @export
#'
SecAct.signalling.direction <- function(SpaCET_obj, gene="TGFB1")
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }

  set.seed(123456)
  st.matrix.data <- SpaCET_obj@input$counts
  vst <- round(sctransform::vst(st.matrix.data, min_cells=5)$y,3)
  weights <- calWeights(colnames(vst), r=3, diag0=TRUE)
  act <- SpaCET_obj@results$SecAct_res$zscore

  act_new <- act[,colnames(weights)] # remove spot island
  vst_new <- vst[,colnames(weights)] # remove spot island

  act_new[act_new<0] <- 0
  vst_new[vst_new<0] <- 0

  weights_new <- weights * vst_new[gene,]
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

  startend[,3] <- startend[,3]*10/max(abs(startend[,3]))
  startend[,4] <- startend[,4]*10/max(abs(startend[,4]))

  startend[,5] <- startend[,1] + startend[,3]
  startend[,6] <- startend[,2] + startend[,4]

  startend[startend[,"vec_len"]<0.1,"vec_len"] <- 0.01
  startend[startend[,"vec_len"]>=0.1,"vec_len"] <- 0.08




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

  startend2[,3] <- startend2[,3]*10/max(abs(startend2[,3]))
  startend2[,4] <- startend2[,4]*10/max(abs(startend2[,4]))

  startend2[,1] <- startend2[,5] - startend2[,3]
  startend2[,2] <- startend2[,6] - startend2[,4]

  startend2[startend2[,"vec_len"]<0.1,"vec_len"] <- 0.01
  startend2[startend2[,"vec_len"]>=0.1,"vec_len"] <- 0.08


  library(ggplot2)
  library(patchwork)


  fig.df <- data.frame(
    x=coordi[,1],
    y=coordi[,2],
    value=vst_new[gene,]
  )

  fig.df[fig.df[,3]>5,3] <- 5

  p1 <- ggplot(fig.df,aes(x=x,y=y))+
    annotation_custom(image$grob)+
    geom_point(aes(colour=value),size=2.5)+
    scale_color_gradient(low="blue",high="green")+
    scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0)) +
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), data=startend, arrow = arrow(length = unit(startend$vec_len, "cm")))+
    ggtitle(paste0(gene," vst (sending)"))+
    theme_void()+
    theme(plot.title = element_text(hjust = 0.5))+coord_flip()


  fig.df2 <- data.frame(x=coordi[,1],y=coordi[,2],value=act_new[gene,])

  p2 <- ggplot(fig.df2,aes(x=x,y=y))+
    annotation_custom(image$grob)+
    geom_point(aes(colour=value),size=2.5)+
    scale_color_gradient(low="blue",high="red")+
    scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0)) +
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), data=startend2, arrow = arrow(length = unit(startend2$vec_len, "cm")))+
    ggtitle(paste0(gene," act (receiving)"))+
    theme_void()+
    theme(plot.title = element_text(hjust = 0.5))+coord_flip()

  p1|p2
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
