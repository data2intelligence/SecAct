#' @title Secreted protein signaling pattern
#' @description Calculate the signaling pattern of secreted proteins based on their activities.
#' @param SpaCET_obj A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param radius Radius cut off.
#' @param k Number of patterns for NMF.
#' @return A ggplot2 object.
#' @examples
#' SpaCET_obj <- SecAct.signaling.pattern(SpaCET_obj, k=3)
#'
#' @rdname SecAct.signaling.pattern
#' @export
#'
SecAct.signaling.pattern <- function(SpaCET_obj, scale.factor = 1e+05, radius=200, k)
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }
  if(is.null(SpaCET_obj @results $SecAct_output $SecretedProteinActivity))
  {
    stop("Please run SecAct.activity.inference first.")
  }

  act <- SpaCET_obj @results $SecAct_output $SecretedProteinActivity $zscore
  act[act<0] <- 0

  print("Step 1. Filtering")

  exp <- SpaCET_obj@input$counts
  rownames(exp) <- transferSymbol(rownames(exp))
  exp <- rm_duplicates(exp)

  # normalize to TPM
  stats <- Matrix::colSums(exp)
  exp <- sweep_sparse(exp,2,stats,"/")
  exp@x <- exp@x * scale.factor

  # transform to log space
  exp@x <- log2(exp@x + 1)

  ## only need SPs
  weights <- calWeights(
    SpaCET_obj@input$spotCoordinates[,c("coordinate_x_um","coordinate_y_um")],
    radius=radius, sigma=100, diagAsZero=TRUE)

  act_new <- act[,colnames(weights)] # remove spot island
  exp_new <- exp[,colnames(weights)] # remove spot island

  exp_new_aggr <- exp_new %*% weights

  corr <- data.frame()
  for(gene in rownames(act_new))
  {
    act_gene <- act_new[gene,]

    if(gene%in%rownames(exp_new))
    {
      exp_gene <- exp_new_aggr[gene,]

      cor_res <- cor.test(act_gene, exp_gene, method="spearman")

      corr[gene,"r"] <- cor_res$estimate
      corr[gene,"p"] <- cor_res$p.value
    }else{
      corr[gene,"r"] <- NA
      corr[gene,"p"] <- NA
    }
  }

  corr <- cbind(corr, padj=p.adjust(corr[,"p"], method="BH") )
  corr_genes <- rownames(corr[!is.na(corr[,"r"])&corr[,"r"]>0.05&corr[,"padj"]<0.01,])

  print(paste0(length(corr_genes),"/",nrow(act_new)," secreted proteins are kept to infer signaling patterns."))



  print("Step 2. NMF")

  suppressPackageStartupMessages({
    library(NMF)
  })
  act_nneg <- nneg(act[corr_genes,])

  if(length(k)==1)
  {
    NMF_res <- nmf(act_nneg, k, seed=123456)
  }else{
    estim.r <- nmf(act_nneg, k, nrun=10, seed=123456)
    v <- estim.r$measures$silhouette.coef

    v_diff <- v[1:(length(v)-1)]-v[2:length(v)]
    maxN <- which( v_diff == max(v_diff) ) +1
    k <- k[maxN]

    print(paste0("The optimal number of factors k = ", k))

    NMF_res <- nmf(act_nneg, k, seed=123456)
  }

  weight.W <- NMF_res@fit@W
  signal.H <- NMF_res@fit@H

  colnames(weight.W) <- as.character(1:k)
  rownames(signal.H) <- as.character(1:k)

  SpaCET_obj @results $SecAct_output $pattern <- list(
    ccc.SP = corr,
    weight.W = weight.W,
    signal.H = signal.H
  )

  SpaCET_obj
}


#' @title Pattern-associated secreted proteins
#' @description Enumerate secreted proteins associated with each signaling pattern.
#' @param SpaCET_obj A SpaCET object.
#' @param n Pattern order.
#' @return A matrix.
#' @examples
#' SpaCET_obj <- SecAct.pattern.gene(SpaCET_obj, n=3)
#'
#' @rdname SecAct.signaling.pattern.gene
#' @export
#'
SecAct.signaling.pattern.gene <- function(SpaCET_obj, n)
{
  # extract coefficient (weight) matrix
  weight.W <- SpaCET_obj@results$SecAct_output$pattern$weight.W

  temp <- weight.W
  temp[,-n] <- 2*temp[,-n] # in case one column

  # identify secreted proteins with pattern n
  res <- weight.W[apply(temp,1,function(x) x[n]==max(x)),]
  res[order(res[,3],decreasing = TRUE),]
}


#' @title Secreted protein signaling velocity
#' @description Calculate the signaling velocity of secreted proteins based on their activities.
#' @param SpaCET_obj A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param gene Gene symbol coding a secreted protein.
#' @param signalMode Mode of signaling velocity, i.e., "receiving", "sending", and "both".
#' @param radius Radius cut off.
#' @param contourMap A logical flag indicating whether to transform as contour map.
#' @param animated A logical flag indicating whether to generate animated figure.
#' @return A ggplot2 object.
#' @details
#' The velocity direction starts from the source cell producing a secreted protein and moves to sink cells receiving the secreted protein signal. The velocity magnitude represents the product between the secreted protein-coding gene expression at source cells and signaling activities at sink cells.
#'
#' @examples
#' SecAct.signaling.velocity.spotST(SpaCET_obj, gene="TGFB1", signalMode="receiving")
#' SecAct.signaling.velocity.spotST(SpaCET_obj, gene="TGFB1", signalMode="sending")
#'
#' @rdname SecAct.signaling.velocity.spotST
#' @export
#'
SecAct.signaling.velocity.spotST <- function(
  SpaCET_obj,
  scale.factor = 1e+05,
  gene,
  signalMode="receiving",
  radius=200,
  contourMap=FALSE,
  coutourBins=11,
  animated=FALSE
)
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }
  if(is.null(SpaCET_obj @results $SecAct_output $SecretedProteinActivity))
  {
    stop("Please run SecAct.activity.inference first.")
  }
  if(contourMap==TRUE&animated==TRUE)
  {
    stop("contourMap and animated can not be TRUE simultaneously.")
  }

  library(ggplot2)
  library(patchwork)

  act <- SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore
  act[act<0] <- 0

  exp <- SpaCET_obj@input$counts
  rownames(exp) <- transferSymbol(rownames(exp))
  exp <- rm_duplicates(exp)

  # normalize to TPM
  stats <- Matrix::colSums(exp)
  exp <- sweep_sparse(exp,2,stats,"/")
  exp@x <- exp@x * scale.factor

  # transform to log space
  exp@x <- log2(exp@x + 1)

  weights <- calWeights(
    SpaCET_obj@input$spotCoordinates[,c("coordinate_x_um","coordinate_y_um")],
    radius=radius, sigma=100, diagAsZero=TRUE)

  act_new <- act[,colnames(weights)] # remove spot island
  exp_new <- exp[,colnames(weights)] # remove spot island

  if(!gene%in%rownames(exp_new))
  {
    weights_new <- weights * 0
  }else{
    weights_new <- weights * exp_new[gene,]
  }

  weights_new <- Matrix::t(Matrix::t(weights_new) * act_new[gene,])

  coordi <- data.frame(
    SpaCET_obj@input$spotCoordinates[colnames(weights_new),1],
    SpaCET_obj@input$spotCoordinates[colnames(weights_new),2])

  image <- SpaCET_obj@input$image
  xDiml <- dim(image$grob$raster)[1] # dim pixel
  yDiml <- dim(image$grob$raster)[2] # dim pixel

  coordi[,1] <- xDiml-coordi[,1]

  colnames(weights_new) <- paste0(coordi[,1],"x",coordi[,2])
  rownames(weights_new) <- paste0(coordi[,1],"x",coordi[,2])


  if(signalMode=="sending")
  {
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

  }else{ # receiving

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

      startend2[startend2[,"vec_len"]<0.1,"vec_len"] <- 0.001
      startend2[startend2[,"vec_len"]>=0.1,"vec_len"] <- 0.08
    }

    startend <- startend2
  }


  if(signalMode=="sending")
  {
    fig.df <- data.frame(
      x=coordi[,1],
      y=coordi[,2],
      value=exp_new[gene,]
    )
    fig.df[fig.df[,3]>5,3] <- 5

  }else{  # receiving
    fig.df <- data.frame(
      x=coordi[,1],
      y=coordi[,2],
      value=act_new[gene,]
    )

  }


  if(animated==TRUE)
  {
    startend_temp <- startend
    startend_temp[,"x_end"] <- startend_temp[,"x_start"]
    startend_temp[,"y_end"] <- startend_temp[,"y_start"]

    startend <- rbind(
      cbind(startend_temp,tim=1),
      cbind(startend,tim=2)
    )
  }

  if(contourMap==TRUE)
  {
    library(akima)
    # Interpolate onto grid
    interp_res <- interp(
      x = fig.df$x,
      y = fig.df$y,
      z = fig.df$value,
      xo = seq(min(fig.df$x), max(fig.df$x), length = 20),  # grid resolution in x
      yo = seq(min(fig.df$y), max(fig.df$y), length = 20)   # grid resolution in y
    )

    # Convert to data frame for ggplot
    grid_df <- data.frame(
      x = rep(interp_res$x, times = length(interp_res$y)),
      y = rep(interp_res$y, each = length(interp_res$x)),
      z = as.vector(interp_res$z)
    )

    # remove arrow out of contour
    grid_df_noNA <- grid_df[!is.na(grid_df$z),]
    startend <- startend[
      startend[,"x_start"]>min(grid_df_noNA$x) & startend[,"x_start"]<max(grid_df_noNA$x) &
      startend[,"x_end"]>min(grid_df_noNA$x) & startend[,"x_end"]<max(grid_df_noNA$x) &
      startend[,"y_start"]>min(grid_df_noNA$y) & startend[,"y_start"]<max(grid_df_noNA$y) &
      startend[,"y_end"]>min(grid_df_noNA$y) & startend[,"y_end"]<max(grid_df_noNA$y)
      ,]

    p <- ggplot(grid_df,aes(x=x,y=y))+
      geom_contour_filled(aes(z=z),bins = coutourBins) +
      scale_fill_brewer(palette = "RdYlGn",direction = -1)+
      #scale_fill_manual(values=c("#b8e186","#de77ae","#c51b7d"))+
      #scale_fill_manual(values=c("#000004FF","#1C1044FF","#4F127BFF","#812581FF","#B5367AFF","#E55064FF","#FB8761FF","#FEC287FF","#FCFDBFFF"))+
      scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0)) +
      ggtitle(paste0(gene," (",signalMode,")"))+
      theme_void()+
      theme(
        plot.title = element_text(hjust = 0.5)
      )+coord_flip()
    p <- p+geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), color="black", data=startend, arrow = arrow(length = unit(startend$vec_len, "cm")),inherit.aes = FALSE)

  }else{
    p <- ggplot(fig.df,aes(x=x,y=y))+
      #annotation_custom(image$grob)+
      geom_point(aes(colour=value))+
      #scale_color_brewer(palette = "RdYlGn",direction = -1)+
      #scale_color_gradientn(colors=c("#a5a6ff","#ff72a1","brown"))+
      scale_color_gradientn(colors=c("#b8e186","#de77ae","#c51b7d"))+
      scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0)) +
      geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), color="black", data=startend, arrow = arrow(length = unit(startend$vec_len, "cm")),inherit.aes = FALSE)+
      ggtitle(paste0(gene," (",signalMode,")"))+
      theme_void()+
      theme(
        plot.title = element_text(hjust = 0.5)
      )+coord_flip()
  }

  if(animated==TRUE)
  {
    library(gganimate)
    p <- animate(p+transition_time(tim),nframes=15)
  }

  p
}


#' @title Secreted protein signaling velocity
#' @description Calculate the signaling velocity of secreted proteins based on their activities.
#' @param SpaCET_obj A SpaCET object.
#' @param sender Sender cell types.
#' @param secretedProtein Secreted proteins.
#' @param receiver Receiver cell types.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param CustomizedAreaCoordinates A vector of four numbers for coordinates of the Customized Area, i.e., x_left, x_right, y_bottom, y_top.
#' @param radius Radius cut off (unit: um).
#' @param show.coordinates Whether to show coordinates.
#' @param colors A vector of colors for spots.
#' @param pointSize Size of spots.
#' @param pointAlpha Alpha transparency scales of spots. Must lie between 0 and 1.
#' @param legend.position The position of the legend. Set it as "none" if you want to remove the legend.
#' @param legend.size The size of the legend.
#' @param arrow.color The color of the arrow.
#' @param arrow.width The width of the arrow.
#' @param arrow.size The size of the arrow.
#' @return A ggplot2 object.
#' @details
#' The velocity direction starts from the source cell producing a secreted protein and moves to sink cells receiving the secreted protein signal. The velocity magnitude represents the product between the secreted protein-coding gene expression at source cells and signaling activities at sink cells.
#'
#' @examples
#' SecAct.signaling.velocity.scST(SpaCET_obj, sender="Fibroblast", secretedProtein="THBS2", receiver="Tumor_boundary", cellType_meta="cellType")
#' SecAct.signaling.velocity.scST(SpaCET_obj, sender="Fibroblast", secretedProtein="THBS2", receiver="Tumor_boundary", cellType_meta="cellType", CustomizedAreaCoordinates=c(8290,8366,1100,1400))
#'
#' @rdname SecAct.signaling.velocity.scST
#' @export
#'
SecAct.signaling.velocity.scST <- function(
    SpaCET_obj,
    sender,
    secretedProtein,
    receiver,
    cellType_meta,
    scale.factor = 1e+05,
    CustomizedAreaCoordinates = NULL,
    radius = 20,
    show.coordinates = TRUE,
    colors = NULL,
    pointSize = 1,
    pointAlpha = 1,
    legend.position = "right",
    legend.size = 1,
    arrow.color = "#ff0099",
    arrow.width = 1,
    arrow.size = 0.3
)
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }
  if(is.null(SpaCET_obj @results $SecAct_output $SecretedProteinActivity))
  {
    stop("Please run SecAct.activity.inference first.")
  }

  coordinate_mat <- SpaCET_obj@input$spotCoordinates
  cellType_vec <- SpaCET_obj@input$metaData[,cellType_meta]


  cellType1_cells <- which(cellType_vec==sender)
  cellType2_cells <- which(cellType_vec==receiver)


  nn_result <- RANN::nn2(coordinate_mat[,c("coordinate_x_um","coordinate_y_um")], k=100, searchtype="radius", radius=radius)

  neighbor_indices <- nn_result$nn.idx
  neighbor_distances <- nn_result$nn.dists

  i <- rep(1:nrow(neighbor_indices), each=ncol(neighbor_indices)) # row indices (cell index)
  j <- as.vector(t(neighbor_indices))
  x <- as.vector(t(neighbor_distances))

  valid <- x<=radius & x>0
  i <- i[valid]          # Keep only valid indices
  j <- j[valid]          # Valid neighbor indices
  x <- x[valid]          # Valid distances


  exp <- SpaCET_obj@input$counts
  rownames(exp) <- transferSymbol(rownames(exp))
  exp <- rm_duplicates(exp)


  act <- SpaCET_obj @results $SecAct_output $SecretedProteinActivity$zscore
  act[act<0] <- 0



  Tmat <- data.frame(i,j)

  # all cell pair
  Tmat_cellTypePair <- Tmat[Tmat[,"i"]%in%cellType1_cells & Tmat[,"j"]%in%cellType2_cells, ,drop=F]


  # expr(1) * act(2)
  CCC_vec <- exp[secretedProtein, Tmat_cellTypePair[,1]] * act[secretedProtein, Tmat_cellTypePair[,2]]

  Tmat_cellTypePair <- Tmat_cellTypePair[CCC_vec>0,]


  fg.df <- data.frame(coordinate_mat, cellType=cellType_vec)

  fg.df[!fg.df[,3]%in%c(sender,receiver),3] <- "Other"


  startend <- data.frame(
    sender=rownames(coordinate_mat)[Tmat_cellTypePair[,1]],
    receiver=rownames(coordinate_mat)[Tmat_cellTypePair[,2]],
    vec_len=1
  )

  startend <- cbind(startend, x_start=coordinate_mat[startend[,1],1] )
  startend <- cbind(startend, y_start=coordinate_mat[startend[,1],2] )
  startend <- cbind(startend, x_end=coordinate_mat[startend[,2],1] )
  startend <- cbind(startend, y_end=coordinate_mat[startend[,2],2] )

  if(is.null(CustomizedAreaCoordinates))
  {
    p <- ggplot(fg.df, aes(coordinate_x_um, coordinate_y_um)) + #sdimx, sdimy
      geom_point(aes(colour=cellType), size=pointSize, alpha=pointAlpha) +
      geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), data=startend, color=arrow.color,
                   arrow = arrow(length = unit(arrow.size, "cm")), linewidth=arrow.width)+
      scale_color_manual(values=colors)+
      ggtitle(" ")+
      xlab(" ")+
      ylab(" ")

    if(show.coordinates==TRUE)
    {
      p <- p+theme_classic()
    }else{
      p <- p+theme_void()
    }

    p+theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = legend.position
    )+
    guides(color = guide_legend(override.aes = list(size = legend.size)))

  }else{
    x.left <- CustomizedAreaCoordinates[1]
    x.right <- CustomizedAreaCoordinates[2]
    y.bottom <- CustomizedAreaCoordinates[3]
    y.top <- CustomizedAreaCoordinates[4]

    fg.df_cut <- fg.df[
      fg.df[,1]> x.left&
        fg.df[,1]< x.right&
        fg.df[,2]> y.bottom&
        fg.df[,2]< y.top,
    ]

    startend_cut <- startend[startend[,1]%in%rownames(fg.df_cut) & startend[,2]%in%rownames(fg.df_cut),]

    p <- ggplot(fg.df_cut, aes(coordinate_x_um, coordinate_y_um)) +
      geom_point(aes(color=cellType), size=pointSize, alpha=pointAlpha) +
      geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), data=startend_cut, color=arrow.color,
                   arrow = arrow(length = unit(arrow.size, "cm")), linewidth=arrow.width)+
      scale_color_manual(values=colors)+
      ggtitle(" ")+
      xlab(" ")+
      ylab(" ")

    if(show.coordinates==TRUE)
    {
      p <- p+theme_classic()
    }else{
      p <- p+theme_void()
    }

    p+theme(
        plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = legend.position
      )+
      guides(color = guide_legend(override.aes = list(size = legend.size)))
  }
}


#' @title Cell-cell communication from spatial data
#' @description Calculate cell-cell communication mediated by secreted proteins from spatial transcriptomics data.
#' @param SpaCET_obj A SpaCET object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param radius Radius cut off (unit: um).
#' @param ratio_cutoff Ratio cut off.
#' @param padj_cutoff Adjusted p value cut off.
#' @param coreNo Core number in parallel computation.
#' @return A Seurat object.
#' @rdname SecAct.CCC.scST
#' @export
#'
SecAct.CCC.scST <- function(
    SpaCET_obj,
    cellType_meta,
    scale.factor = 1000,
    radius = 20,
    ratio_cutoff = 0.2,
    padj_cutoff = 0.01,
    coreNo=6
)
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }
  if(is.null(SpaCET_obj @results $SecAct_output $SecretedProteinActivity))
  {
    stop("Please run SecAct.activity.inference first.")
  }

  coreNoDect <- parallel::detectCores(logical = FALSE)
  if(coreNoDect<coreNo)
  {
    message(paste0("Since the number of your physical cores is ",coreNoDect,", coreNo=",coreNoDect," is used automatically."))
    coreNo <- coreNoDect
  }

  coordinate_mat <- SpaCET_obj@input$spotCoordinates
  cellType_vec <- SpaCET_obj@input$metaData[,cellType_meta]

  print("Step 1. Filtering")

  exp <- SpaCET_obj@input$counts
  rownames(exp) <- transferSymbol(rownames(exp))
  exp <- rm_duplicates(exp)

  # normalize to TPM
  stats <- Matrix::colSums(exp)
  exp <- sweep_sparse(exp,2,stats,"/")
  exp@x <- exp@x * scale.factor

  # transform to log space
  exp@x <- log2(exp@x + 1)


  nn_result <- RANN::nn2(coordinate_mat[,c("coordinate_x_um","coordinate_y_um")], k=100, searchtype="radius", radius=radius)

  neighbor_indices <- nn_result$nn.idx
  neighbor_distances <- nn_result$nn.dists

  i <- rep(1:nrow(neighbor_indices), each=ncol(neighbor_indices)) # row indices (cell index)
  j <- as.vector(t(neighbor_indices))
  x <- as.vector(t(neighbor_distances))

  valid <- x<=radius & x>0
  i <- i[valid]          # Keep only valid indices
  j <- j[valid]          # Valid neighbor indices
  x <- x[valid]          # Valid distances

  # Create the sparse matrix using the 'i', 'j', and 'x' vectors
  library(Matrix)
  distance_mat <- sparseMatrix(i=i, j=j, x=x, dims=c(nrow(neighbor_indices), nrow(neighbor_indices)), repr="T")
  rownames(distance_mat) <- rownames(coordinate_mat)
  colnames(distance_mat) <- rownames(coordinate_mat)


  weights <- distance_mat
  weights@x <- as.numeric(weights@x>0)


  act <- SpaCET_obj @results $SecAct_output $SecretedProteinActivity$zscore
  act[act<0] <- 0

  act_new <- act[,colnames(weights)] # remove spot island
  exp_new <- exp[,colnames(weights)] # remove spot island


  exp_new_aggr <- exp_new %*% weights

  if(is.null(SpaCET_obj @results $SecAct_output $ccc.SP))
  {
    corr <- data.frame()
    for(gene in rownames(act_new))
    {
      act_gene <- act_new[gene,]

      if(gene%in%rownames(exp_new))
      {
        exp_gene <- exp_new_aggr[gene,]

        cor_res <- cor.test(act_gene, exp_gene, method="spearman")

        corr[gene,"r"] <- cor_res$estimate
        corr[gene,"p"] <- cor_res$p.value
      }else{
        corr[gene,"r"] <- NA
        corr[gene,"p"] <- NA
      }
    }
    corr <- cbind(corr, padj=p.adjust(corr[,"p"], method="BH") )
  }else{
    SpaCET_obj @results $SecAct_output $ccc.SP -> corr
  }

  corr_genes <- rownames(corr[!is.na(corr[,"r"])&corr[,"r"]>0.05&corr[,"padj"]<0.01,])

  print(paste0(length(corr_genes),"/",nrow(act_new)," secreted proteins are kept to infer cell-cell communication."))


  print("Step 2. CCC")

  cellTypes <- unique(cellType_vec)
  cellGroups <- split(seq_along(cellType_vec), cellType_vec)

  cellTypePair1 <- rep(cellTypes, each=length(cellTypes))
  cellTypePair2 <- rep(cellTypes, length(cellTypes))

  cellTypePair <- cbind(cellTypePair1, cellTypePair2)
  cellTypePair <- cellTypePair[cellTypePair1>cellTypePair2,]
  rownames(cellTypePair) <- paste0(cellTypePair[,1],"_",cellTypePair[,2])

  olp <- corr_genes

  Tmat <- data.frame(i,j)

  compute_pair <- function(m)
  {
    cellType1 <- cellTypePair[m,1]
    cellType2 <- cellTypePair[m,2]

    cellType1_cells <- cellGroups[[cellType1]]
    cellType2_cells <- cellGroups[[cellType2]]

    # all neighboring cell pair
    Tmat_cellTypePair <- Tmat[Tmat[,"i"]%in%cellType1_cells & Tmat[,"j"]%in%cellType2_cells, ,drop=F]
    n_neighbor <- nrow(Tmat_cellTypePair)

    # ratio neighboring cells vs all cells
    if(
      length(unique(Tmat_cellTypePair[,1]))/length(cellType1_cells) < 0.05 &
      length(unique(Tmat_cellTypePair[,2]))/length(cellType2_cells) < 0.05)
    return(NULL)

    set.seed(123)
    Tmat_background <- data.frame(
      i=sample(cellType1_cells, n_neighbor*1000, replace=T),
      j=sample(cellType2_cells, n_neighbor*1000, replace=T)
    )

    out <- list()

    for(SP in olp)
    {
      exp_SP <- exp[SP,]
      act_SP <- act[SP,]

      # exp(1) * act(2)
      CCC_vec <- exp_SP[Tmat_cellTypePair[,1]] * act_SP[Tmat_cellTypePair[,2]]

      # ratio interacting pairs vs neighboring pairs
      n_communication <- sum(CCC_vec>0)
      posRatio <- n_communication/n_neighbor

      if(posRatio > ratio_cutoff)
      {
        CCC_raw <- mean(CCC_vec)

        CCC1000_vec <- exp_SP[Tmat_background[,1]] * act_SP[Tmat_background[,2]]
        CCC1000 <- colMeans(matrix(CCC1000_vec, nrow = n_neighbor))

        out[[paste0(cellType1,"_",SP,"_",cellType2)]] <- c(
          sender = cellType1,
          secretedProtein = SP,
          receiver = cellType2,
          sender_count = length(cellType1_cells),
          receiver_count = length(cellType2_cells),
          neighboringCellPairs = n_neighbor,
          communicatingCellPairs = n_communication,
          ratio = posRatio,
          #score = CCC_raw/mean(CCC1000),
          pv = (sum(CCC1000 >= CCC_raw) + 1) / 1001
        )
      }


      # exp(2) * act(1)
      CCC_vec <- exp_SP[Tmat_cellTypePair[,2]] * act_SP[Tmat_cellTypePair[,1]]

      n_communication <- sum(CCC_vec>0)
      posRatio <- n_communication/n_neighbor

      if(posRatio > ratio_cutoff)
      {
        CCC_raw <- mean(CCC_vec)

        CCC1000_vec <- exp_SP[Tmat_background[,2]] * act_SP[Tmat_background[,1]]
        CCC1000 <- colMeans(matrix(CCC1000_vec, nrow = n_neighbor))

        out[[paste0(cellType2,"_",SP,"_",cellType1)]] <- c(
          sender = cellType2,
          secretedProtein = SP,
          receiver = cellType1,
          sender_count = length(cellType2_cells),
          receiver_count = length(cellType1_cells),
          neighboringCellPairs = n_neighbor,
          communicatingCellPairs = n_communication,
          ratio = posRatio,
          #score = CCC_raw/mean(CCC1000),
          pv = (sum(CCC1000 >= CCC_raw) + 1) / 1001
        )
      }

    } # SP

    return(out)
  } #m

  results <- parallel::mclapply(1:nrow(cellTypePair), compute_pair, mc.cores = coreNo)
  ccc_list <- do.call(c, results)
  ccc <- do.call(rbind, ccc_list)
  ccc <- as.data.frame(ccc)
  ccc[,4:9] <- apply(ccc[,4:9], 2, as.numeric)

  ccc <- cbind(ccc, pv.adj=p.adjust(ccc[,"pv"], method="BH"))
  ccc <- ccc[ccc[,"pv.adj"]<padj_cutoff,]

  ccc <- ccc[order(ccc[,"pv.adj"]),]

  SpaCET_obj @results $SecAct_output $ccc.SP <- corr
  SpaCET_obj @results $SecAct_output $SecretedProteinCCC  <- ccc

  SpaCET_obj
}


#' @title Cell-cell communication from single cell data
#' @description Calculate condition-specific cell-cell communication mediated by secreted proteins from scRNA-Seq data.
#' @param Seurat_obj A Seurat object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param condition_meta Column name in meta data that includes condition information.
#' @param conditionCase Case condition.
#' @param conditionControl Control condition.
#' @param scale.factor Sets the scale factor for cell-level normalization in step2.
#' @param act_diff_cutoff Cut off for activity change (i.e., z score) in step 1.
#' @param exp_logFC_cutoff Cut off for log fold change in step 2.
#' @param exp_fraction_case_cutoff Cut off for the fraction of cells expressing secreted protein-coding genes in step 2.
#' @param padj_cutoff Adjusted p value cut off.
#' @param sigMatrix Secreted protein signature matrix.
#' @param is.group.sig A logical flag indicating whether to group similar signatures.
#' @param is.group.cor Correlation cutoff of similar signatures.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @return A Seurat object.
#' @rdname SecAct.CCC.scRNAseq
#' @export
#'
SecAct.CCC.scRNAseq <- function(
  Seurat_obj,
  cellType_meta,
  condition_meta,
  conditionCase,
  conditionControl,
  scale.factor = 1e+05,
  act_diff_cutoff = 2,
  exp_logFC_cutoff = 0.2,
  exp_mean_all_cutoff = 2,
  exp_fraction_case_cutoff = 0.1,
  padj_cutoff = 0.01,
  sigMatrix="SecAct",
  is.group.sig=TRUE,
  is.group.cor=0.9,
  lambda=5e+05,
  nrand=1000
)
{
  if(!class(Seurat_obj)[1]=="Seurat")
  {
    stop("Please input a Seurat object.")
  }

  if(class(Seurat_obj@assays$RNA)=="Assay5")
  {
    counts <- Seurat_obj@assays$RNA@layers$counts
    colnames(counts) <- rownames(Seurat_obj@assays$RNA@cells)
    rownames(counts) <- rownames(Seurat_obj@assays$RNA@features)
  }else{
    counts <-  Seurat_obj@assays$RNA@counts
  }

  rownames(counts) <- transferSymbol(rownames(counts))
  counts <- rm_duplicates(counts)

  meta <- Seurat_obj@meta.data

  if(is.null(condition_meta))
  {
    cellTypes <- unique(meta[,cellType_meta])
  }else{
    cellTypes <- intersect(
      meta[meta[,condition_meta]==conditionCase,cellType_meta],
      meta[meta[,condition_meta]==conditionControl,cellType_meta]
    )
  }


  print("Step 1: assessing changes in secreted protein expression.")

  if(sigMatrix=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    SPs <- strsplit(readLines(Xfile, n = 1), "\t")[[1]]
  }else{
    SPs <- strsplit(readLines(sigMatrix, n = 1), "\t")[[1]]
  }

  for(cellType in cellTypes)
  {
    # case
    if(is.null(condition_meta))
    {
      expr <- counts[,meta[,cellType_meta]==cellType,drop=FALSE]
    }else{
      expr <- counts[,meta[,condition_meta]==conditionCase&meta[,cellType_meta]==cellType,drop=FALSE]
    }
    if(ncol(expr)<30) next

    # normalize to TPM
    stats <- Matrix::colSums(expr)
    expr <- sweep_sparse(expr,2,stats,"/")
    expr@x <- expr@x * scale.factor

    # transform to log space
    expr@x <- log2(expr@x + 1)

    expr_case <- expr


    # control
    if(is.null(condition_meta))
    {
      expr <- counts[,meta[,cellType_meta]!=cellType,drop=FALSE]
    }else{
      expr <- counts[,meta[,condition_meta]==conditionControl&meta[,cellType_meta]==cellType,drop=FALSE]
    }

    # normalize to TPM
    stats <- Matrix::colSums(expr)
    expr <- sweep_sparse(expr,2,stats,"/")
    expr@x <- expr@x * scale.factor

    # transform to log space
    expr@x <- log2(expr@x + 1)

    expr_control <- expr


    # case vs control
    genes <- intersect(rownames(expr_case), SPs)

    case_mat <- expr_case[genes, , drop = FALSE]
    control_mat <- expr_control[genes, , drop = FALSE]

    # Compute summary statistics
    exp_mean_case     <- Matrix::rowMeans(case_mat)
    exp_mean_control  <- Matrix::rowMeans(control_mat)
    exp_mean_all      <- (exp_mean_case * ncol(case_mat) + exp_mean_control * ncol(control_mat)) / (ncol(case_mat) + ncol(control_mat))
    exp_logFC         <- exp_mean_case - exp_mean_control

    exp_fraction_case    <- Matrix::rowMeans(case_mat > 0)
    exp_fraction_control <- Matrix::rowMeans(control_mat > 0)

    # Wilcoxon tests row-wise
    pvals <- mapply(
      function(tc, cc) wilcox.test(tc, cc)$p.value,
      split(case_mat, row(case_mat)),
      split(control_mat, row(control_mat))
    )

    # Assemble into one data frame
    smy_deg <- data.frame(
      exp_logFC = exp_logFC,
      exp_mean_all = exp_mean_all,
      exp_mean_case = exp_mean_case,
      exp_mean_control = exp_mean_control,
      exp_fraction_case = exp_fraction_case,
      exp_fraction_control = exp_fraction_control,
      exp_pv = pvals
    )

    # Assign rownames
    rownames(smy_deg) <- genes

    smy_deg <- smy_deg[smy_deg[,"exp_mean_all"]>0,]

    smy_deg[,"exp_pv.adj"] <- p.adjust(smy_deg[,"exp_pv"], method="BH")
    smy_deg <- smy_deg[order(smy_deg[,"exp_pv.adj"]),]

    Seurat_obj @misc $SecAct_output $SecretedProteinExpression [[cellType]] <- smy_deg
  }


  print("Step 2: calculating changes in secreted protein activity.")

  bulk.diff <- data.frame()
  for(cellType in cellTypes)
  {
    # case
    if(is.null(condition_meta))
    {
      expr <- counts[,meta[,cellType_meta]==cellType,drop=FALSE]
    }else{
      expr <- counts[,meta[,condition_meta]==conditionCase&meta[,cellType_meta]==cellType,drop=FALSE]
    }
    if(ncol(expr)<30) next

    expr <- Matrix::rowSums(expr)
    expr <- matrix(expr,ncol=1,dimnames = list(names(expr),"bulk"))

    # normalize to TPM
    expr <- t(t(expr)*1e6/colSums(expr))

    # transform to log space
    expr <- log2(expr + 1)

    expr_case <- expr


    # control
    if(is.null(condition_meta))
    {
      expr <- counts[,meta[,cellType_meta]!=cellType,drop=FALSE]
    }else{
      expr <- counts[,meta[,condition_meta]==conditionControl&meta[,cellType_meta]==cellType,drop=FALSE]
    }

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

  Seurat_obj @misc $SecAct_output $SecretedProteinActivity <-
    SecAct.activity.inference(
      inputProfile = bulk.diff,
      is.differential = TRUE,
      sigMatrix = sigMatrix,
      is.group.sig = is.group.sig,
      is.group.cor = is.group.cor,
      lambda = lambda,
      nrand = nrand
      )


  print("Step 3: linking source and receiver cell types.")

  smy_deg_comb <- data.frame()
  smy_act_comb <- data.frame()

  cellTypes <- names(Seurat_obj @misc $SecAct_output $SecretedProteinExpression)

  if(length(cellTypes)<2) stop("Only one cell type with more than 30 cells.")

  for(cellType in cellTypes)
  {
    smy_deg <- Seurat_obj @misc $SecAct_output $SecretedProteinExpression [[cellType]]

    smy_deg_up <- smy_deg[
      smy_deg[,"exp_logFC"]>exp_logFC_cutoff&
      smy_deg[,"exp_mean_all"]>exp_mean_all_cutoff&
      smy_deg[,"exp_fraction_case"]>exp_fraction_case_cutoff&
      smy_deg[,"exp_pv.adj"]<padj_cutoff, ]

    if(nrow(smy_deg_up)>0)
    {
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"cellType"] <- cellType
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"secretedProtein"] <- rownames(smy_deg_up)
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"exp_logFC"] <- smy_deg_up[,"exp_logFC"]
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"exp_pv"] <- smy_deg_up[,"exp_pv"]
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"exp_pv.adj"] <- smy_deg_up[,"exp_pv.adj"]
    }

    res <- Seurat_obj @misc $SecAct_output $SecretedProteinActivity

    smy_act <- data.frame()
    for(gene in rownames(res$zscore))
    {
      smy_act[gene,"act_diff"] <- res$zscore[gene,cellType]
      smy_act[gene,"act_pv"] <- res$pvalue[gene,cellType]
      smy_act[gene,"act_pv.adj"] <- res$pvalue[gene,cellType]
    }
    smy_act[,"act_pv.adj"] <- p.adjust(smy_act[,"act_pv"], method="BH")
    smy_act <- smy_act[order(smy_act[,"act_pv.adj"]),]

    smy_act_up <- smy_act[
      smy_act[,"act_diff"]>act_diff_cutoff&
      smy_act[,"act_pv.adj"]<padj_cutoff, ]

    if(nrow(smy_act_up)>0)
    {
      smy_act_comb[paste0(cellType,"_",rownames(smy_act_up)),"cellType"] <- cellType
      smy_act_comb[paste0(cellType,"_",rownames(smy_act_up)),"secretedProtein"] <- rownames(smy_act_up)
      smy_act_comb[paste0(cellType,"_",rownames(smy_act_up)),"act_diff"] <- smy_act_up[,"act_diff"]
      smy_act_comb[paste0(cellType,"_",rownames(smy_act_up)),"act_pv"] <- smy_act_up[,"act_pv"]
      smy_act_comb[paste0(cellType,"_",rownames(smy_act_up)),"act_pv.adj"] <- smy_act_up[,"act_pv.adj"]
    }
  }


  # calculate cell-cell communication
  ccc <- data.frame()

  for(i in 1:nrow(smy_deg_comb))
  {
    sender <- smy_deg_comb[i,"cellType"]
    SP <- smy_deg_comb[i,"secretedProtein"]

    smy_act_comb_sub <- smy_act_comb[smy_act_comb[,"secretedProtein"]==SP,]
    if(nrow(smy_act_comb_sub)==0) next
    receiver <- smy_act_comb_sub[,"cellType"]

    ccc[paste0(sender,"_",SP,"_",receiver),"sender"] <- sender
    ccc[paste0(sender,"_",SP,"_",receiver),"secretedProtein"] <- SP
    ccc[paste0(sender,"_",SP,"_",receiver),"receiver"] <- receiver
    ccc[paste0(sender,"_",SP,"_",receiver),"sender_exp_logFC"] <- smy_deg_comb[i,"exp_logFC"]
    ccc[paste0(sender,"_",SP,"_",receiver),"sender_exp_pv"] <- smy_deg_comb[i,"exp_pv"]
    ccc[paste0(sender,"_",SP,"_",receiver),"sender_exp_pv.adj"] <- smy_deg_comb[i,"exp_pv.adj"]
    ccc[paste0(sender,"_",SP,"_",receiver),"receiver_act_diff"] <- smy_act_comb_sub[,"act_diff"]
    ccc[paste0(sender,"_",SP,"_",receiver),"receiver_act_pv"] <- smy_act_comb_sub[,"act_pv"]
    ccc[paste0(sender,"_",SP,"_",receiver),"receiver_act_pv.adj"] <- smy_act_comb_sub[,"act_pv.adj"]
  }

  ccc <- ccc[!ccc[,"sender"]==ccc[,"receiver"],]

  if(ncol(ccc)>0)
  {
    ccc[,"overall_strength"] <- ccc[,"sender_exp_logFC"] * ccc[,"receiver_act_diff"]

    library(metap)
    ccc[ccc[,"sender_exp_pv"]==0,"sender_exp_pv"] <- .Machine$double.xmin # sumlog requires non-zero

    ccc[,"overall_pv"] <- apply(ccc[,c("sender_exp_pv","receiver_act_pv")],1,function(x) sumlog(x)$p)
    ccc[,"overall_pv.adj"] <- p.adjust(ccc[,"overall_pv"], method="BH")
    ccc <- ccc[ccc[,"overall_pv.adj"]<padj_cutoff,]

    ccc <- ccc[order(ccc[,"overall_pv.adj"]),]
  }

  Seurat_obj @misc $SecAct_output $SecretedProteinCCC  <- ccc

  Seurat_obj
}


#' @title Survival regression
#' @description Calculate the risk score of each secreted protein.
#' @param mat Activity matrix.
#' @param surv Survival matrix.
#' @return A matrix.
#' @rdname SecAct.coxph.regression
#' @export
#'
SecAct.coxph.regression <- function(mat, surv)
{
  mat <- t(mat)

  olp <- intersect(rownames(surv),rownames(mat))
  X_olp <- surv[olp,,drop=F]
  Y_olp <- mat[olp,,drop=F]

  smy <- data.frame()
  for(gene in colnames(Y_olp))
  {
    comb <- cbind(X_olp, Act=Y_olp[,gene])

    library(survival)
    errflag <- F
    coxmodel_fit <- tryCatch(
      coxph(Surv(Time, Event) ~ ., data = comb),

      error = function(e){
        errflag <<- T
      },

      warning = function(w){
        errflag <<- T
      }
    )

    smy[gene,"risk score z"] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
    smy[gene,"p value"] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","Pr(>|z|)"])
  }
  as.matrix(smy)
}

