#' @title Secreted protein signaling pattern
#' @description Calculate the signaling pattern of secreted proteins based on their activities.
#' @param SpaCET_obj A SpaCET object.
#' @param k Number of patterns for NMF.
#' @return A ggplot2 object.
#' @examples
#' SpaCET_obj <- SecAct.signaling.pattern(SpaCET_obj, k=3)
#'
#' @rdname SecAct.signaling.pattern
#' @export
#'
SecAct.signaling.pattern <- function(SpaCET_obj, k=3)
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }
  if(is.null(SpaCET_obj @results $SecAct_output $SecretedProteinActivity))
  {
    stop("Please run SecAct.activity.inference first.")
  }

  print("Step 1. Filtering")

  act <- SpaCET_obj @results $SecAct_output $SecretedProteinActivity $zscore

  expr <- SpaCET_obj@input$counts
  rownames(expr) <- transferSymbol(rownames(expr))
  expr <- rm_duplicates(expr)

  # normalize to TPM
  stats <- Matrix::colSums(expr)
  expr <- sweep_sparse(expr,2,stats,"/")
  expr@x <- expr@x * 1e5

  # transform to log space
  expr@x <- log2(expr@x + 1)

  weights <- calWeights(colnames(exp), r=3, diag0=TRUE)
  act_new <- act[,colnames(weights)] # remove spot island
  exp_new <- expr[,colnames(weights)] # remove spot island

  exp_new_aggr <- exp_new %*% weights

  corr <- data.frame()
  for(gene in rownames(act))
  {
    act_gene <- act_new[gene,]

    if(gene%in%rownames(exp))
    {
      exp_gene <- exp_new_aggr[gene,]

      cor_res <- cor.test(act_gene, exp_gene, method="pearson")

      corr[gene,"r"] <- cor_res$estimate
      corr[gene,"p"] <- cor_res$p.value
    }else{
      corr[gene,"r"] <- NA
      corr[gene,"p"] <- NA
    }
  }
  corr <- cbind(corr, padj=p.adjust(corr[,"p"], method="BH") )
  corr_genes <- rownames(corr[!is.na(corr[,"r"])&corr[,"r"]>0.1&corr[,"padj"]<0.01,])

  print(paste0(length(corr_genes),"/",nrow(act)," secreted proteins are kept to infer signaling patterns."))


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


#' @title Secreted protein signaling velocity
#' @description Calculate the signaling velocity of secreted proteins based on their activities.
#' @param SpaCET_obj A SpaCET object.
#' @param gene Gene symbol coding a secreted protein.
#' @param signalMode Mode of signaling velocity, i.e., "receiving", "sending", and "both".
#' @return A ggplot2 object.
#' @details
#' The velocity direction starts from the source cell producing a secreted protein and moves to sink cells receiving the secreted protein signal. The velocity magnitude represents the product between the secreted protein-coding gene expression at source cells and signaling activities at sink cells.
#'
#' @examples
#' SecAct.signaling.velocity.spotST(SpaCET_obj, gene="TGFB1", signalMode="receiving")
#' SecAct.signaling.velocity.spotST(SpaCET_obj, gene="TGFB1", signalMode="sending")
#' SecAct.signaling.velocity.spotST(SpaCET_obj, gene="TGFB1", signalMode="both")
#'
#' @rdname SecAct.signaling.velocity.spotST
#' @export
#'
SecAct.signaling.velocity.spotST <- function(
  SpaCET_obj,
  gene,
  signalMode="both",
  animated=FALSE
)
{
  if(class(SpaCET_obj)!="SpaCET")
  {
    stop("SpaCET object is requried.")
  }

  counts <- SpaCET_obj@input$counts
  rownames(counts) <- transferSymbol(rownames(counts))
  counts <- rm_duplicates(counts)

  exp <- sweep(counts, 2, Matrix::colSums(counts), "/") *1e5
  exp <- log2(exp+1)

  act <- SpaCET_obj@results$SecretedProteinActivity$zscore
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

    startend2[startend2[,"vec_len"]<0.1,"vec_len"] <- 0.001
    startend2[startend2[,"vec_len"]>=0.1,"vec_len"] <- 0.08

    #startend2[,"vec_len"] <- startend2[,"vec_len"]/max(startend2[,"vec_len"])
    #startend2[,"vec_len"] <- startend2[,"vec_len"]/10
  }

  if(animated==TRUE)
  {
    startend2_temp <- startend2
    startend2_temp[,"x_end"] <- startend2_temp[,"x_start"]
    startend2_temp[,"y_end"] <- startend2_temp[,"y_start"]

    startend2 <- rbind(
      cbind(startend2_temp,tim=1),
      cbind(startend2,tim=2)
    )
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
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), color="black", data=startend2, arrow = arrow(length = unit(startend2$vec_len, "cm")))+
    #+
    ggtitle(paste0(gene," (receiving)"))+
    theme_void()+
    theme(
      plot.title = element_text(hjust = 0.5)
    )+coord_flip()

  if(animated==TRUE)
  {
    library(gganimate)
    p2 <- animate(p2+transition_time(tim),nframes=15)
  }

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

transferSymbol <- function(x)
{
  alias2symbol <- read.csv(system.file("extdata", 'NCBI_20230818_gene_result_alias2symbol.csv', package = 'SecAct'),as.is=T)
  alias2symbol[is.na(alias2symbol[,"Alias"]),"Alias"] <- "NA"

  x[x%in%alias2symbol[,1]] <- alias2symbol[
    match(
      x[x%in%alias2symbol[,1]],
      alias2symbol[,1]
    ), 2]

  x
}

rm_duplicates <- function(mat){
  dupl <- duplicated(rownames(mat))
  if (sum(dupl) > 0){
    dupl_genes <- unique(rownames(mat)[dupl])
    mat_dupl <- mat[rownames(mat) %in% dupl_genes,,drop=F]
    mat_dupl_names <- rownames(mat_dupl)
    mat <- mat[!dupl,,drop=F]

    for(gene in dupl_genes){
      mat_dupl_gene <- mat_dupl[mat_dupl_names == gene,]
      dupl_sum <- apply(mat_dupl_gene,1,sum)
      max_flag <- which(dupl_sum==max(dupl_sum))
      mat[gene,] <- mat_dupl_gene[max_flag[1],] # in case two values are max
    }
  }
  return(mat)
}



#' @title Cell-cell communication from single cell data
#' @description Calculate condition-specific cell-cell communication mediated by secreted proteins from scRNA-Seq data.
#' @param data A Seurat object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param condition_meta Column name in meta data that includes condition information.
#' @param conditionCase Case condition.
#' @param conditionControl Control condition.
#' @param act_diff_cutoff Cut off for activity change (i.e., z score) in step 1.
#' @param exp_logFC_cutoff Cut off for log fold change in step 2.
#' @param exp_fraction_case_cutoff Cut off for the fraction of cells expressing secreted protein-coding genes in step 2.
#' @param padj_cutoff Adjusted p value cut off.
#' @return A Seurat object.
#' @rdname SecAct.CCC.scRNAseq
#' @export
#'
SecAct.CCC.scRNAseq <- function(
  data,
  cellType_meta,
  condition_meta,
  conditionCase,
  conditionControl,
  act_diff_cutoff=2,
  exp_logFC_cutoff=0.2,
  exp_fraction_case_cutoff=0.4,
  padj_cutoff=0.01
)
{
  counts <-  data@assays$RNA@counts
  meta <- data@meta.data

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

  data @misc $SecAct_output $SecretedProteinActivity <- SecAct.signaling.inference(bulk.diff)

  print("Step 2: assessing changes in secreted protein expression.")

  for(cellType in cellTypes)
  {
    # case
    expr <- counts[,meta[,condition_meta]==conditionCase&meta[,cellType_meta]==cellType]

    # normalize to TPM
    #expr <- sweep(expr, 2, Matrix::colSums(expr), "/") * 1e5

    stats <- Matrix::colSums(expr)
    expr <- sweep_sparse(expr,2,stats,"/")
    expr@x <- expr@x * 1e5

    # transform to log space
    expr@x <- log2(expr@x + 1)

    expr_case <- expr


    # control
    expr <- counts[,meta[,condition_meta]==conditionControl&meta[,cellType_meta]==cellType]

    # normalize to TPM
    #expr <- sweep(expr, 2, Matrix::colSums(expr), "/") * 1e5

    stats <- Matrix::colSums(expr)
    expr <- sweep_sparse(expr,2,stats,"/")
    expr@x <- expr@x * 1e5

    # transform to log space
    expr@x <- log2(expr@x + 1)

    expr_control <- expr


    smy_deg <- data.frame()
    genes <- intersect(rownames(expr_case), rownames(data @misc $SecAct_output $SecretedProteinActivity $zscore))
    for(gene in genes)
    {
      T_vec <- expr_case[gene,]
      C_vec <- expr_control[gene,]

      wTest <- wilcox.test(T_vec,C_vec)

      smy_deg[gene,"exp_logFC"] <- mean(T_vec) - mean(C_vec)
      smy_deg[gene,"exp_mean_all"] <- mean(c(T_vec,C_vec))
      smy_deg[gene,"exp_mean_case"] <- mean(T_vec)
      smy_deg[gene,"exp_mean_control"] <- mean(C_vec)
      smy_deg[gene,"exp_fraction_case"] <- sum(T_vec>0)/length(T_vec)
      smy_deg[gene,"exp_fraction_control"] <- sum(C_vec>0)/length(C_vec)
      smy_deg[gene,"exp_pv"] <- wTest$p.value
      smy_deg[gene,"exp_pv.adj"] <- wTest$p.value
    }
    smy_deg <- smy_deg[smy_deg[,"exp_mean_all"]>0,]

    smy_deg[,"exp_pv.adj"] <- p.adjust(smy_deg[,"exp_pv"], method="BH")
    smy_deg <- smy_deg[order(smy_deg[,"exp_pv.adj"]),]

    data @misc $SecAct_output $SecretedProteinExpression [[cellType]] <- smy_deg
  }


  print("Step 3: linking sender and receiver cell types.")

  smy_deg_comb <- data.frame()
  smy_act_comb <- data.frame()

  for(cellType in cellTypes)
  {
    smy_deg <- data @misc $SecAct_output $SecretedProteinExpression [[cellType]]

    smy_deg_up <- smy_deg[
      smy_deg[,"exp_logFC"]>exp_logFC_cutoff&
      smy_deg[,"exp_fraction_case"]>exp_fraction_case_cutoff&
      smy_deg[,"exp_pv.adj"]<adjp_cutoff, ]

    if(nrow(smy_deg_up)>0)
    {
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"cellType"] <- cellType
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"secretedProtein"] <- rownames(smy_deg_up)
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"exp_logFC"] <- smy_deg_up[,"exp_logFC"]
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"exp_pv"] <- smy_deg_up[,"exp_pv"]
      smy_deg_comb[paste0(cellType,"_",rownames(smy_deg_up)),"exp_pv.adj"] <- smy_deg_up[,"exp_pv.adj"]
    }

    res <- data @misc $SecAct_output $SecretedProteinActivity

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
      smy_act[,"act_pv.adj"]<adjp_cutoff, ]

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

  ccc[,"overall_strength"] <- ccc[,"sender_exp_logFC"] * ccc[,"receiver_act_diff"]

  library(metap)
  ccc[,"overall_pv"] <- apply(ccc[,c("sender_exp_pv","receiver_act_pv")],1,function(x) sumlog(x)$p)
  ccc[,"overall_pv.adj"] <- p.adjust(ccc[,"overall_pv"], method="BH")
  ccc <- ccc[ccc[,"overall_pv.adj"]<adjp_cutoff,]

  ccc <- ccc[order(ccc[,"overall_pv.adj"]),]

  data @misc $SecAct_output $SecretedProteinCCC  <- ccc

  data
}


#' @title Cell-cell communication from spatial data
#' @description Calculate cell-cell communication mediated by secreted proteins from spatial transcriptomics data.
#' @param data A SpaCET object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param padj_cutoff Adjusted p value cut off.
#' @return A Seurat object.
#' @rdname SecAct.CCC.scST
#' @export
#'
SecAct.CCC.scST <- function(
    data,
    cellType_meta,
    radius=0.02,
    padj_cutoff=0.01
)
{
  counts <-  data@input$counts
  meta <- data@input$spotCoordinates
  coordinate_mat <- meta[,1:2]
  cellType_vec <- meta[,3]

  print("Step 1: ")

  # extract count matrix
  expr <- data@input$counts

  # normalize to TPM
  stats <- Matrix::colSums(expr)
  expr <- sweep_sparse(expr,2,stats,"/")
  expr@x <- expr@x * 1e5

  # transform to log space
  expr@x <- log2(expr@x + 1)

  # normalized with the control samples
  expr.diff <- expr - Matrix::rowSums(expr)

  Sys.time()
  data @results $SecAct_output $SecretedProteinActivity <- SecAct.signaling.inference(expr.diff, sigFilter=TRUE)
  Sys.time()


  print("Step 2: ")

  nn_result <- RANN::nn2(coordinate_mat, k=100, searchtype="radius", radius=radius)

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

  exp <- expr
  act <- data @results $SecAct_output $SecretedProteinActivity$zscore
  act[act<0] <- 0

  act_new <- act[,colnames(weights)] # remove spot island
  exp_new <- exp[,colnames(weights)] # remove spot island

  # Option 1
  exp_new_aggr <- exp_new %*% weights

  corr <- data.frame()
  for(gene in rownames(act))
  {
    act_gene <- act_new[gene,]

    if(gene%in%rownames(exp))
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
  corr <- cbind(corr, padj=p.adjust(corr[,"p"], method="BH"))
  corr_genes <- rownames(corr[!is.na(corr[,"r"])&corr[,"r"]>0&corr[,"padj"]<0.5,])

  print(paste0(length(corr_genes),"/",nrow(act)," secreted proteins are filtered to infer cell-cell communication."))


  print("Step 3: ")

  cellTypes <- unique(cellType_vec)

  cellTypePair1 <- rep(cellTypes, each=length(cellTypes))
  cellTypePair2 <- rep(cellTypes, length(cellTypes))

  cellTypePair <- cbind(cellTypePair1, cellTypePair2)
  cellTypePair <- cellTypePair[cellTypePair1>cellTypePair2,]
  rownames(cellTypePair) <- paste0(cellTypePair[,1],"_",cellTypePair[,2])

  olp <- corr_genes

  ccc <- data.frame()
  for(m in 1:nrow(cellTypePair))
  {
    print(m)
    cellType1 <- cellTypePair[m,1]
    cellType2 <- cellTypePair[m,2]

    cellType1_cells <- which(cellType_vec==cellType1)
    cellType2_cells <- which(cellType_vec==cellType2)


    Tmat <- data.frame(i,j)

    # all cell pair
    Tmat_cellTypePair <- Tmat[Tmat[,"i"]%in%cellType1_cells & Tmat[,"j"]%in%cellType2_cells, ,drop=F]

    n_neighbor <- nrow(Tmat_cellTypePair)

    if(n_neighbor ==0 ) next

    set.seed(123)
    Tmat_background <- data.frame(
      i=sample(cellType1_cells, n_neighbor*1000, replace=T),
      j=sample(cellType2_cells, n_neighbor*1000, replace=T)
    )


    for(SP in olp)
    {
      print(SP)

      # expr(1) * act(2)
      CCC_vec <- expr[SP, Tmat_cellTypePair[,1]] * act[SP, Tmat_cellTypePair[,2]]

      posRatio <- sum(CCC_vec>0)/length(CCC_vec)

      if(posRatio <= 0.05) next

      CCC1000_vec <- expr[SP, Tmat_background[,1]] * act[SP, Tmat_background[,2]]

      CCC_raw <- mean(CCC_vec)
      CCC1000 <- sapply(1:1000, function(x) mean(CCC1000_vec[((x-1)*n_neighbor+1):(x*n_neighbor)]) )

      score <- CCC_raw/mean(CCC1000)
      pv <- (sum(CCC1000>=CCC_raw)+1)/1001

      if(pv < 0.05 & posRatio > 0.05)
      {
        ccc[paste0(cellType1,"_",SP,"_",cellType2),"sender"] <- cellType1
        ccc[paste0(cellType1,"_",SP,"_",cellType2),"secretedProtein"] <- SP
        ccc[paste0(cellType1,"_",SP,"_",cellType2),"receiver"] <- cellType2
        ccc[paste0(cellType1,"_",SP,"_",cellType2),"posRatio"] <- posRatio
        ccc[paste0(cellType1,"_",SP,"_",cellType2),"score"] <- score
        ccc[paste0(cellType1,"_",SP,"_",cellType2),"pv"] <- pv
      }

      if(cellType1==cellType2) next

      # expr(2) * act(1)
      CCC_vec <- expr[SP, Tmat_cellTypePair[,2]] * act[SP, Tmat_cellTypePair[,1]]

      posRatio <- sum(CCC_vec>0)/length(CCC_vec)

      if(posRatio <= 0.05) next

      CCC1000_vec <- expr[SP, Tmat_background[,2]] * act[SP, Tmat_background[,1]]

      CCC_raw <- mean(CCC_vec)
      CCC1000 <- sapply(1:1000, function(x) mean(CCC1000_vec[((x-1)*n_neighbor+1):(x*n_neighbor)]) )

      score <- CCC_raw/mean(CCC1000)
      pv <- (sum(CCC1000>=CCC_raw)+1)/1001

      if(pv < 0.05 & posRatio > 0.05)
      {
        ccc[paste0(cellType2,"_",SP,"_",cellType1),"sender"] <- cellType2
        ccc[paste0(cellType2,"_",SP,"_",cellType1),"secretedProtein"] <- SP
        ccc[paste0(cellType2,"_",SP,"_",cellType1),"receiver"] <- cellType1
        ccc[paste0(cellType2,"_",SP,"_",cellType1),"posRatio"] <- posRatio
        ccc[paste0(cellType2,"_",SP,"_",cellType1),"score"] <- score
        ccc[paste0(cellType2,"_",SP,"_",cellType1),"pv"] <- pv
      }

    }

  }

  ccc[,"pv.adj"] <- p.adjust(ccc[,"pv"], method="BH")
  ccc <- ccc[ccc[,"pv.adj"]<adjp_cutoff,]

  ccc <- ccc[order(ccc[,"pv.adj"]),]

  data @results $SecAct_output $SecretedProteinCCC  <- ccc

  data
}


#' @title Survival regression
#' @description Calculate the risk score of each secreted protein.
#' @param mat Activity matrix.
#' @param surv Survival matrix.
#' @return A matrix.
#' @rdname SecAct.coxph
#' @export
#'
SecAct.coxph <- function(mat, surv)
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

    smy[gene,"risk z"] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","z"])
    smy[gene,"p value"] <- ifelse(errflag, NA, summary(coxmodel_fit)$coefficients["Act","Pr(>|z|)"])
  }
  as.matrix(smy)
}
