#' @title Cell-cell communication heatmap
#' @description Draw a heatmap of cell-cell communicationss mediated by secreted proteins.
#' @param data A SpaCET object or a Seurat object.
#' @param row.sorted Whether to sort rows.
#' @param column.sorted Whether to sort columns.
#' @param colors_cellType Colors for cell types.
#' @return A ggplot2 object.
#' @rdname SecAct.CCC.heatmap
#' @export
#'
SecAct.CCC.heatmap <- function(data, row.sorted=FALSE, column.sorted=FALSE, colors_cellType)
{
  if(class(data)[1]=="SpaCET")
  {
    ccc <- data @results $SecAct_output $SecretedProteinCCC
  }
  if(class(data)[1]=="Seurat")
  {
    ccc <- data @misc $SecAct_output $SecretedProteinCCC
  }

  ccc <- cbind(ccc, communication=1)
  ccc <- cbind(ccc, senderReceiver=paste0(ccc[,"sender"],"-",ccc[,"receiver"]))

  ccc <- ccc[!ccc[,"sender"]==ccc[,"receiver"],]

  mat = reshape2::acast( ccc[,c("sender","receiver","communication")], sender~receiver, length, value.var="communication")

  cellTypes <- sort(unique(c(rownames(mat),colnames(mat))))
  for(cellType in cellTypes)
  {
    if(cellType%in%rownames(mat)&cellType%in%colnames(mat))
    {
      mat[cellType,cellType] <- NA
    }
  }

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
  })

  if(row.sorted==TRUE) mat <- mat[order(apply(mat,1,function(x) sum(x,na.rm=T)),decreasing=T),,drop=F]
  if(column.sorted==TRUE) mat <- mat[,order(apply(mat,2,function(x) sum(x,na.rm=T)),decreasing=T),drop=F]

  row_ha <- rowAnnotation(
    Count = anno_barplot(rowSums(mat,na.rm=T),gp = gpar(fill = colors_cellType[rownames(mat)]) )
  )
  column_ha <- columnAnnotation(
    Count = anno_barplot(colSums(mat,na.rm=T),gp = gpar(fill = colors_cellType[colnames(mat)]) ),
    annotation_name_side = "left"
  )

  Heatmap(as.matrix(mat),
     name = "Count",
     col = circlize::colorRamp2(c(-50, 0,50), c("green", "white", "red")),
     row_names_side = "left",
     column_title_side = "bottom",
     row_title = "Sender",
     column_title = "Receiver",
     show_heatmap_legend = FALSE,
     top_annotation = column_ha,
     right_annotation = row_ha,
     cluster_rows = FALSE,
     cluster_columns = FALSE,
     cell_fun = function(j, i, x, y, width, height, fill) {grid.text(mat[i, j], x, y)}
  )
}


#' @title Cell-cell communication circle plot
#' @description Draw a circle plot of cell-cell communicationss mediated by secreted proteins.
#' @param data A SpaCET object or a Seurat object.
#' @param colors_cellType Colors for cell types.
#' @return A ggplot2 object.
#' @rdname SecAct.CCC.circle
#' @export
#'
SecAct.CCC.circle <- function(data, colors_cellType)
{
  if(class(data)[1]=="SpaCET")
  {
    ccc <- data @results $SecAct_output $SecretedProteinCCC
  }
  if(class(data)[1]=="Seurat")
  {
    ccc <- data @misc $SecAct_output $SecretedProteinCCC
  }

  ccc <- cbind(ccc, communication=1)
  ccc <- cbind(ccc, senderReceiver=paste0(ccc[,"sender"],"-",ccc[,"receiver"]))

  ccc <- ccc[!ccc[,"sender"]==ccc[,"receiver"],]


  mat = reshape2::acast( ccc[,c("sender","receiver","communication")], sender~receiver, length, value.var="communication")

  for(cellType in cellTypes)
  {
    if(cellType%in%rownames(mat)&cellType%in%colnames(mat))
    {
      mat[cellType,cellType] <- NA
    }
  }

  suppressPackageStartupMessages({
    library(circlize)
  })

  mat = log(mat+1)

  chordDiagram(
    mat,
    directional = 1,
    grid.col = colors_cellType,
    annotationTrack = c("name", "grid"),
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow"
  )
}


#' @title Draw a heat map
#' @description Draw a heat map of secreted proteins.
#' @param fg.mat A matrix of values.
#' @param title The title for plot.
#' @param colors Colors.
#' @return A ggplot2 object.
#' @rdname SecAct.heat.map
#' @export
#'
SecAct.heat.map <- function(fg.mat, title=NULL, colors=c("#03c383","#aad962","#fbbf45","#ef6a32"))
{
  fg.df <- reshape2::melt(fg.mat)
  colnames(fg.df)[3] <- "Activity"

  fg.df[["Var1"]] <- factor(fg.df[["Var1"]], levels=rev(rownames(fg.mat)))
  fg.df[["Var2"]] <- factor(fg.df[["Var2"]], levels=colnames(fg.mat))

  library(ggplot2)
  ggplot(fg.df, aes(Var2, Var1, fill=Activity)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colours = colors)+
    ggtitle(title)+
    theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_blank()
    )
}


#' @title Draw a bar plot
#' @description Draw a bar plot of secreted proteins.
#' @param fg.vec A vector of values.
#' @param title The title for plot.
#' @param colors Colors.
#' @return A ggplot2 object.
#' @rdname SecAct.bar.plot
#' @export
#'
SecAct.bar.plot <- function(fg.vec, title=NULL, colors=c("#91bfdb","#fc8d59"))
{
  fg.df <- data.frame(gene=names(fg.vec), value=fg.vec)
  fg.df <- cbind(fg.df, dir=ifelse(fg.vec<0,"down","up"))
  fg.df <- cbind(fg.df, y=ifelse(fg.vec<0,0.1,-0.1))
  fg.df <- cbind(fg.df, hjust=ifelse(fg.vec<0,0,1))
  fg.df[["gene"]] <- factor(fg.df[["gene"]], levels=names(sort(fg.vec)))

  library(ggplot2)
  ggplot(fg.df, aes(gene, value, label=gene)) +
    geom_col(aes(fill=dir), width = .88, color = "white") +
    geom_text(aes(y = y, hjust=hjust), angle = 0) +
    scale_fill_manual(values=colors)+
    geom_hline(yintercept=0)+
    ggtitle(title)+
    theme_classic()+
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(color="black", vjust=0.5),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "none"
    )+
    coord_flip()+
    scale_y_continuous(position = "right")
}


#' @title Draw a lollipop plot
#' @description Draw a lollipop plot of secreted proteins.
#' @param fg.vec A vector of values.
#' @param title The title for plot.
#' @return A ggplot2 object.
#' @rdname SecAct.lollipop.plot
#' @export
#'
SecAct.lollipop.plot <- function(fg.vec, title=NULL)
{
  fg.df <- data.frame(gene=names(fg.vec), value=fg.vec)
  fg.df <- cbind(fg.df, dir=ifelse(fg.vec<0,"down","up"))
  fg.df <- cbind(fg.df, y=ifelse(fg.vec<0,0.1,-0.1))
  fg.df <- cbind(fg.df, hjust=ifelse(fg.vec<0,0,1))
  fg.df[["gene"]] <- factor(fg.df[["gene"]], levels=names(sort(fg.vec)))

  library(ggplot2)
  ggplot(fg.df, aes(gene, value, label=gene)) +
    geom_segment(aes(x = gene, xend = gene, y = 0, yend = value), color = "grey") +
    geom_point(color = "#619CFF", size=3)+
    geom_text(aes(y = y, hjust=hjust), angle = 0) +
    geom_hline(yintercept=0)+
    ggtitle(title)+
    theme_classic()+
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(color="black", vjust=0.5),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "none"
    )+
    coord_flip()+
    scale_y_continuous(position = "right")
}


#' @title Draw a survival plot
#' @description Draw a survival plot of secreted proteins.
#' @param mat Activity matrix.
#' @param surv Survival matrix.
#' @param gene Gene symbol.
#' @param x.title Title for x axis.
#' @return A ggplot2 object.
#' @rdname SecAct.survival.plot
#' @export
#'
SecAct.survival.plot <- function(mat, surv, gene, x.title="Time")
{
  mat <- t(mat)

  olp <- intersect(rownames(surv),rownames(mat))
  X_olp <- surv[olp,,drop=F]
  Y_olp <- mat[olp,,drop=F]


  data <- Y_olp[,gene,drop=F]
  survival <- X_olp
  margin <- 5
  # align matrix names
  common = Reduce(intersect, list(rownames(data),rownames(survival)))
  # sprintf("%s samples", length(common))

  data = data[common,,drop=F]
  survival = survival[common,,drop=F]

  # stop at low death rate
  death_rate = sum(survival[,2])/dim(survival)[1]
  if(length(death_rate) < 0.1) q()

  # split up survival and background
  surv = Surv(survival[,1], survival[,2])

  if(dim(survival)[2] > 2){
    B = survival[,3:dim(survival)[2], drop=F]
  }else{
    B = survival[,c(), drop=F]
  }

  # build up regression data space
  B = cbind(B, rep(0, dim(data)[1]))
  B = as.data.frame(B)
  N_B = ncol(B)
  colnames(B)[N_B] = "pivot"

  # iterate over features
  features = colnames(data)
  N = length(features)

  result = NULL

  step = round(max(N/100,1))
  for (i in 1:N)
  {
    # progress report
    if((i %% step) == 0){
      sprintf("%s", round(100 * i/N, 2))
    }

    fid = features[i]

    # part 1: overall regression
    arr = B[,N_B] = data[,i]

    arr_result = CoxPH_best_separation(B, surv, margin)

    if(sum(is.na(arr_result)) > 0){
      warning(paste0('Jump with failed continuous regression ', fid))
      next
    }

    if(is.null(result)){
      result = matrix(nrow = N, ncol=length(arr_result))
      colnames(result) = names(arr_result)
      rownames(result) = features
    }

    result[i,] = arr_result
  }

  mean.value = colMeans(data)
  result = cbind(result, mean.value)

  N = rep(length(common), dim(result)[1])
  result = cbind(result, N)

  cutoff <- round(result[gene,"thres.opt"],3)

  groups <- as.character(data[,gene]>cutoff)

  library("survminer")
  surv.df <- cbind(survival,groups)
  fit <- survfit(Surv(Time, Event) ~ groups, data = surv.df)

  ggsurvplot(fit,
    data=surv.df,
    palette = c("blue","red"),
    legend.labs = c(paste0("Low (n=",sum(groups=="FALSE"),")"),paste0("High (n=",sum(groups=="TRUE"),")"))
    )$plot+
    xlab(x.title)+
    ylab("Percentage")
}


CoxPH_best_separation = function(X, Y, margin)
{
  # part 1: continuous regression
  errflag = F

  coxph.fit = tryCatch(
    coxph(Y~., data=X),
    error = function(e) errflag <<- T,
    warning = function(w) errflag <<- T)

  if(errflag) return (NA)

  n_r = nrow(X)
  n_c = ncol(X)

  arr_result = summary(coxph.fit)$coef[n_c,]

  if(is.na(arr_result["z"])) return (NA)

  # no need to find optimal threshold
  if(is.null(margin)) return (arr_result)

  # part 2: find the optimal threshold
  arr = X[, n_c]
  vthres = sort(arr)

  # these are missing values, not NULL not existing values
  zscore_opt = thres_opt = NA

  for(i in (margin+1):(n_r-margin))
  {
    X[, n_c] = as.numeric(arr >= vthres[i])

    errflag = F
    coxph.fit = tryCatch(
      coxph(Y~., data=X),
      error = function(e) errflag <<- T,
      warning = function(w) errflag <<- T)

    if(errflag) next

    z = summary(coxph.fit)$coef[n_c, "z"]
    if(is.na(z)) next

    if (is.na(zscore_opt)){
      zscore_opt = z
      thres_opt= vthres[i]

    }else if(arr_result['z'] > 0){
      if(z > zscore_opt){
        zscore_opt = z
        thres_opt = vthres[i]
      }

    }else{ # arr_result['z'] <= 0
      if(z < zscore_opt){
        zscore_opt = z
        thres_opt = vthres[i]
      }
    }
  }

  arr_result['thres.opt'] = thres_opt
  arr_result['z.opt'] = zscore_opt

  return (arr_result)
}
