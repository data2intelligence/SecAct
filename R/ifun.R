sweep_sparse <- function(m, margin, stats, fun)
{
  f <- match.fun(fun)

  if(margin==1)
  {
    idx <- m@i + 1
  }else{
    if(class(m)[1]=="dgCMatrix")
    {
      idx <- rep(1:m@Dim[2], diff(m@p))
    }else{
      idx <- x@j + 1
    }
  }

  m@x <- f(m@x, stats[idx])
  m
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

rm_duplicates <- function(mat)
{
  dupl <- duplicated(rownames(mat))
  if (sum(dupl) > 0){
    dupl_genes <- unique(rownames(mat)[dupl])
    mat_dupl <- mat[rownames(mat) %in% dupl_genes,,drop=F]
    mat_dupl_names <- rownames(mat_dupl)
    mat <- mat[!dupl,,drop=F]

    for(gene in dupl_genes){
      mat_dupl_gene <- mat_dupl[mat_dupl_names == gene,]
      dupl_sum <- Matrix::rowSums(mat_dupl_gene)
      max_flag <- which(dupl_sum==max(dupl_sum))
      mat[gene,] <- mat_dupl_gene[max_flag[1],] # in case two values are max
    }
  }
  return(mat)
}

rm_duplicates_sparse <- function(mat)
{
  gene_count <- table(rownames(mat))
  gene_dupl <- names(gene_count)[gene_count>1]

  if(length(gene_dupl) > 0){
    gene_unique <- names(gene_count)[gene_count==1]
    gene_unique_index <- which(rownames(mat)%in%gene_unique)

    gene_dupl_index <- c()
    for(gene in gene_dupl)
    {
      gene_dupl_index_gene <- which(rownames(mat)%in%gene)
      mat_dupl_gene <- mat[gene_dupl_index_gene,]
      dupl_sum <- Matrix::rowSums(mat_dupl_gene)
      max_flag <- which(dupl_sum==max(dupl_sum))
      gene_dupl_index <- c(gene_dupl_index,gene_dupl_index_gene[max_flag[1]])
    }

    mat <- mat[sort(c(gene_unique_index,gene_dupl_index)),]
  }

  return(mat)
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

expand_rows <- function(mat) {
  new_rows <- lapply(1:nrow(mat), function(i) {
    names <- strsplit(rownames(mat)[i], "\\|")[[1]]
    do.call(rbind, replicate(length(names), mat[i, , drop = FALSE], simplify = FALSE)) |>
      `rownames<-`(names)
  })
  do.call(rbind, new_rows)
}
