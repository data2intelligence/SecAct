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
  alias2symbol <- read.csv(system.file("extdata", 'NCBI_20251008_gene_result_alias2symbol.csv.gz', package = 'SecAct'),as.is=T)
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


calWeights <- function(spotCoordinates, radius, sigma=100, diagAsZero=TRUE)
{
  nn_result <- RANN::nn2(spotCoordinates, searchtype="radius", radius=radius, k=nrow(spotCoordinates))

  neighbor_indices <- nn_result$nn.idx
  neighbor_distances <- nn_result$nn.dists

  i <- rep(1:nrow(neighbor_indices), each=ncol(neighbor_indices)) # row indices (cell index)
  j <- as.vector(t(neighbor_indices))
  x <- as.vector(t(neighbor_distances))

  valid <- x<=radius & x>0
  i <- i[valid]          # Keep only valid indices
  j <- j[valid]          # Valid neighbor indices
  x <- x[valid]          # Valid distances

  # transform distance to weight
  x <- exp( -x^2 / (2*sigma^2) )

  # Create the sparse matrix using the 'i', 'j', and 'x' vectors
  library(Matrix)
  W <- sparseMatrix(i=i, j=j, x=x, dims=c(nrow(neighbor_indices), nrow(neighbor_indices)), repr="T")
  rownames(W) <- rownames(spotCoordinates)
  colnames(W) <- rownames(spotCoordinates)

  if(diagAsZero==FALSE) diag(W) <- 1

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

expand_rows <- function(mat)
{
  new_rows <- lapply(1:nrow(mat), function(i) {
    names <- strsplit(rownames(mat)[i], "\\|")[[1]]
    do.call(rbind, replicate(length(names), mat[i, , drop = FALSE], simplify = FALSE)) |>
      `rownames<-`(names)
  })
  do.call(rbind, new_rows)
}
