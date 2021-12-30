dist.corr <- function(data, d2=NULL, method=c("pearson", "kendall", "spearman", "euclidean", "manhattan")){

  if (is.null(data)){
    stop("'dataset' can not be NULL")
  }

  if (class(data) != 'data.frame'){
    stop("'dataset' must be a data.frame object.")
  }

  method <- match.arg(method)

  if(!is.null(d2)){
    if(method %in% c("pearson", "spearman", "kendall")){
      dists <- 1 - cor(t(data), t(d2), method = method)
    }
    else{
      dists <- Rfast::dista(data, d2, type = method, square = TRUE)
    }
    return(apply(dists, 1, which.min))
  }
  else{
    if(method %in% c("pearson", "spearman", "kendall")){
      dists <- 1 - cor(t(data), method = method)
    }
    else{
      dists <- dist(data, method = method, diag = FALSE, upper = FALSE)
    }
    return(dists)
  }
}
