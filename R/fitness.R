# source("./R/cluster.R")
# source("./R/dist.R")

fitness.cor <- function(individual, penalty.function = NULL, k, data, known_lfc, method, ...){

  dims <- length(individual)/k
  fitness.value <- NA

  m.individual <- matrix(individual, nrow = k, ncol = dims)

  which.dists <- dist.corr(data, d2=m.individual, method = "pearson")

  which.dists <- as.factor(which.dists)

  # to avoid a convergence for configurations different from user-specified k
  if (length(unique(which.dists)) < k) {
    # maximum penalty
    fitness.value = -1

  } else {
    cluster_lfc <- compute.l2fc(data, which.dists)
    fitness.value <- abs(cor(cluster_lfc, known_lfc))
    cat(sprintf("Fitness: %f\n", fitness.value))
  }

  return(fitness.value)
}
