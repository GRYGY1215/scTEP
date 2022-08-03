
clustering <- function(data) {
  cl <- parallel::makeCluster(16, outfile = "/dev/null")
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(tidyverse)
    library(scDHA)
    library(matrixStats)
    library(doParallel)
    library(SummarizedExperiment)
    library(SingleCellExperiment)
    library(Matrix)
    library(foreach)
  })
  allCluster <- NULL
  for(k in 5:10){
    scDHARes <- scDHA(data$expr, k = k, do.clus = T, gen_fil = T, ncores = 10, seed = 1)
    allCluster = c(allCluster, list(scDHARes$cluster))
  }
  allCluster <- lapply(allCluster, as.integer)
  parallel::stopCluster(cl)
  return(allCluster)
}


