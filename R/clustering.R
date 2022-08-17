
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
  allCluster = foreach (k = 5:10) %dopar% {
    tmp = (function(){
      scDHARes <- scDHA(data$expr, k = k, do.clus = T, gen_fil = T, ncores = 10, seed = 1)
      scDHARes$cluster
    })()
  }
  parallel::stopCluster(cl)
  allCluster <- lapply(allCluster, as.integer)
  return(allCluster)
}







