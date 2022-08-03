getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

cluster_dist = function(latent, data){
  cl = data$cl
  label = data$label
  cell.stages = data$stages
  start = data$startStage

  #Careful about cl is numeric factor or character factor
  latent = cbind(latent, cl)
  latent = as.data.frame(latent)
  cell_stages = factor(cl)
  #Calculate R according to distance to starting group
  x = latent
  total_col_num = ncol(x)
  #Euclidian center
  clus_cent =NULL
  #Select the last column, cell type/clustering result, which has to be represent by numeric id.
  for (i in 1:length(unique(x[,total_col_num]))){
    temp = x[which(x[,total_col_num] == i),]
    temp_cent = apply(temp, 2, mean)
    clus_cent = rbind(clus_cent, temp_cent)
  }

  #Move the start cluster to the first row.
  #Load clustering results from last column.
  #getmode, find the cluster(data[,total_col_num]) that has the most start cell type.
  #Think about it, should I try to select start cluster by percentage of start cell.
  start_clus = getmode(latent[,total_col_num][which(label==start)])
  clus_cent_row_names = paste("clus", seq(1:length(unique(x[,total_col_num]))))
  #Number of clusters need to be larger than 2.
  clus_cent = rbind(clus_cent[start_clus, ] , clus_cent[-start_clus,])
  clus_cent_row_names = c(clus_cent_row_names[start_clus], clus_cent_row_names[-start_clus])
  rownames(clus_cent) = clus_cent_row_names

  clus_dist = as.matrix(dist(clus_cent[,1:(total_col_num -1)], method = "euclidian"))[,1]
  clus_cent = as.data.frame(clus_cent)

  dist_order =NULL
  for (i in 1:length(latent[,total_col_num])){
    dist_order = c(dist_order, clus_dist[which(clus_cent$cl == latent[,total_col_num][i])])
  }

  out = list(center = clus_cent, clus_dist = clus_dist, dist_order = dist_order)
}

sort_igraph_object = function(g, dist_vector, start_point){
  #plot(g)
  next_nodes = g %>% igraph::neighbors(start_point) %>% names()

  if (length(next_nodes) == 0){
    return(list(g = g, s = start_point))
  }else{
    g_vertices = igraph::V(g) %>% names()
    if (dist_vector[start_point] != min(dist_vector[g_vertices])){

      min_node = names(dist_vector)[dist_vector == min(dist_vector[g_vertices])]

      # 互相交换两个节点的边，如果它们相邻的话
      bool_adjancent = igraph::are_adjacent(g, min_node, start_point)
      if (bool_adjancent && length(g_vertices) == 2){
        start_point = min_node
      }else{
        # 这里 为除了节点本身的临节点
        neighbour_nodes_start = g %>% igraph::neighbors(start_point) %>% names()
        neighbour_nodes_start = neighbour_nodes_start[!neighbour_nodes_start == min_node]
        neighbour_nodes_min = g %>% igraph::neighbors(min_node) %>% names()
        neighbour_nodes_min = neighbour_nodes_min[!neighbour_nodes_min == start_point]

        # 这里没有考虑到互相交换的两个node相邻的情况
        g = g %>% igraph::delete_vertices(c(start_point, min_node))

        edges_start = start_point %>% rep(length(neighbour_nodes_min)) %>%
          rbind(neighbour_nodes_min) %>% c()
        edges_min = min_node %>% rep(length(neighbour_nodes_start)) %>%
          rbind(neighbour_nodes_start) %>% c()
        if (bool_adjancent){
          edges_start = edges_start %>% c(c(start_point, min_node))
        }
        g = g %>% igraph::add_vertices(2, name = c(start_point, min_node)) %>%
          add_edges(edges_start) %>%
          #add_vertices(1, name = edges_min) %>%
          add_edges(edges_min)
        start_point = min_node
      }

    }
    next_nodes = g %>% igraph::neighbors(start_point) %>% names()
    next_g = g %>% igraph::delete_vertices(start_point) %>% igraph::decompose.graph()

    sorted_sub_graph = sort_igraph_object(next_g[[1]], dist_vector, intersect(next_nodes, igraph::V(next_g[[1]]) %>% names()))
    sorted_g = sorted_sub_graph$g %>% igraph::add_vertices(1, name = start_point) %>%
      igraph::add_edges(c(start_point, sorted_sub_graph$s))
    if (length(next_g) > 1){
      for (id in 2:length(next_g)){
        temp_sub_graph = sort_igraph_object(next_g[[id]], dist_vector, intersect(next_nodes, igraph::V(next_g[[id]]) %>% names()))
        temp_sub_graph$g = temp_sub_graph$g %>% igraph::add_vertices(1, name = start_point) %>%
          igraph::add_edges(c(start_point, temp_sub_graph$s))
        sorted_g = igraph::union(sorted_g, temp_sub_graph$g)
      }
    }
  }
  return(list(g = sorted_g, s = start_point))
}

#Load data
preprocessing <- function(data){
  expr <- as.matrix(Matrix::Matrix(t(SummarizedExperiment::assay(data)), sparse = TRUE))
  expr <- expr[, colSums(expr) > 0]
  colnames(expr) <- toupper(colnames(expr))
  if (max(expr) > 100) expr <- log2(expr + 1)

  data = list(
    expr = expr,
    label = as.character(data$label),
    startStage = as.character(data@metadata$cell.stages[1]),
    stages = data@metadata$cell.stages
  )
  return(data)
}

fa = function(data, genesets, data_org = 'hsa'){
  cl <- parallel::makeCluster(16, outfile = "/dev/null")
  doParallel::registerDoParallel(cl, cores = 8)
  parallel::clusterEvalQ(cl, {library(psych)})

  datByGeneset <- lapply(genesets[[data_org]], function(genes){
    commonGenes <- intersect(genes, colnames(data$expr))
    if (length(commonGenes) < 10) return(NULL)
    as.matrix(data$expr[, commonGenes])
  })

  datByGeneset <- datByGeneset[!sapply(datByGeneset, is.null)]

  faData <- foreach(dat = datByGeneset) %dopar% {
    r <- suppressMessages(suppressWarnings(psych::fa(dat, nfactors = 2, warnings = F)$scores))
    r[r > 5] <- 5
    r[r < -5] <- -5
    r
  } %>% do.call(what = cbind) %>% `-`(min(.))
  parallel::stopCluster(cl)
  return(faData)
}


# id here is the name of dataset without '.rds'
trajectoryinference<- function(data, data_org, fa, allCluster, ncores = 10L, seed = NULL) {

  scDHA_res <- scDHA(fa, do.clus = T, gen_fil = T, ncores = 16, seed = 1)

  tmp <- rep(0, nrow(data$expr))
  latent = scDHA_res$latent
  df <- as.data.frame(latent)

  for (clus in allCluster) {
    # Distance method
    data$cl = clus
    out = cluster_dist(latent, data)
    pseudotime_tmp = out$dist_order
    pseudotime_tmp <- pseudotime_tmp/max(pseudotime_tmp)
    tmp <- tmp + pseudotime_tmp

  }
  pseudotime = tmp

  # process for graph
  cl <- scDHA_res$cluster
  # Assign clusters with only 1 cells to closest group
  clus_id_remove = which(table(cl) == 1)
  for (i in clus_id_remove){
    temp_index = which(cl == i)
    close_cell_index = which.min(abs(pseudotime_tmp[-temp_index]-pseudotime_tmp[temp_index]))
    if (close_cell_index >= temp_index){close_cell_index = close_cell_index +1}
    cl[temp_index] = cl[close_cell_index]
  }

  for (i in 1:length(clus_id_remove)){
    temp_clus_id = clus_id_remove[i]
    cl[which(cl>temp_clus_id)] = cl[which(cl>temp_clus_id)] - 1
    clus_id_remove = clus_id_remove-1
  }

  data_clus_cent = NULL
  for (cl_id in unique(cl)){
    cl_index = which(cl == cl_id)
    clus_pseudotime = mean(pseudotime_tmp[cl_index])
    tmp_clus_cent = scDHA_res$latent[cl_index, ] %>%  colMeans() %>% c(clus_pseudotime)
    data_clus_cent = rbind(data_clus_cent, tmp_clus_cent)
  }

  data_clus_cent = data_clus_cent %>% as.matrix()
  rownames(data_clus_cent) = unique(cl) %>% as.character()


  distM <- data_clus_cent[, 1:15] %>% dist() %>% as.matrix() # distance matrix for cluster centers

  # Build a tree consistent with pseudotime 4
  sorted_graph <- igraph::graph_from_adjacency_matrix(distM, weighted = TRUE, mode = "undirected") %>%
    igraph::mst() %>% sort_igraph_object(data_clus_cent[, 16], as.character(1))

  # assign weight according to distM
  E(sorted_graph$g)$weight = sorted_graph$g %>% get.edgelist() %>% apply(1, function(x) distM[as.numeric(x)[1], as.numeric(x)[2]])

  # Calculate the length of edge according to the euclidean distance.
  # mst %>% igraph::as_data_frame()
  edges = sorted_graph$g %>% get.edgelist()
  edges = edges[nrow(edges):1, ]
  for (i in 1:nrow(edges)){
    start = edges[i, 1]
    end = edges[i, 2]
    if (data_clus_cent[, 16][start] > data_clus_cent[, 16][end]){
      edges[i, 1] = end
      edges[i, 2] = start
    }
  }
  edge_length = E(sorted_graph$g)$weight %>% rev()

  milestone_network <- tibble::tibble(
    from = as.character(edges[,1]),
    to = as.character(edges[,2]),
    length = edge_length,
    directed = TRUE
  )

  out = list(
    pseudotime = pseudotime,
    cluster = cl,
    data_clus_cent = data_clus_cent,
    milestone_network = milestone_network,
    sorted_g = sorted_graph
  )
  return(out)
}
