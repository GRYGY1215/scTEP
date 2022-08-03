# Draw the figure for ground truth
Draw_figure_pseudotime = function(dataset, data, pseudotime, r, cols){
  data = as.data.frame(data)
  # Calculate cluster centroid and create dataframe saving start and end points.
  figure <- ggplot2::ggplot(data, ggplot2::aes(x = V1, y = V2)) +
    ggplot2::labs(x = paste0("UMAP1"), y = paste0("UMAP2"),title = paste("R =", r)) +
    ggplot2::geom_point(ggplot2::aes(color = pseudotime), size = 2) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
          plot.title = ggplot2::element_text(size=20),
          axis.text= ggplot2::element_text(size=20),
          axis.title=ggplot2::element_text(size=20),
          legend.key.size =ggplot2::unit(1, 'cm'), legend.text=ggplot2::element_text(size=20), legend.title=ggplot2::element_text(size=20))
  return(figure)
}

# Draw the figure for ground truth
Draw_figure_lable = function(dataset, data, label, cell.stages, r, cols){
  cell_stages = factor(label, levels = cell.stages)
  label_2 = as.numeric(factor(label, levels = cell.stages))
  data = cbind(data, label_2)
  data = as.data.frame(data)

  # Calculate cluster centroid and create dataframe saving start and end points.
  clus_center = lapply(1:length(unique(label_2)), function(cl) colMeans(data[label_2 == cl, ])) %>%
    do.call(what = rbind)
  colnames(clus_center) = c("x","y", 'cluster_id')

  figure <- ggplot2::ggplot(data, ggplot2::aes(x = V1, y = V2, color = cell_stages)) +
    ggplot2::geom_point()+
    # Plot the centroid
    ggplot2::geom_point(data = as.data.frame(clus_center), ggplot2::aes(x= x,y= y), color = 'black', size = 2)+
    ggplot2::labs(x = paste0("UMAP1"), y = paste0("UMAP2"),title = paste("R =", r)) +
    ggplot2::theme_classic()  +
    ggplot2::scale_color_manual(values = cols) +
    #annotate(geom = "table", x=min(tsne_original$V1),y=min(tsne_original$V2), label = silhou)+
    # scale_color_npg()+
    ggsci::scale_fill_npg() +
    ggplot2::theme(legend.position = "right",
          plot.title = ggplot2::element_text(size=20),
          axis.text= ggplot2::element_text(size=20),
          axis.title=ggplot2::element_text(size=20),
          legend.key.size =ggplot2::unit(1, 'cm'), legend.text=ggplot2::element_text(size=20), legend.title=ggplot2::element_text(size=20))
  return(figure)
}

#Draw the plots
Draw_figure_clus = function(dataset, data, cluster, cell.stages, edge_nodes, r, cols){
  cell_stages = factor(cluster)

  data = cbind(data, cluster)
  data = as.data.frame(data)

  # Calculate cluster centroid and create dataframe saving start and end points.
  # the order is 1 to length
  clus_center = lapply(1:length(unique(cluster)), function(cl) colMeans(data[cluster == cl, ])) %>%
    do.call(what = rbind)

  colnames(clus_center) = c("x","y", 'cluster_id')
  clulines = NULL
  for (edge_index in 1:nrow(edge_nodes)){
    start_cluster = edge_nodes[edge_index, 1] %>% as.numeric()
    end_cluster = edge_nodes[edge_index, 2] %>% as.numeric()
    temp_points = c(clus_center[, 'x'][start_cluster], clus_center[, 'y'][start_cluster],
                    clus_center[, 'x'][end_cluster], clus_center[, 'y'][end_cluster])
    clulines = rbind(clulines, temp_points)
  }
  colnames(clulines) = c("x","y", "xend","yend")
  clulines = as.data.frame(clulines)
  clulines$x.mid <-(clulines$x+clulines$xend)/2
  clulines$y.mid <-(clulines$y+clulines$yend)/2


  figure <- ggplot2::ggplot(data, ggplot2::aes(x = V1, y = V2, color = cell_stages)) +
    ggplot2::geom_point()+
    # Plot the centroid
    ggplot2::geom_point(data = as.data.frame(clus_center), ggplot2::aes(x= x,y= y), color = 'black', size = 2)+
    ggplot2::labs(x = paste0("UMAP1"), y = paste0("UMAP2"),title = paste("R =", r)) +
    # Add Sting to the figure
    ggplot2::geom_segment(ggplot2::aes_string(x = "x.mid", xend = "xend", y = "y.mid", yend = "yend", size = NULL)
                 , data = clulines, color="black", size=1.25)+
    # Add arrow to string
    ggplot2::geom_segment(arrow = ggplot2::arrow(length=ggplot2::unit(0.3,"cm"), type = "closed", ends = 'last'),
                          ggplot2::aes_string(x = "x", xend = "x.mid", y = "y", yend = "y.mid", size = NULL)
                 , data = clulines, color="black", size=1.25) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = cols) +
    #annotate(geom = "table", x=min(tsne_original$V1),y=min(tsne_original$V2), label = silhou)+
    # scale_color_npg()+
    ggsci::scale_fill_npg() +
    ggplot2::theme(legend.position = "right",
          plot.title = ggplot2::element_text(size=20),
          axis.text= ggplot2::element_text(size=20),
          axis.title=ggplot2::element_text(size=20),
          legend.key.size =ggplot2::unit(1, 'cm'), legend.text=ggplot2::element_text(size=20), legend.title=ggplot2::element_text(size=20))
  return(figure)
}

Draw_dev_pseudotime = function(dataset, label, cell.stages, pseudotime, r, cols){
  group = label
  raw <- as.data.frame(group)
  raw$time = pseudotime
  raw$group <- factor(raw$group , levels = cell.stages)

  j1 <- ggplot2::ggplot(raw, ggplot2::aes(y = group, x = time, color = group)) +
    #geom_boxplot() +
    ggplot2::geom_jitter() +
    ggplot2::labs(x = "Pseudo Time", y = '', title =paste0(' R=', r)) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
          #panel.background = element_blank(),
          plot.title = ggplot2::element_text(size=20),
          axis.text= ggplot2::element_text(size=20),
          axis.title=ggplot2::element_text(size=20),
          legend.key.size =ggplot2::unit(1, 'cm'), legend.text=ggplot2::element_text(size=20), legend.title=ggplot2::element_text(size=20),
    )
  # axis.ticks = element_blank()

  return(j1)
}

scTEP_plot = function(dataset_id, data, scDHA_res, out, fig_type, dimension_method){
  cols <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF",
            "#EE4C97FF","#000075", "#a9a9a9", "#8DD3C7", "#C8EABC", "#FBFBB4", "#D9D7C9", "#C3B4D0" ,
            "#E39699", "#E9877F", "#A9A0B2", "#97B1BD", "#D9B382",
            "#EBBD63", "#C4D367", "#C7D98C", "#EED0CD", "#F0D1E1",
            "#DED7DA", "#CDB7CE", "#BE88BF", "#C2ADC0", "#CBE5C4",
            "#E4EB9C", "#FFED6F", '#CCFF99', '#33FF00', '#FFDB6D', '#33CC33', '#003300', '#00CC99',
            '#FFFF00', '#CC9900', '#FFFFCC', '#CCFFCC',
            "#2059BB", "#16489E", "#0F3980", "#0E2F68", "#8A3ABF",
            "#7330A0", "#5A3870", '#DE1A64', '#C21A59', '#6D1234', '#3D3135', '#2ABAA4', '#5C9F95',
            '#335650', '#D55C31', '#B07966', '#E2D8D4')

  label = data$label
  cell.stages =data$stages
  startStage = as.character(data$startStage)
  cluster = out$cluster
  edge_nodes = out$milestone_network[,1:2] %>% as.data.frame()
  pseudotime = out$pseudotime
  r <- pseudotime %>% cor(as.numeric(factor(label, levels = cell.stages))) %>% round(digits = 2)

  latent = scDHA_res$latent
  if (dimension_method == "tsne"){
    latent = Rtsne::Rtsne(latent, dims = 2, initial_dims = 100, max_iter = 5000, check_duplicates = FALSE)
    latent = latent[["Y"]]
  }
  if (dimension_method == "umap"){
    latent = uwot::umap(latent)
  }

  if (fig_type == 'landscape_gt'){
    clus.stages = cluster %>% unique() %>% factor()
    Proposed_label_fig = Draw_figure_lable(dataset_id, latent, label, cell.stages, r, cols)
    plot(Proposed_label_fig)
  }

  if (fig_type == 'landscape_scTEP'){
    Proposed_clus_fig = Draw_figure_clus(dataset_id, latent, cluster, clus.stages, edge_nodes, r, cols)
    plot(Proposed_clus_fig)
  }

  if (fig_type == 'pseudotime_scTEP'){
    Proposed_pseudotime_fig = Draw_figure_pseudotime(dataset_id, latent, pseudotime, r, cols)
    plot(Proposed_pseudotime_fig)
  }

  if (fig_type == 'dev_pseudotime'){
    proposed_dev_pseudotime_fig = Draw_dev_pseudotime(dataset_id, label, cell.stages, pseudotime, r, cols)
    plot(proposed_dev_pseudotime_fig)
  }
}





