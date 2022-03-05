gg_color_ggplot <- function(n){
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

HGC.PlotDendrogram <- function(tree, k = 5, plot.label = FALSE, labels){
  checkDataStructure(tree, "tree", "hclust")
  ncell <- length(tree$order)
  if(ncell == 0){
    stop("There is no element in the tree.")
  }
  checkInteger(k, "k")
  if(k > ncell){
    stop("k should be smaller than cell number.")
  }
  if(plot.label){
    if(nrow(labels) != ncell){
      stop("The labels do not match the data.")
    }
  }
  clus.dendrogram <- stats::as.dendrogram(tree)
  dend <- clus.dendrogram %>%
    dendextend::set("branches_k_color", k=k) %>%
    dendextend::set("branches_lwd", 1.2) %>%
    dendextend::set("labels_colors", "white")
  if(plot.label == FALSE){
    plot(dend)
    return(1)
  }else{
    n_label <- ncol(labels)
    n_diff_labels <- rep(0, n_label)
    for(i in seq_len(n_label)){
      temp_n_label <- length(unique(labels[,i]))
      n_diff_labels[i] <- temp_n_label
    }
    T.colourbar <- matrix(nrow = ncell, ncol = n_label)
    colnames(T.colourbar) = colnames(labels)
    for(i in seq_len(n_label)){
      temp_label <- labels[,i]
      temp_u_label <- unique(labels[,i])
      temp_n_label <- n_diff_labels[i]
      temp_colourbar <- gg_color_ggplot(temp_n_label)
      for(j in seq_len(temp_n_label)){
        T.colourbar[which(temp_label ==
                            temp_u_label[j]),i] = temp_colourbar[j]
      }
    }
    plot(dend)
    dendextend::colored_bars(colors = T.colourbar[
      order.dendrogram(dend),])
    return(1)
  }
}