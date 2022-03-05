
n_neighbors <- 15
min_dist <- 0.01

library(dendextend)
library(caret)


## TODO: simulation data currently loaded in genefishing.R

## 
HClust_with_start_clusters <- 
  function(clusters, 
           subg_inds,
           dissim_mat, 
           cell_ids = NULL,
           hclust_method = "ward.D2") {
  
  n_subgs <- max(as.numeric(subg_inds))
  ## set distance to zero for cells within the same cluster.
  for(i in seq_len(length(unique(clusters)))) {
    cell_inds <- which(clusters==i)
    
    for(a in cell_inds) {
      for(b in cell_inds) {
        dissim_mat[a, b] = 0
      }
    }
  }
  
  ## hierarchical clustering of preprocessed matrix
  hclust_res <- stats::hclust(d = as.dist(dissim_mat), 
                              method = hclust_method)
  
  
  cut_res_list <- vector(mode = "list")
  cut_res_list[[1]] <- cutree(hclust_res, k = length(unique(clusters)))
  names(cut_res_list) <- "c" %p% length(unique(clusters))
  
  cut_res_list[[2]] <- cutree(hclust_res, k = n_subgs)  
  names(cut_res_list)[2] <- "c" %p% n_subgs
  

  if(is.null(cell_ids)) {
    cell_ids <- seq_len(ncol(dissim_mat))
  }
  
  clust_res_df <- 
    data.frame(cell_index = cell_ids,
               cluster_low  = cut_res_list[[1]],
               cluster_high = cut_res_list[[2]],
               cell_group_low  = clusters,
               cell_group_high = subg_inds)
  
  return(list(hclust_res = hclust_res,
              cut_res_list = cut_res_list,
              clust_res_df = clust_res_df))
  
}


PlotHClustResNew <- function(hclust_res_df,
                             dist_mat_for_p,
                             size = 0.5,
                             n_neighbors = 15,
                             min_dist = 0.01) {
  
  ## get HClust results
  
  ## Get Umap Embedding
  umap_res <-
    uwot::umap(as.dist(dist_mat_for_p),
               n_neighbors = n_neighbors,
               min_dist = min_dist) %>% data.frame()
  
  
  if(length(unique(hclust_res_df$cell_group_low)) == nrow(hclust_res_df)) {
    p <- 
      ggarrange( 
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cell_group_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("Actual Subgroups") +
          theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cluster_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("HClust Subgroups") +
          theme(title = element_text(size=8)),
        ncol = 2
      )
  } else{
    p <- 
      ggarrange( 
        # umap_res %>% 
        #   mutate(celltype = factor(hclust_res_df$cell_group_low)) %>%
        #   ggplot(aes(x = X1, y = X2, col = celltype)) + 
        #   geom_point(size = size) + 
        #   theme_bw() +
        #   labs(x="UMAP1", y="UMAP2") +
        #   theme(legend.position = "none") +
        #   ggtitle("Actual Low-Level Groups") +
        #   theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cell_group_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("Actual High-Level Groups") +
          theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cluster_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("HClust Subgroups") +
          theme(title = element_text(size=8)),
        ncol = 2
      )
  }
  
  return(p)
  
}




RunHClustPipeline <- 
  function(clusters, 
           subg_inds, 
           hclust_dissim_mat, 
           p_dissim_mat,
           cell_ids = NULL) {
    
  hclust_res <- 
    HClust_with_start_clusters(clusters   = clusters, 
                               subg_inds  = subg_inds,
                               dissim_mat = hclust_dissim_mat,
                               cell_ids = cell_ids)
  
  
  p <- PlotHClustResNew(hclust_res_df = hclust_res$clust_res_df,
                        dist_mat_for_p = p_dissim_mat)
  
  
  ## Hierarchical Rand Index
  
  adj_rand_index_high = 
    mclust::adjustedRandIndex(hclust_res$clust_res_df$cell_group_high, 
                              hclust_res$clust_res_df$cluster_high)
  
  adj_rand_index_low = 
    mclust::adjustedRandIndex(hclust_res$clust_res_df$cell_group_low, 
                              hclust_res$clust_res_df$cluster_low)
  
  return(list(hclust_res = hclust_res, 
              p = p, 
              adj_rand_index_high = adj_rand_index_high,
              adj_rand_index_low = adj_rand_index_low))
}
