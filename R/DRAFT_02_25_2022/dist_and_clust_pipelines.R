source(here("R/spectral_dist_funcs.R"))

selectFullDistIdx <- function(cell_id_vec,
                              dist,
                              meta_data) {
  
  meta_data <- 
    meta_data %>% mutate(idx = seq_len(nrow(meta_data)))
  
  idx_subset <- meta_data[cell_vec, ] %>% pull(idx)
  
  return(dist[idx_subset, idx_subset])
}


getGroupLabels <- function(cell_id_vec, meta_data) {
  if(!("group_label" %in% names(meta_data))) {
    stop("group_label variable is not in DF meta_data.")}
  group_labels <- 
    as.factor(as.character(meta_data[cell_id_vec, ]$group_label))
  
  return(group_labels)
}



TabMurisProcSubset <- 
  function(full_data, meta_data, filename=NULL,
           var_features,
           return_seurat = FALSE) {
    if(!is.null(filename)) {
      cell_inds <- which(meta_data$file_source == filename)
    } else{ 
      cell_inds <- seq_len(ncol(full_data))
    }
    
    full_data <- full_data[, cell_inds]
    
    full_data <-
      CreateSeuratObject(full_data, project = "SeuratProject", 
                         assay = "RNA")
    
    full_data <- FindVariableFeatures(full_data, nfeatures = var_features)
    full_data <- NormalizeData(full_data)
    full_data <- ScaleData(full_data)
    
    full_data_subset_variable <- full_data@assays$RNA@scale.data
    
    if(return_seurat) {
      return(list(full_data, cell_inds))
    } else{
      return(list(full_data_subset_variable, cell_inds))
    }
    
  }


splatDistResPipeline <- 
  function(full_data, group_labels, 
           spectral_dims = 10, 
           n_clust = 20,
           n_cores = 1,
           R = 100,
           G = c(300),
           PCs = c(25),
           var_features = 3000,
           do_spectral_pca = TRUE) {
    
    ## reset n_clust for small data
    n_clust = min(floor(ncol(full_data) / 3), n_clust)
    
    ### Preproc 
    var_feature_df <- FindVariableFeatures(full_data, nfeatures = var_features)
    var_feature_df <- var_feature_df %>% 
      mutate(gene_idx = seq_len(nrow(var_feature_df))) %>% 
      arrange(-vst.variance.standardized) %>% pull(gene_idx)
    
    var_feature_df <- var_feature_df[1:var_features]
    
    full_data <- ScaleData(full_data[sort(var_feature_df), ])
    
    print("Dim Data: " %p% dim(full_data)[1] %p% ", " %p% dim(full_data)[2])
    
    ## Run dim Reds
    dist_res <-
      runDimReds(counts = full_data, 
                 group_labels = group_labels, 
                 r = R, g = G, pcs = PCs,
                 n_clust = n_clust, n_cores = n_cores)
    
    ###  Appending Spectral Results 
    print("Spectral distance computation")
    dist_res <- addSpectralDistRes(dist_res,
                                   spectral_dims = spectral_dims,
                                   do_spectral_pca = do_spectral_pca)
    
    return(dist_res)
    
  }





TabMurisDistResPipeline <- 
  function(full_data, meta_data, filename = NULL, 
           spectral_dims = 10, 
           n_clust = 20,
           n_cores = 1,
           R = 150,
           G = c(300),
           PCs = c(25),
           var_features = 3000, 
           do_spectral_pca = TRUE) {
    
    sub_res <- TabMurisProcSubset(full_data, meta_data, filename,
                                  var_features = var_features) 
    
    full_data_subset_variable <- sub_res[[1]]
    cell_inds <- sub_res[[2]]
    
    ## reset n_clust for small data
    n_clust = min(floor(ncol(full_data_subset_variable) / 3), n_clust)
    
    
    ## Run dim Reds
    dist_res <-
      runDimReds(counts = full_data_subset_variable, 
                 group_labels = meta_data$cell_type[cell_inds], 
                 r = R, g = G, pcs = PCs,
                 n_clust = n_clust, n_cores = n_cores)
    
    ###  Appending Spectral Results 
    print("Spectral distance computation")
    dist_res <- addSpectralDistRes(dist_res, 
                                   spectral_dims = spectral_dims,
                                   do_spectral_pca = do_spectral_pca)
    
    return(dist_res)
    
  }


TabMurisClustResPipeline <- 
  function(full_data, meta_data, filename,
           spectral_dims = 10, 
           n_cores = 1, 
           n_clust = 20,
           R = 150, G = c(50, 150, 300), 
           PCs = c(3, 25),
           var_features = 3000) {
  
    sub_res <- TabMurisProcSubset(full_data, meta_data, filename,
                                  var_features = var_features) 
    
    full_data_subset_variable <- sub_res[[1]]
    cell_inds <- sub_res[[2]]
    cell_names <- colnames(full_data_subset_variable)
    data_full_subset <- full_data[, cell_inds]
    
    
    dist_res <- 
      TabMurisDistResPipeline(full_data, meta_data, filename, 
                              spectral_dims = spectral_dims, 
                              n_cores = n_cores,
                              n_clust = n_clust,
                              R = R,
                              G = G,
                              PCs = PCs,
                              var_features = var_features)
  
    ### Low-level clustering pipeline:
    data_seurat <-
      CreateSeuratObject(data_full_subset, project = "SeuratProject", 
                         assay = "RNA")
    
    data_seurat@assays$RNA@scale.data <- 
      full_data_subset_variable
    
    clusts <- list()
    modularities <- list()
    avg_sils <- list()
    
    for(x in names(dist_res$dist_list)) {
      print(x)
      dist = dist_res$dist_list[[x]]
      
      rownames(dist) = colnames(data_seurat)
      colnames(dist) = colnames(data_seurat)
      
      
      nn_graphs <- FindNeighbors(dist, distance.matrix=TRUE)
      
      data_seurat@graphs$RNA_nn <- nn_graphs$nn
      data_seurat@graphs$RNA_snn <- nn_graphs$snn
      
      data_seurat <- FindClusters(data_seurat)
      
      diag(nn_graphs$nn) = 0
      diag(nn_graphs$snn) = 0
      
      m_nn = 0.5*sum(nn_graphs$nn)
      m_snn = 0.5*sum(nn_graphs$snn)
      
      k_nn = rowSums(nn_graphs$nn)
      k_snn = rowSums(nn_graphs$snn)
      
      Q_nn = 0
      Q_snn = 0
      for(i in 1:ncol(nn_graphs$nn)) {
        for(j in 1:ncol(nn_graphs$nn)) {
          if(data_seurat@meta.data$seurat_clusters[i] == 
             data_seurat@meta.data$seurat_clusters[j] & i != j) {
            Q_nn = Q_nn + 1/(2*m_nn) * (nn_graphs$nn[i,j] - k_nn[i]*k_nn[j]/(2*m_nn))
            Q_snn = Q_snn + 1/(2*m_snn) * (nn_graphs$snn[i,j] - k_snn[i]*k_snn[j]/(2*m_snn))
          }
        }
      }
      
      modularities[[x]] <- c(Q_nn, Q_snn)
      
      clusts[[x]] <- data_seurat@meta.data$seurat_clusters
      names(clusts[[x]]) <- cell_names
      
      ## compute average sil width
      if(length(unique(clusts[[x]])) == 1) {
        avg_sils[[x]] = NA
      } else{
        sil_scores <- 
          cluster::silhouette(as.integer(factor(clusts[[x]])),
                              dist = stats::dist(dist, 
                                                 method = "euclidean"))[,3]
        avg_sils[[x]] <- mean(sil_scores)
      }
      
    }
    
    return(list(clusts=clusts, modularities = modularities, avg_sils = avg_sils,
                umaps=dist_res$umap_list))
    
  }


splatUMAPPlot <- function(data, title, celltypes, color_mapping = NULL,
                          pcs = NULL, plot_sil_score = FALSE) {
  
  if(is.null(pcs)) {
    umap_res <-
      uwot::umap(data,
                 n_neighbors = n_neighbors,
                 min_dist = min_dist) 
  } else{
    umap_res <-
      uwot::umap(data,
                 n_neighbors = n_neighbors,
                 min_dist = min_dist,
                 pca = pcs) 
  }
  
  
  
  ## Group Avg Sil Scores:
  if(plot_sil_score) {
    sil_scores <- cluster::silhouette(as.integer(factor(celltypes)),
                                      dist = stats::dist(umap_res, 
                                                         method = "euclidean"))[,3]
    avg_sils <- sapply(1:length(unique(celltypes)), function(i){mean(sil_scores[as.integer(factor(celltypes)) == i])})
    
    
    levels(celltypes) <- paste0(levels(celltypes), " (", round(avg_sils, 2), ")") 
  }
  
  p <-
    data.frame(umap_res) %>% 
    mutate(celltype = celltypes) %>%
    ggplot(aes(x = X1, y = X2, col = celltype)) + 
    geom_point() + 
    labs(x="UMAP1", y="UMAP2", 
         col = ifelse(plot_sil_score, 
                      "Cell Type (Sil. Score)",
                      "Cell Type")) + 
    theme_bw() + 
    ggtitle(title)
  
  if(!is.null(color_mapping)) {
    cols = brewer.pal(max(color_mapping), "Set1")
    cols = cols[color_mapping]
    p <- p + 
      scale_colour_manual(values=cols)
  }
  
  return(p)
}

