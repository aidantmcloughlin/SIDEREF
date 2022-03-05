TabMurisProcSubset <- 
  function(full_data, meta_data, filename) {
    cell_inds <- which(meta_data$file_source == filename)
    full_data_subset <- full_data[, cell_inds]
    
    full_data_subset <-
      CreateSeuratObject(full_data_subset, project = "SeuratProject", 
                         assay = "RNA")
    
    full_data_subset <- FindVariableFeatures(full_data_subset, nfeatures = var_features)
    full_data_subset <- NormalizeData(full_data_subset)
    full_data_subset <- ScaleData(full_data_subset)
    
    full_data_subset_variable <- full_data_subset@assays$RNA@scale.data
    
    return(list(full_data_subset_variable, cell_inds))
  }

TabMurisDistResPipeline <- 
  function(full_data, meta_data, filename, 
           spectral_dims = 10, 
           cores = 1,
           R = 150,
           G = c(50, 150, 300),
           var_features = 3000) {
    
    sub_res <- TabMurisProcSubset(full_data, meta_data, filename) 
    
    full_data_subset_variable <- sub_res[[1]]
    cell_inds <- sub_res[[2]]
    
    ## Run dim Reds
    dist_res_tab_muris_0.5pct_3000g <-
      runDimReds(full_data_subset_variable, 
                 group_labels = meta_data$cell_type[cell_inds], 
                 r = R, g=G, 
                 n_clust = 20, n_cores = cores)
    
    ###  Appending Spectral Results
    for(dist_name in names(dist_res_tab_muris_0.05pct_3000g$dist_list[
      !grepl("pca", names(dist_res_tab_muris_0.05pct_3000g$dist_list))])) {
      print(dist_name)
      n <- ncol(dist_res_tab_muris_0.05pct_3000g$dist_list[[1]])
      
      l <- length(dist_res_tab_muris_0.05pct_3000g$dist_list)
      
      dissim <- dist_res_tab_muris_0.05pct_3000g$dist_list[[dist_name]]
      
      dissim <- dist(spectralEmbed(dissim, ndim = spectral_dims))
      
      mat <- matrix(0, nrow = n, ncol = n)
      
      mat[lower.tri(mat)]  <- c(dissim)
      
      mat[upper.tri(mat)] <-
        t(mat)[upper.tri(mat)]
      
      dist_res_tab_muris_0.05pct_3000g$dist_list[[l+1]] <- mat
      names(dist_res_tab_muris_0.05pct_3000g$dist_list)[l+1] <- 
        dist_name %p% "_spectral_d" %p% spectral_dims %p% "_dist"
    }
    
    return(dist_res_tab_muris_0.05pct_3000g)
    
  }


TabMurisClustResPipeline <- 
  function(full_data, meta_data, filename,
           spectral_dims = 10, cores = 1, 
           R = 150, G = c(50, 150, 300), 
           var_features = 3000) {
  
    sub_res <- TabMurisProcSubset(full_data, meta_data, filename) 
    
    data_full_variable <- sub_res[[1]]
    data_full_subset <- full_data[, cell_inds]
    cell_inds <- sub_res[[2]]
    
    dist_res <- 
      TabMurisDistResPipeline(full_data, meta_data, filename, 
                              spectral_dims = spectral_dims, 
                              cores = cores,
                              R = R,
                              G = G,
                              var_features = var_features)
  
    ### Low-level clustering pipeline:
    data_seurat <-
      CreateSeuratObject(data_full_subset, project = "SeuratProject", 
                         assay = "RNA")
    
    data_seurat@assays$RNA@scale.data <- 
      full_data_subset_variable
    
    clusts <- list()
    
    for(x in names(dist_res_tab_muris_0.05pct_3000g$dist_list)) {
      print(x)
      dist = dist_res_tab_muris_0.05pct_3000g$dist_list[[x]]
      
      rownames(dist) = colnames(data_seurat)
      colnames(dist) = colnames(data_seurat)
      
      
      nn_graphs <- FindNeighbors(dist, distance.matrix=TRUE)
      
      data_seurat@graphs$RNA_nn <- nn_graphs$nn
      data_seurat@graphs$RNA_snn <- nn_graphs$snn
      
      data_seurat <- FindClusters(data_seurat)
      
      clusts[[x]] <- data_seurat@meta.data$seurat_clusters
      
    }
    
    return(clusts)
    
  }
