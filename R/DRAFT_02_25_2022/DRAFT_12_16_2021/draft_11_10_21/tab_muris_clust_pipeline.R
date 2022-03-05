TabMurisProcSubset <- 
  function(full_data, meta_data, filename=NULL,
           var_features) {
    if(!is.null(filename)) {
      cell_inds <- which(meta_data$file_source == filename)
    } else{
      cell_inds = seq_len(ncol(full_data))
    }
    
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
  function(full_data, meta_data, filename=NULL, 
           spectral_dims = 10, 
           n_clust = 20,
           n_cores = 1,
           R = 150,
           G = c(50, 150, 300),
           PCs = c(3, 25),
           var_features = 3000) {
    
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
    for(dist_name in names(dist_res$dist_list[
      !grepl("pca", names(dist_res$dist_list))])) {
      n <- ncol(dist_res$dist_list[[1]])
      
      l <- length(dist_res$dist_list)
      
      dissim <- dist_res$dist_list[[dist_name]]
      
      spec_embed_res <- spectralEmbed(dissim, ndim = spectral_dims) 
      
      ## standard:
      dissim <- dist(spec_embed_res$U)
      
      mat <- matrix(0, nrow = n, ncol = n)
      
      mat[lower.tri(mat)]  <- c(dissim)
      
      mat[upper.tri(mat)] <-
        t(mat)[upper.tri(mat)]
      
      dist_res$dist_list[[l+1]] <- mat
      names(dist_res$dist_list)[l+1] <- 
        dist_name %p% "_spectral_d" %p% spectral_dims %p% "_dist"
      
      dist_res$umap_list[[l+1]] <- 
        uwot::umap(as.dist(dist_res$dist_list[[l+1]]),
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
      
      names(dist_res$umap_list)[l+1] <- names(dist_res$dist_list)[l+1]
      ## weighted:
      U_wtd = spec_embed_res$U
      for(c in seq_len(ncol(spec_embed_res$U))) {
        U_wtd[,c] <- U_wtd[,c] * spec_embed_res$eig_weights[c]
      }
      
      dissim <- dist(U_wtd)
      
      mat <- matrix(0, nrow = n, ncol = n)
      
      mat[lower.tri(mat)]  <- c(dissim)
      
      mat[upper.tri(mat)] <-
        t(mat)[upper.tri(mat)]
      
      dist_res$dist_list[[l+2]] <- mat
      names(dist_res$dist_list)[l+2] <- 
        dist_name %p% "_spectral_wtd_d" %p% spectral_dims %p% "_dist"
      
      dist_res$umap_list[[l+2]] <- 
        uwot::umap(as.dist(dist_res$dist_list[[l+2]]),
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
      
      names(dist_res$umap_list)[l+2] <- names(dist_res$dist_list)[l+2]
      

    }
    
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
