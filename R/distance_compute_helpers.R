source(here("R/computeAllDists.R"))
source(here("R/spectral_dist_funcs.R"))


selectFullDistIdx <- function(cell_id_vec,
                              dist,
                              meta_data) {
  
  if(is.null(cell_id_vec)) {return(dist)
  } else{
    meta_data <- 
      meta_data %>% mutate(idx = seq_len(nrow(meta_data)))
    
    idx_subset <- meta_data[cell_id_vec, ] %>% pull(idx)
    
    return(dist[idx_subset, idx_subset])    
  }
  
}


getGroupLabels <- function(meta_data, cell_id_vec = NULL) {
  if(!("group_label" %in% names(meta_data))) {
    stop("group_label variable is not in DF meta_data.")}
  
  if(!is.null(cell_id_vec)) {
    group_labels <- 
      as.factor(as.character(meta_data[cell_id_vec, ]$group_label))  
  } else{
    group_labels <- as.factor(as.character(meta_data$group_label))
  }
  
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



TabMurisDistResPipeline <- 
  function(full_data, meta_data, filename = NULL, 
           spectral_dims = 10, 
           n_cores = 1,
           R,
           G,
           PCs = c(25),
           var_features = 3000, 
           save_loc = NULL,
           ## additional Niche params:
           do_spectral_pca = TRUE,
           elbow_pca = TRUE,
           ## SIDEREF hypers:
           n_pcs = NULL,
           n_clust = NULL,
           pca_elbow_chg_thres = PC_ELBOW_CHG_THRES,
           kmeans_elbow_chg_thres = KM_ELBOW_CHG_THRES,
           rafsil_nrep = RAFSIL_NREP,
           n_neighbors = 15,
           min_dist = 0.01) {
    
    sub_res <- TabMurisProcSubset(full_data, meta_data, filename,
                                  var_features = var_features)
    
    full_data_subset_variable <- sub_res[[1]]
    cell_inds <- sub_res[[2]]
    
    
    ## Run distance computes
    dist_res <-
      computeAllDists(counts = full_data_subset_variable,
                      group_labels = meta_data$cell_type[cell_inds],
                      r = R, g = G,
                      pcs = PCs,
                      simlr_dims = spectral_dims,
                      elbow_pca = elbow_pca,
                      n_pcs = n_pcs, n_clust = n_clust,
                      pca_elbow_chg_thres = pca_elbow_chg_thres,
                      kmeans_elbow_chg_thres = kmeans_elbow_chg_thres,
                      rafsil_nrep = rafsil_nrep,
                      n_cores = n_cores,
                      save_loc = save_loc,
                      n_neighbors = n_neighbors,
                      min_dist = min_dist)
    
    
    ###  Appending Spectral Results 
    
    print("Spectral distance computation")
    ## We run on SIDEREF and PCA (weighted)
    set.seed(1)
    
    spectral_dist_res <- dist_res$main_dist_output
    
    spectral_dist_names <-
      names(dist_res$main_dist_output$dist_list)[
        grepl("pca_wtd|side", names(dist_res$main_dist_output$dist_list))]
    
    spectral_umap_names <-
      names(dist_res$main_dist_output$umap_list)[
        grepl("pca_wtd|side", names(dist_res$main_dist_output$umap_list))]
    
    spectral_dist_res$dist_list <- 
      spectral_dist_res$dist_list[spectral_dist_names]
    spectral_dist_res$umap_list <- 
      spectral_dist_res$umap_list[spectral_umap_names]
    
    spectral_dist_res <- addSpectralDistRes(
      spectral_dist_res,
      spectral_dims   = spectral_dims,
      do_spectral_pca = do_spectral_pca,
      save_loc = save_loc)
    
    
    return(list(main_dist_res = dist_res$main_dist_output, 
                other_dist_res = dist_res$othr_dist_output,
                spectral_dist_res = spectral_dist_res))
    
  }


splatDistResPipeline <- 
  function(full_data, group_labels, 
           spectral_dims = 10, 
           n_cores = 1,
           R = 150,
           G = c(300),
           PCs = c(25),
           var_features = 3000, 
           save_loc = NULL,
           ## additional Niche params:
           do_spectral_pca = TRUE,
           elbow_pca = TRUE,
           ## SIDEREF hypers:
           n_pcs = NULL,
           n_clust = NULL,
           pca_elbow_chg_thres = PC_ELBOW_CHG_THRES,
           kmeans_elbow_chg_thres = KM_ELBOW_CHG_THRES,
           rafsil_nrep = RAFSIL_NREP,
           n_neighbors = 15, 
           min_dist = 0.01) {
    
    ### Preproc 
    var_feature_df <- FindVariableFeatures(full_data, nfeatures = var_features)
    var_feature_df <- var_feature_df %>% 
      mutate(gene_idx = seq_len(nrow(var_feature_df))) %>% 
      arrange(-vst.variance.standardized) %>% pull(gene_idx)
    
    var_feature_df <- var_feature_df[1:var_features]
    
    full_data <- ScaleData(full_data[sort(var_feature_df), ])
    
    print("Dim Data: " %p% dim(full_data)[1] %p% ", " %p% dim(full_data)[2])
    
    ## Run distance computes.
    dist_res <-
      computeAllDists(counts = full_data, 
                      group_labels = group_labels, 
                      r = R, g = G, 
                      pcs = PCs,
                      simlr_dims = spectral_dims,
                      elbow_pca = elbow_pca,
                      n_pcs = n_pcs, n_clust = n_clust, 
                      pca_elbow_chg_thres = pca_elbow_chg_thres,
                      kmeans_elbow_chg_thres = kmeans_elbow_chg_thres,
                      rafsil_nrep = rafsil_nrep,
                      n_cores = n_cores,
                      save_loc = save_loc,
                      n_neighbors = n_neighbors, 
                      min_dist = min_dist)
    
    ###  Appending Spectral Results 
    print("Spectral distance computation")
    ## We run on SIDEREF and PCA (weighted)
    spectral_dist_res <- dist_res$main_dist_output
    
    spectral_dist_names <-
      names(dist_res$main_dist_output$dist_list)[
        grepl("pca_wtd|side", names(dist_res$main_dist_output$dist_list))]
    
    spectral_umap_names <-
      names(dist_res$main_dist_output$umap_list)[
        grepl("pca_wtd|side", names(dist_res$main_dist_output$umap_list))]
    
    spectral_dist_res$dist_list <- 
      spectral_dist_res$dist_list[spectral_dist_names]
    spectral_dist_res$umap_list <- 
      spectral_dist_res$umap_list[spectral_umap_names]
    
    spectral_dist_res <- addSpectralDistRes(
      spectral_dist_res,
      spectral_dims   = spectral_dims,
      do_spectral_pca = do_spectral_pca,
      save_loc = save_loc)
    
    return(list(main_dist_res = dist_res$main_dist_output, 
                other_dist_res = dist_res$othr_dist_output,
                spectral_dist_res = spectral_dist_res))
    
  }




getKNNList <- function(dist_mat, k) {
  knn_list <- 
    lapply(seq_len(dim(dist_mat)[1]),
           function(i) {
             return(order(dist_mat[,i])[-1][1:k])
           })
}
