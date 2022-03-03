source(here("R/SIMLR_no_clust.R"))

computeAllDists <- function(counts, group_labels, 
                            r = 100, g = 300, 
                            pcs = c(25),
                            simlr_dims = c(10),
                            elbow_pca = TRUE,
                            n_pcs = NULL, n_clust = NULL, 
                            pca_elbow_chg_thres = 0.025,
                            kmeans_elbow_chg_thres = 0.025,
                            rafsil_nrep = 20,
                            n_cores = detectCores()-2,
                            save_loc = NULL,
                            ## umap hypers
                            n_neighbors = 15, 
                            min_dist = 0.01) {
  
  ### 1. SIDEREF --------------------------------------------------------------
  side_ref = list()
  side_ref_dist = list()
  side_ref_umap = list()
  
  for(i in seq_len(length(g))) {
    cat("sideref rand" %p% i %p% "...")
    
    side_ref[[i]] <- 
      SIDEREF(expr_matrix = counts, 
              n_top_genes = g[i],
              R = r,
              selection_method = "cell_embed_sample",
              diff_expr_method = "diff_expr_norm",
              similarity_method = "n_intersect",
              n_pcs = n_pcs,
              n_clust = n_clust,
              pca_elbow_chg_thres = pca_elbow_chg_thres,
              kmeans_elbow_chg_thres = kmeans_elbow_chg_thres,
              ## parallelization parameters
              parallelize = TRUE,
              n_cores = n_cores,
              ## other
              verbose = TRUE)
    
    side_ref_dist[[i]] <- side_ref[[i]]
    
    side_ref_umap[[i]] <-
      uwot::umap(as.dist(side_ref_dist[[i]]),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist)

  }
  
  ### 2. PCA ------------------------------------------------------------------
  pca_dist  = list()
  pca_umap  = list()
  
  
  pca_wtd_dist  = list()
  pca_wtd_umap  = list()
  
  
  pca_res_full <- prcomp(t(counts))
  
  ## include quick elbow PCA result if specified.
  if(elbow_pca) {
    elb_pcs <- quick_elbow(varpc = pca_res_full$sdev**2,
                           low = pca_elbow_chg_thres,
                           max_expl = 1)
    
    if(!elb_pcs %in% pcs) {
      pcs <- append(pcs, elb_pcs)  
    }
  }
  
  for(i in seq_len(length(pcs))) {
    cat(pcs[i] %p% " PCs...")
    
    ## principal component stdevs
    pca_sdevs <- pca_res_full$sdev[1:(min(dim(counts)[2], pcs[i]))]
    
    ## Embedding
    pca_res <- pca_res_full$x[, 1:(min(dim(counts)[2], 
                                       pcs[i]))]
    
    pca_wtd_res <- pca_res
    for(c in seq_len(pcs[i])) {
      pca_wtd_res[,c] <- pca_res[,c] * pca_sdevs[c]
    }
    
    pca_res     <- dist(pca_res)
    pca_wtd_res <- dist(pca_wtd_res)
    
    
    mat <- matrix(0, nrow = ncol(counts),
                  ncol = ncol(counts))
    
    ## PCA (regular)
    mat[lower.tri(mat)]  <- c(pca_res)
    
    mat[upper.tri(mat)] <- 
      t(mat)[upper.tri(mat)]
    
    pca_dist[[i]] <- mat
    
    ## PCA (weighted)
    mat[lower.tri(mat)]  <- c(pca_wtd_res)
    
    mat[upper.tri(mat)] <- 
      t(mat)[upper.tri(mat)]
    
    pca_wtd_dist[[i]] <- mat
    
    
    
    pca_umap[[i]] <-
      uwot::umap(as.dist(pca_dist[[i]]),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist) 
    
    pca_wtd_umap[[i]] <-
      uwot::umap(as.dist(pca_wtd_dist[[i]]),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist) 
    
    
    
  }
  
  
  ### 3. Euclid UMAP ------------------------------------------------------------
  
  cat("euclidean and correlation ...")
  
  euclid_dist <- dist(t(counts))
  
  mat <- matrix(0, nrow = ncol(counts),
                ncol = ncol(counts))
  
  mat[lower.tri(mat)]  <- c(euclid_dist)
  
  mat[upper.tri(mat)] <- 
    t(mat)[upper.tri(mat)]
  
  euclid_dist <- mat
  
  
  euclid_umap <-
    uwot::umap(as.dist(euclid_dist),
               n_neighbors = n_neighbors,
               min_dist = min_dist) 
  
  
  
  ### 4. Pearson UMAP -----------------------------------------------------------
  
  pear_dist <- 1 - abs(cor(counts))
  
  pear_umap <-
    uwot::umap(as.dist(pear_dist),
               n_neighbors = n_neighbors,
               min_dist = min_dist)
  
  
  
  ### 5. Spearman UMAP ----------------------------------------------------------
  
  spearman_dist <- 1 - abs(cor(counts, method = "spearman"))
  
  spearman_umap <-
    uwot::umap(as.dist(spearman_dist),
               n_neighbors = n_neighbors,
               min_dist = min_dist)
  
  
  
  names_core <- 
    c("side_ref_g" %p% g,
      "pca_" %p% pcs,
      "pca_wtd" %p% pcs,
      "euclid", "pear", "spearman")
  
  names_dist <- names_core %p% "_dist"
  names_umap <- names_core %p% "_umap"
  
  output <-
    list(
      dist_list = append(side_ref_dist, 
                         append(pca_dist,
                                append(pca_wtd_dist,
                                       list(euclid_dist,
                                            pear_dist,
                                            spearman_dist)))),
      umap_list = append(side_ref_umap,
                         append(pca_umap,
                                append(pca_wtd_umap,
                                       list(euclid_umap,
                                            pear_umap,
                                            spearman_umap))))
    )
  names(output$dist_list) = names_dist 
  names(output$umap_list) = names_umap 
  
  
  
  ## SAVE of main distance measures, if specified.
  if(!is.null(save_loc)) {
    save(output, file = save_loc %p% "_main_dist_comps.RData")
  } 
  
  
  
  
  ### scRNA-seq specific alternatives to SIDEREF:
  
  othr_dist_output <- vector(mode = "list", length = 1 + length(simlr_dims))
  
  ### 1: RAFSIL
  othr_dist_output[[1]] <-
    RAFSIL(t(counts), NumC=NULL, nrep=rafsil_nrep, 
           method="RAFSIL1", gene_filter = FALSE,
           parallelize=TRUE, n_cores = n_cores)$D
  
  
  
  ### 2: SIMLR
  for(i in seq_len(length(simlr_dims))) {
    simlr_res <- SIMLR_no_clust(counts, c = simlr_dims[i])
    ## convert to distance
    dist_simlr <- max(simlr_res$S) - simlr_res$S
    diag(dist_simlr) <- 0
    
    ## store to list
    othr_dist_output[[1+i]] <- 
      dist_simlr
    
  }
  
  names(othr_dist_output) <- 
    c("RAFSIL", "SIMLR_" %p% simlr_dims %p% "_dims")
  
  
  if(!is.null(save_loc)) {
    save(othr_dist_output, file = save_loc %p% "_othr_scrna_seq_dist_comps.RData")
  } 
  
  return(list(output=output, othr_dist_output=othr_dist_output))
  
}


