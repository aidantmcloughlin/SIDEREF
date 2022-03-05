###############################################################################
### Primary Plotting Functions 
###############################################################################


### Splatter:
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






snareUMAPPlot <- function(data, title, celltypes, size = 0.9,
                          pcs = NULL,
                          plot_sil_score = FALSE) {
  
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
    geom_point(size = size) + 
    labs(x="UMAP1", y="UMAP2", 
         col = ifelse(plot_sil_score, 
                      "Cell Type (Sil. Score)",
                      "Cell Type")) + 
    theme_bw() + 
    ggtitle(title) +
    guides(colour = guide_legend(override.aes = list(size=2)))
  return(p)
}



runDimRedsAndPlot <- function(sim_counts, group_labels, r = 50,
                              return_dist_mats = FALSE) {
  ### PCA => UMAP
  
  cat("PCA plots...")
  
  umap_res_pc_75 <-
    uwot::umap(t(as.matrix(sim_counts)),
               n_neighbors = n_neighbors,
               min_dist = min_dist,
               pca = 75) 
  
  p_splat_pca_umap_75 <-
    splatUMAPPlot(umap_res_pc_75, title = "Euclidean Applied to 75 PCs", 
                  celltypes = group_labels)
  
  
  umap_res_pc_25 <-
    uwot::umap(t(as.matrix(sim_counts)),
               n_neighbors = n_neighbors,
               min_dist = min_dist,
               pca = 25) 
  
  p_splat_pca_umap_25 <-
    splatUMAPPlot(umap_res_pc_25, title = "Euclidean Applied to 25 PCs", 
                  celltypes = group_labels)
  
  
  umap_res_pc_5 <-
    uwot::umap(t(as.matrix(sim_counts)),
               n_neighbors = n_neighbors,
               min_dist = min_dist,
               pca = 5) 
  
  p_splat_pca_umap_5 <-
    splatUMAPPlot(umap_res_pc_5, title = "Euclidean Applied to 5 PCs", 
                  celltypes = group_labels)
  
  umap_res_pc_3 <-
    uwot::umap(t(as.matrix(sim_counts)),
               n_neighbors = n_neighbors,
               min_dist = min_dist,
               pca = 3) 
  
  p_splat_pca_umap_3 <-
    splatUMAPPlot(umap_res_pc_3, title = "Euclidean Applied to 3 PCs", 
                  celltypes = group_labels)
  
  
  
  
  ### 3. Euclid UMAP ------------------------------------------------------------
  
  cat("euclidean and correlation plots...")
  
  euclid_dist = dist(t(as.matrix(sim_counts)))
  
  p_splat_euclid <-
    splatUMAPPlot(t(as.matrix(sim_counts)), title = "Euclidean Distance", 
                  celltypes = group_labels)
  
  
  ### 4. Pearson UMAP -----------------------------------------------------------
  
  pear_dist <- 1 - abs(cor(sim_counts))
  
  p_splat_pear <-
    splatUMAPPlot(as.dist(pear_dist), title = "Pearson Distance", 
                  celltypes = group_labels)
  
  
  ### 5. Spearman UMAP ----------------------------------------------------------
  
  spearman_dist <- 1 - abs(cor(sim_counts, method = "spearman"))
  
  p_splat_spear <-
    splatUMAPPlot(as.dist(spearman_dist), title = "Spearman Distance", 
                  celltypes = group_labels)
  
  
  
  ### 6. SIDERef Rand -----------------------------------------------------------
  
  cat("sideref rand...")
  
  side_seq_ref_rand_50 <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = 50,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)
  
  plot_list_splat_rand_50 <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand_50$dissim_final),
                  title = "SIDEREF (random) 50 Cells 50 Top Genes",
                  celltypes = group_labels)
  
  ## Spectral Embedding of the Similarity scores
  
  U <- spectralEmbed(side_seq_ref_rand_50$dissim_final, ndim = 10)
  
  plot_list_splat_rand_50_spectral <- 
    splatUMAPPlot(data  = dist(U),
                  title = "SIDEREF (100 Genes) => Spectral(10Dims) => Euclid => UMAP",
                  celltypes = group_labels)
  
  ## Spearman
  spear_sim <- abs(cor(sim_counts, method = "spearman"))
  diag(spear_sim) = 0
  ## spectral embedding
  U <- spectralEmbed(spear_sim,
                     is_dissim = FALSE,
                     ndim = 10)
  
  plot_list_splat_spear_spectral <- 
    splatUMAPPlot(data  = dist(U),
                  title = "Spearman => Spectral(10Dims) => Euclid => UMAP",
                  celltypes = group_labels)
  
  
  side_seq_ref_rand_150 <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = 150,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)
  
  plot_list_splat_rand_150 <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand_150$dissim_final),
                  title = "SIDEREF (random) 50 Cells 150 Top Genes",
                  celltypes = group_labels)
  
  side_seq_ref_rand_300 <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = 300,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)
  
  plot_list_splat_rand_300 <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand_300$dissim_final),
                  title = "SIDEREF (random) 50 Cells 300 Top Genes",
                  celltypes = group_labels)
  if(return_dist_mats == FALSE) {
    return(plot_list = list(
      plot_list_splat_rand, plot_list_splat_rand_150, plot_list_splat_rand_300,
      p_splat_pca_umap_3, p_splat_pca_umap_25, p_splat_pca_umap_75,
      p_splat_euclid, p_splat_pear, p_splat_spear)) 
  } else{
    return(list(plot_list = 
                  list(
                    plot_list_splat_rand, plot_list_splat_rand_150, plot_list_splat_rand_300,
                    p_splat_pca_umap_3, p_splat_pca_umap_25, p_splat_pca_umap_75,
                    p_splat_euclid, p_splat_pear, p_splat_spear),
                dist_list = 
                  list(umap_res_pc_75        = umap_res_pc_75,
                       umap_res_pc_25        = umap_res_pc_25,
                       umap_res_pc_3         = umap_res_pc_3,
                       euclid_dist           = euclid_dist, 
                       pear_dist             = pear_dist,
                       spearman_dist         = spearman_dist,
                       side_seq_ref_rand_50  = side_seq_ref_rand_50$dissim_final,
                       side_seq_ref_rand_150 = side_seq_ref_rand_150$dissim_final,
                       side_seq_ref_rand_300 = side_seq_ref_rand_300$dissim_final)))
  }
}



runDimReds <- function(sim_counts, group_labels, r = 50) {
  ### PCA => UMAP
  
  cat("PCA...")
  
  pca_res <-
    prcomp(t(as.matrix(sim_counts)))$x[,1:75]
  
  pc_75_dist <- dist(pca_res)
  pc_25_dist <- dist(pca_res[, 1:25])
  pc_3_dist  <- dist(pca_res[, 1:3])
  
  ### 3. Euclid UMAP ------------------------------------------------------------
  
  cat("euclidean and correlation...")
  
  euclid_dist = dist(t(as.matrix(sim_counts)))
  
  ### 4. Pearson UMAP -----------------------------------------------------------
  
  pear_dist <- 1 - abs(cor(sim_counts))
  
  ### 5. Spearman UMAP ----------------------------------------------------------
  
  spearman_dist <- 1 - abs(cor(sim_counts, method = "spearman"))
  
  
  ### 6. SIDERef Rand -----------------------------------------------------------
  
  cat("sideref rand...")
  
  side_seq_ref_rand_50 <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = 50,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)
  

  side_seq_ref_rand_150 <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = 150,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)
  
  side_seq_ref_rand_300 <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = 300,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)


    return(dist_list = 
             list(pc_75_dist            = pc_75_dist,
                  pc_25_dist            = pc_25_dist,
                  pc_3_dist             = pc_3_dist,
                  euclid_dist           = euclid_dist, 
                  pear_dist             = pear_dist,
                  spearman_dist         = spearman_dist,
                  side_seq_ref_rand_50  = side_seq_ref_rand_50$dissim_final,
                  side_seq_ref_rand_150 = side_seq_ref_rand_150$dissim_final,
                  side_seq_ref_rand_300 = side_seq_ref_rand_300$dissim_final))
}









