################################################################################
### Running Splatter simulations in functional manner
################################################################################

### Plot Function:
splatUMAPPlot <- function(data, title, celltypes,
                          pcs = NULL, plot_sil_score = TRUE) {
  
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
    ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
    labs(x="UMAP1", y="UMAP2", 
         col = ifelse(plot_sil_score, 
                      "Cell Type (Sil. Score)",
                      "Cell Type")) + 
    theme_bw() + 
    ggtitle(title)
  
  return(p)
}


runDimRedsAndPlot <- function(sim_counts, group_labels,
                              g = 50, r = 50) {
  ### PCA => UMAP
  
  cat("PCA plots...")
  
  p_splat_pca_umap_75 <-
    splatUMAPPlot(t(as.matrix(sim_counts)), title = "Euclidean Applied to 75 PCs", 
                  celltypes = group_labels,
                  pcs = 75,
                  plot_sil_score = FALSE)
  
  
  p_splat_pca_umap_25 <-
    splatUMAPPlot(t(as.matrix(sim_counts)), title = "Euclidean Applied to 25 PCs", 
                  celltypes = group_labels,
                  pcs = 25,
                  plot_sil_score = FALSE)
  
  
  p_splat_pca_umap_3 <-
    splatUMAPPlot(t(as.matrix(sim_counts)), title = "Euclidean Applied to 3 PCs", 
                  celltypes = group_labels,
                  pcs = 3,
                  plot_sil_score = FALSE)
  
  
  
  
  ### 3. Euclid UMAP ------------------------------------------------------------
  
  cat("euclidean and correlation plots...")
  
  p_splat_euclid <-
    splatUMAPPlot(t(as.matrix(sim_counts)), title = "Euclidean Distance", 
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  
  ### 4. Pearson UMAP -----------------------------------------------------------
  
  pear_dist <- 1 - abs(cor(sim_counts))
  
  p_splat_pear <-
    splatUMAPPlot(as.dist(pear_dist), title = "Pearson Distance", 
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  
  ### 5. Spearman UMAP ----------------------------------------------------------
  
  spearman_dist <- 1 - abs(cor(sim_counts, method = "spearman"))
  
  p_splat_spear <-
    splatUMAPPlot(as.dist(spearman_dist), title = "Spearman Distance", 
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  
  
  ### 6. SIDERef Rand -----------------------------------------------------------
  
  cat("sideref rand...")
  
  side_seq_ref_rand <- 
    SIDEseqRefSet(expr_matrix = sim_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = g,
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
  
  plot_list_splat_rand <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand$dissim_final),
                  title = "SIDEREF (random) 50 Cells 50 Top Genes",
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  
  
  side_seq_ref_rand <- 
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
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand$dissim_final),
                  title = "SIDEREF (random) 50 Cells 150 Top Genes",
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  side_seq_ref_rand <- 
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
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand$dissim_final),
                  title = "SIDEREF (random) 50 Cells 300 Top Genes",
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  return(plot_list = list(
    plot_list_splat_rand, plot_list_splat_rand_150, plot_list_splat_rand_300,
    p_splat_pca_umap_3, p_splat_pca_umap_25, p_splat_pca_umap_75,
    p_splat_euclid, p_splat_pear, p_splat_spear)
  )
}


###############################################################################
###############################################################################
### Load assorted simulations and give it a go.
###############################################################################
###############################################################################

## load simulation data
load(here("output/splatter_sideseq/splatter_sim1.RData"))
## Extract log counts
sim1 <- logNormCounts(sim1)
sim1_counts <- as.matrix(sim1@assays@data$logcounts)

splat_group_labels <- factor(sim1@colData@listData$Group)


levels(splat_group_labels) <- 
  c("D.E. Prob 2.0%, Low Var.", "D.E. Prob 2.0%, High Var.",
    "D.E. Prob 2.5%, Low Var.", "D.E. Prob 3.0%, Low Var.", 
    "D.E. Prob 4.0%, Low Var.")




plot_list_og <- 
  runDimRedsAndPlot(sim_counts = sim1_counts, 
                    group_labels = splat_group_labels)


ggarrange(plotlist = plot_list_og, ncol = 2, nrow = 4)


### Sim 2
## Extract log counts
sim2 <- logNormCounts(sim2)
sim2_counts <- as.matrix(sim2@assays@data$logcounts)

splat_group_labels <- factor(sim2@colData@listData$Group)


plot_list_shared_de_123 <-
  runDimRedsAndPlot(sim_counts = sim2_counts, 
                    group_labels = splat_group_labels)



### something to do: examine PC loadings of first 3 PCs

ggarrange(plotlist = plot_list_shared_de_123, ncol = 2, nrow = 4)




###############################################################################
###############################################################################
###############################################################################

sim3 <- logNormCounts(sim3)
sim3_counts <- as.matrix(sim3@assays@data$logcounts)

splat_group_labels <- factor(sim3@colData@listData$Group)


levels(splat_group_labels) <- 
  c("1: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 2, 3)", 
    "2: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 1, 3)",
    "3: D.E. Prob 3.0%\nShared D.E. Prob 2.0% (with Grp. 1, 2)",
    "4: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 5)",
    "5: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 4)",
    "6: D.E. Prob 4.0%\nNo Shared D.E.")


plot_list_shared_de_123_56 <-
  runDimRedsAndPlot(sim_counts = sim3_counts, 
                    group_labels = splat_group_labels)


ggarrange(plot_list_shared_de_123_56[[1]], plot_list_shared_de_123_56[[4]], 
          plot_list_shared_de_123_56[[2]], plot_list_shared_de_123_56[[5]], 
          plot_list_shared_de_123_56[[3]], plot_list_shared_de_123_56[[6]],
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/splatter_shared_de.png"),
       width = 8.5, height = 11)





###############################################################################
###############################################################################
###############################################################################


### Noisy Simulation



sim4 <- logNormCounts(sim4)
sim4_counts <- as.matrix(sim4@assays@data$logcounts)

splat_group_labels <- factor(sim4@colData@listData$Group)

#group.prob = c(0.1, 0.15, 0.3, 0.15, 0.1, 0.2)
#de.prob = c(0.01, 0.01, 0.03, 0.03, 0.03, 0.04)
#de.facLoc = de_fac_locs
#de.facScale = de_fac_scales
#de.downProb = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

#group_share_list <- list(c(1, 2, 3), c(4, 5))
#group_share_probs <- c(0.07, 0.04)
#group_share_fac_locs <- c(0.9, 0.9)
#group_share_fac_scales <- c(0.12, 0.12)
#group_share_downProbs <- c(0.5, 0.5)


levels(splat_group_labels) <- 
  c("1: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 2, 3)", 
    "2: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 1, 3)",
    "3: D.E. Prob 3.0%\nShared D.E. Prob 2.0% (with Grp. 1, 2)",
    "4: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 5)",
    "5: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 4)",
    "6: D.E. Prob 4.0%\nNo Shared D.E.")


plot_list_shared_de_123_56_noisy <-
  runDimRedsAndPlot(sim_counts = sim4_counts, 
                    group_labels = splat_group_labels)


ggarrange(plot_list_shared_de_123_56_noisy[[1]], plot_list_shared_de_123_56_noisy[[4]], 
          plot_list_shared_de_123_56_noisy[[2]], plot_list_shared_de_123_56_noisy[[5]], 
          plot_list_shared_de_123_56_noisy[[3]], plot_list_shared_de_123_56_noisy[[6]],
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/splatter_noisy.png"),
       width = 8.5, height = 11)



## TODO: get rid
(var(c(as.matrix(sim4@assays@data$TrueCounts))) - var(c(as.matrix(sim3@assays@data$TrueCounts))))  / 
  var(c(as.matrix(sim3@assays@data$TrueCounts)))


