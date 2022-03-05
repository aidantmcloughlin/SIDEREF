###############################################################################
###
### Main Splatter Results
###
###############################################################################


set.seed(3984201)

n_neighbors <- 15
min_dist <- 0.01
def_n_top_genes <- 150


## load simulation data
load(here("output/splatter_sim1.RData"))
## Extract log counts
sim1 <- logNormCounts(sim1)
sim1_counts <- as.matrix(sim1@assays@data$logcounts)
## remove genes with no expression
sim1_counts <- sim1_counts[rowSums(sim1_counts) > 0, ]


################################################################################
### Formal Plots of the Simulation
################################################################################

splat_group_labels <- factor(sim1@colData@listData$Group)


levels(splat_group_labels) <- 
  c("D.E. Prob 2.0%, Low Var.", "D.E. Prob 2.0%, High Var.",
    "D.E. Prob 2.5%, Low Var.", "D.E. Prob 3.0%, Low Var.", 
    "D.E. Prob 4.0%, Low Var.")







### 1. FULL SIDESEQ -----------------------------------------------------------

splatter_sideseq_full_res <-
  SIDEseqSimult(expr_matrix = sim1_counts, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = def_n_top_genes,
                parallelize = TRUE,
                n_cores = 5)


save(splatter_sideseq_full_res, 
     file = here("output/splatter_sideseq/splatter_sideseq_full.RData"))

dissim_sideseq <- stats::as.dist(splatter_sideseq_full_res)


p_splat_side_full <- 
  splatUMAPPlot(dissim_sideseq, title = "SIDEseq (Top 150 Genes)", 
                celltypes = splat_group_labels)


splatter_sideseq_full_res_1000 <-
  SIDEseqSimult(expr_matrix = sim1_counts, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 1000,
                parallelize = TRUE,
                n_cores = 5)

dissim_sideseq <- stats::as.dist(splatter_sideseq_full_res_1000)


p_splat_side_full_1000 <- 
  splatUMAPPlot(dissim_sideseq, title = "SIDEseq (Top 1000 Genes)", 
                celltypes = splat_group_labels)



splatter_sideseq_full_res_500 <-
  SIDEseqSimult(expr_matrix = sim1_counts, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 500,
                parallelize = TRUE,
                n_cores = 5)

dissim_sideseq <- stats::as.dist(splatter_sideseq_full_res_500)


p_splat_side_full_500 <- 
  splatUMAPPlot(dissim_sideseq, title = "SIDEseq (Top 500 Genes)", 
                celltypes = splat_group_labels)




splatter_sideseq_full_res_50 <-
  SIDEseqSimult(expr_matrix = sim1_counts, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 50,
                parallelize = TRUE,
                n_cores = 5)

dissim_sideseq <- stats::as.dist(splatter_sideseq_full_res_50)


p_splat_side_full_50 <- 
  splatUMAPPlot(dissim_sideseq, title = "SIDEseq (Top 50 Genes)", 
                celltypes = splat_group_labels)



### 2. PCA => UMAP ------------------------------------------------------------


p_splat_pca_umap_75 <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Applied to 75 PCs", 
                celltypes = splat_group_labels,
                pcs = 75)


p_splat_pca_umap_25 <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Applied to 25 PCs", 
                celltypes = splat_group_labels,
                pcs = 25)


p_splat_pca_umap_3 <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Applied to 3 PCs", 
                celltypes = splat_group_labels,
                pcs = 3)

### 3. Euclid UMAP ------------------------------------------------------------

p_splat_euclid <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Distance", celltypes = splat_group_labels)


### 4. Pearson UMAP -----------------------------------------------------------

pear_dist <- 1 - abs(cor(sim1_counts))

p_splat_pear <-
  splatUMAPPlot(as.dist(pear_dist), title = "Pearson Distance", celltypes = splat_group_labels)


### 5. Spearman UMAP ----------------------------------------------------------

spearman_dist <- 1 - abs(cor(sim1_counts, method = "spearman"))

p_splat_spear <-
  splatUMAPPlot(as.dist(spearman_dist), title = "Spearman Distance", celltypes = splat_group_labels)



### 6. PCA --------------------------------------------------------------------


p_splat_pca <-
  prcomp(t(sim1_counts))$x[,1:2] %>% data.frame()

## get sil
sil_scores <- cluster::silhouette(as.integer(factor(splat_group_labels)),
                                  dist = stats::dist(p_splat_pca, 
                                                     method = "euclidean"))[,3]
avg_sils <- sapply(1:length(unique(splat_group_labels)), 
                   function(i){mean(sil_scores[as.integer(factor(splat_group_labels)) == i])})

pca_types <- splat_group_labels

levels(pca_types) <- paste0(levels(splat_group_labels), " (", round(avg_sils, 2), ")") 

p_splat_pca$pca_types <- pca_types

p_splat_pca <- p_splat_pca %>%  
  ggplot(aes(x = PC1, y = PC2, col = pca_types)) + geom_point() + 
  labs(col = "Cell Type") + 
  theme_bw() + 
  ggtitle("Standard PCA (No UMAP)")



### Arrange together: comparing our methods:
ggarrange(p_splat_side_full, p_splat_euclid,
          p_splat_pca, p_splat_pear, 
          p_splat_pca_umap_25, p_splat_spear,
          nrow = 3, ncol = 2,
          common.legend = TRUE,
          legend = "right")


ggsave(file = here("output/figures/splatter_comparison.png"), 
       width = 11, height = 10)



### Sensitivity of methods to main hyperparameter:
ggarrange(p_splat_side_full_50, p_splat_pca_umap_3,
          p_splat_side_full, p_splat_pca_umap_25,
          p_splat_side_full_500, p_splat_pca_umap_75, 
          
          nrow = 3, ncol = 2,
          common.legend = TRUE, legend = "right")

ggsave(file = "output/figures/splatter_side_pca_hypers.png", 
       width = 8, height = 8)





### SIDEREF Rand and Cell Embed Plots -----------------------------------------


plot_list_splat_cell_embed_5_clust  <- list()
plot_list_splat_cell_embed_10_clust <- list()
plot_list_splat_rand  <- list()


plot_r_list <- c(5, 10, 25, 50, 100, 150)

i <- 1
for(r in plot_r_list) {
  load(here("output/splatter_sideseq/run1_" %p% r %p% 
              "cell_embed_ref_set_splatter.RData")) 

  ## stratified sampling 5 cluster
  plot_list_splat_cell_embed_5_clust[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed$dissim_final),
                  title = r %p% " Cells of 5 Subgroups Reference Set", 
                  celltypes = splat_group_labels)
  
  ## stratified sampling
  if(r >= 10) {
    
    load(here("output/splatter_sideseq/run1_" %p% r %p% 
                "cell_embed_ref_set_splatter_10_clust.RData")) 
    
    
    plot_list_splat_cell_embed_10_clust[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed$dissim_final),
                    title = r %p% " Cells of 10 Subgroups Reference Set", 
                    celltypes = splat_group_labels) 
  }
  
  ## random cell reference
  load(here("output/splatter_sideseq/run1_" %p% r %p% "rand_ref_set_splatter.RData"))
  plot_list_splat_rand[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand$dissim_final),
                  title = r %p% " Random Cells Reference Set", 
                  celltypes = splat_group_labels)
  
  i <- i + 1
  
}

## arrange and save.

### All together now
ggarrange(plot_list_splat_cell_embed_5_clust[[1]], NULL, plot_list_splat_rand[[1]], 
          plot_list_splat_cell_embed_5_clust[[2]], plot_list_splat_cell_embed_10_clust[[2]], plot_list_splat_rand[[2]], 
          plot_list_splat_cell_embed_5_clust[[3]], plot_list_splat_cell_embed_10_clust[[3]], plot_list_splat_rand[[3]], 
          plot_list_splat_cell_embed_5_clust[[4]], plot_list_splat_cell_embed_10_clust[[4]], plot_list_splat_rand[[4]], 
          plot_list_splat_cell_embed_5_clust[[5]], plot_list_splat_cell_embed_10_clust[[5]], plot_list_splat_rand[[5]], 
          nrow = 5, ncol = 3, common.legend = TRUE, legend = "bottom")


ggsave("output/figures/splatter_sideref_ref_sizes.png", height = 18, width = 14)



# ggarrange(plot_list_splat_cell_embed_5_clust[[2]], plot_list_splat_cell_embed_10_clust[[2]], plot_list_splat_rand[[2]], 
#           plot_list_splat_cell_embed_5_clust[[3]], plot_list_splat_cell_embed_10_clust[[3]], plot_list_splat_rand[[3]], 
#           plot_list_splat_cell_embed_5_clust[[4]], plot_list_splat_cell_embed_10_clust[[4]], plot_list_splat_rand[[4]], 
#           plot_list_splat_cell_embed_5_clust[[5]], plot_list_splat_cell_embed_10_clust[[5]], plot_list_splat_rand[[5]], 
#           plot_list_splat_cell_embed_5_clust[[6]], plot_list_splat_cell_embed_10_clust[[6]], plot_list_splat_rand[[6]], 
#           nrow = 5, ncol = 3, common.legend = TRUE, legend = "bottom")
# 
# 
# ggsave("output/figures/splatter_sideref_ref_sizes_150_bot.png", height = 18, width = 14)







### 8. SIDEGREF Rand and Cell Embed  (Averaged) --------------


plot_list_splat_cell_embed_5_clust_avg  <- list()
plot_list_splat_cell_embed_10_clust_avg <- list()
plot_list_splat_rand_avg  <- list()


plot_r_list <- c(5, 10, 25, 50, 100, 150)

i <- 1
side_seq_ref_cell_embed_5_avg  <- 0
side_seq_ref_cell_embed_10_avg <- 0
side_seq_ref_rand_avg          <- 0

for(r in plot_r_list) {
  for(s in 1:5) {
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "cell_embed_ref_set_splatter.RData")) 
    
    side_seq_ref_cell_embed_5 <- side_seq_ref_cell_embed
    side_seq_ref_cell_embed_5_avg <-  
      1/(s+1) * (s * side_seq_ref_cell_embed_5_avg + 
                   side_seq_ref_cell_embed_5$dissim_final)
    
    print(mean(side_seq_ref_cell_embed_5_avg))
    
    ## stratified sampling
    if(r >= 10) {
      
      load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                  "cell_embed_ref_set_splatter_10_clust.RData")) 
      
      side_seq_ref_cell_embed_10 <- side_seq_ref_cell_embed
      side_seq_ref_cell_embed_10_avg <-  
        1/(s+1) * (s * side_seq_ref_cell_embed_10_avg + 
                     side_seq_ref_cell_embed_10$dissim_final)
      
      print(mean(side_seq_ref_cell_embed_10_avg))
    }
    
    ## random cell reference
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "rand_ref_set_splatter.RData"))
    
    side_seq_ref_rand_avg <-  
      1/(s+1) * (s * side_seq_ref_rand_avg + 
                   side_seq_ref_rand$dissim_final)
    
    
  }
  
  
  ## stratified sampling 5 cluster
  plot_list_splat_cell_embed_5_clust_avg[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed_5_avg),
                  title = r %p% " Cells of 5 Subgroups Reference Set", 
                  celltypes = splat_group_labels)
  
  if(r >= 10) {
    plot_list_splat_cell_embed_10_clust_avg[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed_10_avg),
                    title = r %p% " Cells of 10 Subgroups Reference Set", 
                    celltypes = splat_group_labels)  
  }
  
  plot_list_splat_rand_avg[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand_avg),
                  title = r %p% " Random Cells Reference Set", 
                  celltypes = splat_group_labels) 
  
  i <- i + 1
}



### All together now
ggarrange(plot_list_splat_cell_embed_5_clust_avg[[1]], NULL, plot_list_splat_rand_avg[[1]], 
          plot_list_splat_cell_embed_5_clust_avg[[2]], plot_list_splat_cell_embed_10_clust_avg[[2]], plot_list_splat_rand_avg[[2]], 
          plot_list_splat_cell_embed_5_clust_avg[[3]], plot_list_splat_cell_embed_10_clust_avg[[3]], plot_list_splat_rand_avg[[3]], 
          plot_list_splat_cell_embed_5_clust_avg[[4]], plot_list_splat_cell_embed_10_clust_avg[[4]], plot_list_splat_rand_avg[[4]], 
          plot_list_splat_cell_embed_5_clust_avg[[5]], plot_list_splat_cell_embed_10_clust_avg[[5]], plot_list_splat_rand_avg[[5]], 
          nrow = 5, ncol = 3, common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/splatter_sideref_avgs.png"), height = 18, width = 14)



side_seq_ref_cell_embed <- 
  SIDEseqRefSet(expr_matrix = sim1_counts, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = def_n_top_genes,
                ## cell reference set selection parameters
                selection_method = "cell_embed_sample",
                size_ref_set = 50,
                n_clust = 5,
                B = 0,
                D = 1,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-2,
                ## other
                verbose = TRUE)


plot_list_splat_pca_umap <- 
  splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed$dissim_final),
                title = "50 Cells PCA UMAP Single Reference Set", 
                celltypes = splat_group_labels) 
