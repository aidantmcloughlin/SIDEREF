###############################################################################
###############################################################################
###
### Splatter Results
###
###
###############################################################################
###############################################################################
###############################################################################


library(here)
source(here("R/side_seq/master.R"))


n_neighbors <- 15
min_dist <- 0.01

def_n_top_genes <- 150



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



## load simulation data
load(here("output/splatter_sideseq/splatter_sim1.RData"))
## Extract log counts
sim1 <- logNormCounts(sim1)
sim1_counts <- as.matrix(sim1@assays@data$logcounts)

## Different Reference Set Sizes of Interest
ref_set_sizes <- c(5, 10, 25, 50, 75, 100, 150) 


## Loop of reference set simulations  
for(s in 1:5) {
  for(r in ref_set_sizes) {
    print(r)
    side_seq_ref_rand <- 
      SIDEseqRefSet(expr_matrix = sim1_counts, 
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = def_n_top_genes,
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
    
    
    save(side_seq_ref_rand, 
         file = here(paste0("output/splatter_sideseq/run", s, 
                            "_", r, "rand_ref_set_splatter.RData")))
    
    side_seq_ref_cell_embed <- 
      SIDEseqRefSet(expr_matrix = sim1_counts, 
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = def_n_top_genes,
                    ## cell reference set selection parameters
                    selection_method = "cell_embed_sample",
                    size_ref_set = r,
                    n_clust = 5,
                    B = 3,
                    D = 1,
                    ## parallelization parameters
                    parallelize = TRUE,
                    n_cores = detectCores()-2,
                    ## other
                    verbose = TRUE)
    
    
    save(side_seq_ref_cell_embed, 
         file = here(paste0("output/splatter_sideseq/run", s,
                            "_", r, "cell_embed_ref_set_splatter.RData")))
    
    if(r >= 10) {
      side_seq_ref_cell_embed <- 
        SIDEseqRefSet(expr_matrix = sim1_counts, 
                      diff_expr_method = "diff_expr_norm",
                      similarity_method = "n_intersect",
                      n_top_genes = def_n_top_genes,
                      ## cell reference set selection parameters
                      selection_method = "cell_embed_sample",
                      size_ref_set = r,
                      n_clust = 10,
                      B = 3,
                      D = 1,
                      ## parallelization parameters
                      parallelize = TRUE,
                      n_cores = detectCores()-2,
                      ## other
                      verbose = TRUE)
      
      
      save(side_seq_ref_cell_embed, 
           file = here(paste0("output/splatter_sideseq/run", s,
                              "_", r, "cell_embed_ref_set_splatter_10_clust.RData"))) 
    }
    
  }
  
  
}

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
                celltypes = splat_group_labels,
                plot_sil_score = FALSE)


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
                celltypes = splat_group_labels,
                plot_sil_score = FALSE)



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
                celltypes = splat_group_labels,
                plot_sil_score = FALSE)




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
                celltypes = splat_group_labels,
                plot_sil_score = FALSE)



### 2. PCA => UMAP ------------------------------------------------------------


p_splat_pca_umap_75 <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Applied to 75 PCs", 
                celltypes = splat_group_labels,
                pcs = 75,
                plot_sil_score = FALSE)


p_splat_pca_umap_25 <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Applied to 25 PCs", 
                celltypes = splat_group_labels,
                pcs = 25,
                plot_sil_score = FALSE)


p_splat_pca_umap_3 <-
  splatUMAPPlot(t(as.matrix(sim1_counts)), title = "Euclidean Applied to 3 PCs", 
                celltypes = splat_group_labels,
                pcs = 3,
                plot_sil_score = FALSE)

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


### 6. ZINBwave ---------------------------------------------------------------

## identify most variable genes for computational reasons

## Not included, performs poorly without apriori specifying groups

# gene_vars <- sim1_counts %>% rowVars
# 
# gene_vars <- order(gene_vars, decreasing = TRUE)
# 
# sim1_raw_counts_var_genes <- sim1@assays@data$counts[gene_vars[1:1e3], ]
# 
# as_experiment <- SummarizedExperiment(sim1_raw_counts_var_genes)
# 
# zinb_wave_embed <- zinbwave(as_experiment, K = 2, epsilon = 1e6)
# 
# W <- reducedDim(zinb_wave_embed)
# 
# p_splat_zinb <-
#   splatUMAPPlot(W, title = "ZINB-WaVE", celltypes = splat_group_labels)
# 
# ggplot(data.frame(W), aes(x=W1, y=W2)) + geom_point(aes(col = splat_group_labels))

### 7. PCA --------------------------------------------------------------------


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
          nrow = 3, ncol = 2)


ggsave(file = here("output/figures/splatter_comparison.png"), 
       width = 13, height = 10)



### Sensitivity of methods to main hyperparameter:
ggarrange(p_splat_side_full_50, p_splat_pca_umap_3,
          p_splat_side_full, p_splat_pca_umap_25,
          p_splat_side_full_500, p_splat_pca_umap_75, 
          
          nrow = 3, ncol = 2,
          common.legend = TRUE, legend = "right")

ggsave(file = "output/figures/splatter_side_pca_hypers.png", 
       width = 8, height = 8)

### 7. SIDEGREF Rand and Cell Embed  ----------------------------------------------------------------




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
                  celltypes = splat_group_labels,
                  plot_sil_score = FALSE)
  
  ## stratified sampling
  if(r >= 10) {
    
    load(here("output/splatter_sideseq/run1_" %p% r %p% 
                "cell_embed_ref_set_splatter_10_clust.RData")) 
    
    
    plot_list_splat_cell_embed_10_clust[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed$dissim_final),
                    title = r %p% " Cells of 10 Subgroups Reference Set", 
                    celltypes = splat_group_labels,
                    plot_sil_score = FALSE) 
  }
  
  ## random cell reference
  load(here("output/splatter_sideseq/run1_" %p% r %p% "rand_ref_set_splatter.RData"))
  plot_list_splat_rand[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand$dissim_final),
                  title = r %p% " Random Cells Reference Set", 
                  celltypes = splat_group_labels,
                  plot_sil_score = FALSE)
  
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

# ggarrange(plot_list_splat_cell_embed_10_clust[[2]],
#           plot_list_splat_cell_embed_10_clust[[3]],
#           plot_list_splat_cell_embed_10_clust[[4]],
#           plot_list_splat_cell_embed_10_clust[[5]],
#           nrow = 4, common.legend = TRUE, legend = "bottom")

ggsave("output/figures/splatter_sideref_ref_sizes.png", height = 18, width = 14)



ggarrange(plot_list_splat_cell_embed_5_clust[[2]], plot_list_splat_cell_embed_10_clust[[2]], plot_list_splat_rand[[2]], 
          plot_list_splat_cell_embed_5_clust[[3]], plot_list_splat_cell_embed_10_clust[[3]], plot_list_splat_rand[[3]], 
          plot_list_splat_cell_embed_5_clust[[4]], plot_list_splat_cell_embed_10_clust[[4]], plot_list_splat_rand[[4]], 
          plot_list_splat_cell_embed_5_clust[[5]], plot_list_splat_cell_embed_10_clust[[5]], plot_list_splat_rand[[5]], 
          plot_list_splat_cell_embed_5_clust[[6]], plot_list_splat_cell_embed_10_clust[[6]], plot_list_splat_rand[[6]], 
          nrow = 5, ncol = 3, common.legend = TRUE, legend = "bottom")

# ggarrange(plot_list_splat_cell_embed_10_clust[[2]],
#           plot_list_splat_cell_embed_10_clust[[3]],
#           plot_list_splat_cell_embed_10_clust[[4]],
#           plot_list_splat_cell_embed_10_clust[[5]],
#           nrow = 4, common.legend = TRUE, legend = "bottom")

ggsave("output/figures/splatter_sideref_ref_sizes_150_bot.png", height = 18, width = 14)







### 7. SIDEGREF Rand and Cell Embed  (Averaged) --------------


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
    
    print(mean(side_seq_ref_rand_avg))
    
  }
  
  
  ## stratified sampling 5 cluster
  plot_list_splat_cell_embed_5_clust_avg[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed_5_avg),
                  title = r %p% " Cells of 5 Subgroups Reference Set", 
                  celltypes = splat_group_labels,
                  plot_sil_score = FALSE)
  
  if(r >= 10) {
    plot_list_splat_cell_embed_10_clust_avg[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed_10_avg),
                    title = r %p% " Cells of 10 Subgroups Reference Set", 
                    celltypes = splat_group_labels,
                    plot_sil_score = FALSE)  
  }
  
  plot_list_splat_rand_avg[[i]] <- 
    splatUMAPPlot(data  = as.dist(side_seq_ref_rand_avg),
                  title = r %p% " Random Cells Reference Set", 
                  celltypes = splat_group_labels,
                  plot_sil_score = FALSE) 
  
  i <- i + 1
}

### And we can examine from 1 to 5 iterations of averaging (for 10 and 25)

plot_list_splat_rand_avg_iters <- list()
plot_list_splat_cell_embed_5_clust_avg_iters <- list()
plot_list_splat_cell_embed_10_clust_avg_iters <- list()

plot_avg_over_iters_r_list <- c(10, 25)
i <- 1
for(r in plot_avg_over_iters_r_list) {
  for(s in 1:5) {
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "cell_embed_ref_set_splatter.RData")) 
    
    side_seq_ref_cell_embed_5 <- side_seq_ref_cell_embed
    side_seq_ref_cell_embed_5_avg <-  
      1/(s+1) * (s * side_seq_ref_cell_embed_5_avg + 
                   side_seq_ref_cell_embed_5$dissim_final)
    
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "cell_embed_ref_set_splatter_10_clust.RData")) 
    
    side_seq_ref_cell_embed_10 <- side_seq_ref_cell_embed
    side_seq_ref_cell_embed_10_avg <-  
      1/(s+1) * (s * side_seq_ref_cell_embed_10_avg + 
                   side_seq_ref_cell_embed_10$dissim_final)
  
    ## random cell reference
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "rand_ref_set_splatter.RData"))
    
    side_seq_ref_rand_avg <-  
      1/(s+1) * (s * side_seq_ref_rand_avg + 
                   side_seq_ref_rand$dissim_final)
    
    
    
    ## stratified sampling 5 cluster
    plot_list_splat_cell_embed_5_clust_avg_iters[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed_5_avg),
                    title = r %p% " Cells of 5 Subgroups Reference Set\n" %p% 
                      s %p% "Iteration Average", 
                    celltypes = splat_group_labels,
                    plot_sil_score = FALSE)
    
    plot_list_splat_cell_embed_10_clust_avg_iters[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_cell_embed_10_avg),
                    title = r %p% " Cells of 10 Subgroups Reference Set\n" %p% 
                      s %p% "Iteration Average", 
                    celltypes = splat_group_labels,
                    plot_sil_score = FALSE)  
    
    
    plot_list_splat_rand_avg_iters[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_seq_ref_rand_avg),
                    title = r %p% " Random Cells Reference Set\n" %p% 
                      s %p% "Iteration Average", 
                    celltypes = splat_group_labels,
                    plot_sil_score = FALSE) 
    
    i <- i + 1
}}



## arrange and save.

### All together now
ggarrange(plot_list_splat_rand_avg_iters[[1]], plot_list_splat_rand_avg_iters[[2]],
          plot_list_splat_rand_avg_iters[[3]], plot_list_splat_rand_avg_iters[[4]],
          plot_list_splat_rand_avg_iters[[5]], plot_list_splat_rand_avg_iters[[6]], 
          plot_list_splat_rand_avg_iters[[7]], plot_list_splat_rand_avg_iters[[8]], 
          plot_list_splat_rand_avg_iters[[9]], plot_list_splat_rand_avg_iters[[10]],
          nrow = 2, ncol = 5, common.legend = TRUE, legend = "bottom")



ggarrange(plot_list_splat_cell_embed_5_clust_avg_iters[[1]], plot_list_splat_cell_embed_5_clust_avg_iters[[2]],
          plot_list_splat_cell_embed_5_clust_avg_iters[[3]], plot_list_splat_cell_embed_5_clust_avg_iters[[4]],
          plot_list_splat_cell_embed_5_clust_avg_iters[[5]], plot_list_splat_cell_embed_5_clust_avg_iters[[6]], 
          plot_list_splat_cell_embed_5_clust_avg_iters[[7]], plot_list_splat_cell_embed_5_clust_avg_iters[[8]], 
          plot_list_splat_cell_embed_5_clust_avg_iters[[9]], plot_list_splat_cell_embed_5_clust_avg_iters[[10]],
          nrow = 2, ncol = 5, common.legend = TRUE, legend = "bottom")

ggarrange(plot_list_splat_cell_embed_10_clust_avg_iters[[1]], plot_list_splat_cell_embed_10_clust_avg_iters[[2]],
          plot_list_splat_cell_embed_10_clust_avg_iters[[3]], plot_list_splat_cell_embed_10_clust_avg_iters[[4]],
          plot_list_splat_cell_embed_10_clust_avg_iters[[5]], plot_list_splat_cell_embed_10_clust_avg_iters[[6]], 
          plot_list_splat_cell_embed_10_clust_avg_iters[[7]], plot_list_splat_cell_embed_10_clust_avg_iters[[8]], 
          plot_list_splat_cell_embed_10_clust_avg_iters[[9]], plot_list_splat_cell_embed_10_clust_avg_iters[[10]],
          nrow = 2, ncol = 5, common.legend = TRUE, legend = "bottom")






ggsave(here("output/figures/splatter_sideref_avgs.png"), height = 18, width = 14)





################################################################################
### Brief Timing Study
################################################################################

ptm <- proc.time()
side_seq_ref_rand <- 
  SIDEseqRefSet(expr_matrix = sim1_counts, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = def_n_top_genes,
                ## cell reference set selection parameters
                selection_method = "random",
                size_ref_set = 50,
                n_clust = 5,
                B = 0,
                D = 1,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-2,
                ## other
                verbose = TRUE)

ct_elap1 <- (proc.time() - ptm)[3]


splat_time_df <- data.frame(method = "random", cells = 0, D = 0,
                            run = 0, time_elap = 0)
r= 10
method = "random"
for(d in 1:5) {
  print(d)
  for(s in 1:5) {
  ptm <- proc.time()
  side_seq_ref_rand <- 
    SIDEseqRefSet(expr_matrix = sim1_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = def_n_top_genes,
                  ## cell reference set selection parameters
                  selection_method = method,
                  size_ref_set = r,
                  n_clust = 5,
                  B = 0,
                  D = d,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = detectCores()-2,
                  ## other
                  verbose = TRUE)
  
  ct_elap <- (proc.time() - ptm)[3]
  splat_time_df <- splat_time_df %>% 
    bind_rows(data.frame(method = method, cells = r, D = d, run = s, time_elap = ct_elap))
  }
}


splat_time_df <- splat_time_df[-1, ]


splat_time_df %>% 
  group_by(method, cells, D) %>% 
  summarise(time_elap = mean(time_elap)) %>% 
  ggplot(aes(x = D, y = time_elap)) + geom_point()


## Versus all at once
splat_time_df_simult <- data.frame(method = "random", cells = 0,
                                   run = 0, time_elap = 0)

for(r in c(20, 30, 40, 50)) {
  print(r)
  for(s in 1:5) {
    ptm <- proc.time()
    side_seq_ref_rand <- 
      SIDEseqRefSet(expr_matrix = sim1_counts, 
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = def_n_top_genes,
                    ## cell reference set selection parameters
                    selection_method = method,
                    size_ref_set = r,
                    n_clust = 5,
                    B = 0,
                    D = 1,
                    ## parallelization parameters
                    parallelize = TRUE,
                    n_cores = detectCores()-2,
                    ## other
                    verbose = TRUE)
    
    ct_elap <- (proc.time() - ptm)[3]
    splat_time_df_simult <- splat_time_df_simult %>% 
      bind_rows(data.frame(method = method, cells = r, run = s, time_elap = ct_elap))
  }
}


splat_time_df_simult <- splat_time_df_simult[-1, ]


full_df <- 
  splat_time_df %>% 
  group_by(method, D) %>% 
  summarise(time_elap = median(time_elap)) %>% 
  mutate(cells = D * 10) %>% dplyr::select(-D) %>% 
  mutate(compute = "Average over samples of 10") %>%
  bind_rows(
    splat_time_df_simult %>% 
      group_by(method, cells) %>% 
      summarise(time_elap = median(time_elap)) %>% 
      mutate(compute = "All in one sample")
  )
  
full_df %>% 
  ggplot(aes(x = cells, y = time_elap, color = compute)) + 
  geom_point() + 
  theme_bw()


