### snare rna plots


library(here)
source(here("R/side_seq/master.R"))

load(here("output/snareseq_final.RData"))
DefaultAssay(snare) <- "RNA"

## TODO: this belongs in a master script.
n_neighbors = 15
min_dist = 0.01
size_default = 0.9

## TODO: just save the sample somewhere
set.seed(209827)
snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))

save(snare_sample,
     file = here("output/snare_sideseq/snare_sample_1500.RData"))

snare_test <-
  snare[, setdiff(colnames(snare), colnames(snare_sample))] 

save(snare_test,
     file = here("output/snare_sideseq/snare_test.RData"))



snare_expr_mat <- snare_sample@assays$RNA@scale.data



snareUMAPPlot <- function(data, title, celltypes, size = 0.6,
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




plot_list_snare_cell_embed <- list()

i=1
for(r in c(25, 50, 75, 100, 150, 200, 300, 400)) {
  load(here("output/snare_sideseq/run1_" %p% r %p% "cell_embed_ref_set_snare_1500.RData"))
  plot_list_snare_cell_embed[[i]] <- 
    snareUMAPPlot(data  = as.dist(side_seq_ref_cell_embed$dissim_final), 
                  title = r %p% " Cells of 25 Subgroups Reference Set",
                  celltypes = snare_sample@meta.data$celltype,
                  plot_sil_score = FALSE,
                  size = 0.9) +
    theme(legend.position = "right")
  i=i+1
}


ggarrange(plotlist = plot_list_snare_cell_embed, nrow = 4, ncol = 2)

#ggsave("output/figures/snare_cell_embed_sidegref.png", height = 18, width = 14)





plot_list_snare_rand <- list()

i=1
for(r in c(25, 50, 75, 100, 150, 200, 300, 400)) {
  load(here("output/snare_sideseq/run1_" %p% r %p% "rand_ref_set_snare_1500.RData"))
  plot_list_snare_rand[[i]] <- 
    snareUMAPPlot(data  = as.dist(side_seq_ref_cell_embed$dissim_final), 
                  title = r %p% " Random Cells Reference Set",
                  celltypes = snare_sample@meta.data$celltype,
                  plot_sil_score = FALSE,
                  size = 0.9) + 
    theme(legend.position = "right")
  i=i+1
}

ggarrange(plotlist = plot_list_snare_rand, nrow = 4, ncol = 2)
#ggsave("output/figures/snare_rand_sidegref.png", height = 18, width = 14)



### Arrange with random and embed against each other
ggarrange(plot_list_snare_cell_embed[[1]], plot_list_snare_rand[[1]],
          plot_list_snare_cell_embed[[2]], plot_list_snare_rand[[2]],
          plot_list_snare_cell_embed[[3]], plot_list_snare_rand[[3]],
          plot_list_snare_cell_embed[[4]], plot_list_snare_rand[[4]],
          nrow = 4, ncol = 2,
          common.legend = TRUE, legend = "bottom")

ggsave("output/figures/snare_sideref_small_comps.png", height = 18, width = 11)


ggarrange(plot_list_snare_cell_embed[[5]], plot_list_snare_rand[[5]],
          plot_list_snare_cell_embed[[6]], plot_list_snare_rand[[6]],
          plot_list_snare_cell_embed[[7]], plot_list_snare_rand[[7]],
          plot_list_snare_cell_embed[[8]], plot_list_snare_rand[[8]],
          nrow = 4, ncol = 2,
          common.legend = TRUE, legend = "bottom")

ggsave("output/figures/snare_sideref_big_comps.png", height = 18, width = 11)



###############################################################################
### PCA => UMAP Comparisons
###############################################################################


p_snare_10_pca_umap <-
  snareUMAPPlot(t(as.matrix(snare_expr_mat)),
                title = "Euclidean Applied to 10 PCs",
                celltypes = snare_sample@meta.data$celltype,
                pcs = 10) +
  theme(legend.position = "right")

p_snare_25_pca_umap <-
  snareUMAPPlot(t(as.matrix(snare_expr_mat)),
                title = "Euclidean Applied to 25 PCs",
                celltypes = snare_sample@meta.data$celltype,
                pcs = 25) +
  theme(legend.position = "right")
p_snare_25_pca_umap

p_snare_50_pca_umap <-
  snareUMAPPlot(t(as.matrix(snare_expr_mat)), 
                title = "Euclidean Applied to 50 PCs", 
                celltypes = snare_sample@meta.data$celltype,
                size = test_set_point_size,
                pcs = 50) + 
  theme(legend.position = "right")

p_snare_100_pca_umap <-
  snareUMAPPlot(t(as.matrix(snare_expr_mat)),
                title = "Euclidean Applied to 100 PCs",
                celltypes = snare_sample@meta.data$celltype,
                pcs = 100) +
  theme(legend.position = "right")


ggarrange(p_snare_10_pca_umap, p_snare_25_pca_umap, p_snare_50_pca_umap, p_snare_100_pca_umap,
          ncol = 2, nrow = 2, 
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/snare_pca_compare.png"),
       width = 8, height = 8)

###############################################################################
###############################################################################
###
### TEST SET RESULTS (6582 SNAREseq cells)
###
###############################################################################
###############################################################################

load(here("output/snare_test.RData"))

snare_expr_mat <- snare_test@assays$RNA@scale.data

snare_test_cell_types <- snare_test@meta.data$celltype

test_set_point_size <- 0.7

### 1. FULL SIDESEQ -----------------------------------------------------------


## TODO: I think I need this later.
# splatter_sideseq_full_res <-
#   SIDEseqSimult(expr_matrix = sim1_counts, 
#                 diff_expr_method = "diff_expr_norm",
#                 similarity_method = "n_intersect",
#                 n_top_genes = 30,
#                 parallelize = TRUE)


load(here("output/snare_sideseq/side_seq_test_100_cells.RData"))
dissim_sideref <- stats::as.dist(side_seq_ref_rand$dissim_final)


p_snare_sideref <-
  snareUMAPPlot(dissim_sideref, title = "SIDEREF (random) 100 Cells, 300 Top Genes", 
                celltypes = snare_test_cell_types,
                size = 0.7)

# ggsave(here("output/figures/full_side_seq_snare_1500.png"),
#        width = 5, height = 3.5)



### 2. PCA => UMAP ------------------------------------------------------------

# 
# p_snare_10_pca_umap <-
#   snareUMAPPlot(t(as.matrix(snare_expr_mat)),
#                 title = "Euclidean Applied to 10 PCs",
#                 celltypes = snare_test_cell_types,
#                 pcs = 10) +
#   theme(legend.position = "right")
# 
# p_snare_25_pca_umap <-
#   snareUMAPPlot(t(as.matrix(snare_expr_mat)),
#                 title = "Euclidean Applied to 25 PCs",
#                 celltypes = snare_test_cell_types,
#                 pcs = 25) +
#   theme(legend.position = "right")
# p_snare_25_pca_umap
# 
p_snare_50_pca_umap <-
  snareUMAPPlot(t(as.matrix(snare_expr_mat)),
                title = "Euclidean Applied to 50 PCs",
                celltypes = snare_test_cell_types,
                size = test_set_point_size,
                pcs = 50) +
  theme(legend.position = "right")
# 
# p_snare_100_pca_umap <-
#   snareUMAPPlot(t(as.matrix(snare_expr_mat)),
#                 title = "Euclidean Applied to 100 PCs",
#                 celltypes = snare_test_cell_types,
#                 pcs = 100) +
#   theme(legend.position = "right")
# 
# 
# ggarrange(p_snare_10_pca_umap, p_snare_25_pca_umap, p_snare_50_pca_umap,
#           p_snare_100_pca_umap,
#           ncol = 2, nrow = 2, 
#           common.legend = TRUE, legend = "bottom")


### 3. Euclid UMAP ------------------------------------------------------------

p_snare_euclid <-
  snareUMAPPlot(t(as.matrix(snare_expr_mat)), 
                title = "Euclidean Distance", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)


### 4. Pearson UMAP -----------------------------------------------------------

pear_dist <- 1 - abs(cor(snare_expr_mat))

p_snare_pear <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson Distance", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)


### 5. Spearman UMAP ----------------------------------------------------------

spearman_dist <- 1 - abs(cor(snare_expr_mat, method = "spearman"))

p_snare_spear <-
  snareUMAPPlot(as.dist(spearman_dist), 
                title = "Spearman Distance", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)


### 6. PCA --------------------------------------------------------------------


p_snare_pca <-
  prcomp(t(snare_expr_mat))$x[,1:2] %>% data.frame()

## get sil
sil_scores <- cluster::silhouette(as.integer(factor(snare_test_cell_types)),
                                  dist = stats::dist(p_snare_pca, 
                                                     method = "euclidean"))[,3]
avg_sils <- sapply(1:length(unique(snare_test_cell_types$celltype)), 
                   function(i){mean(sil_scores[as.integer(factor(snare_test_cell_types)) == i])})

pca_types <- snare_test_cell_types

#levels(pca_types) <- paste0(levels(snare_test_cell_types), " (", round(avg_sils, 2), ")") 

p_snare_pca$pca_types <- pca_types

p_snare_pca <- p_snare_pca %>%  
  ggplot(aes(x = PC1, y = PC2, col = pca_types)) + 
  geom_point(size = test_set_point_size) + 
  labs(col = "Cell Type") + 
  theme_bw() + 
  ggtitle("Standard PCA (No UMAP)")



### Arrange together
ggarrange(p_snare_sideref,
          p_snare_euclid,
          p_snare_pca, 
          p_snare_pear,
          p_snare_50_pca_umap,
          p_snare_spear,
          nrow = 3, ncol = 2,
          common.legend = TRUE, 
          legend = "bottom")


ggsave(file = "output/figures/snare_comparison.png", 
       width = 12.5, height = 12)



###############################################################################
###############################################################################
### More or less genes? #######################################################
###############################################################################
###############################################################################

p_snare_side_gene_list <- list()
avg_sils <- c()


gene_list_sizes <- c(50, 150, 300, 450, 600, 900) 

p_snare_side_gene_list <- list()

i <- 1
for(g in gene_list_sizes) {
  if(g != 150) {
    load(here("output/snare_sideseq/side_seq_rand_" %p% g %p% "_genes.RData"))  
    
    dissim_sideseq <- stats::as.dist(side_seq_ref_rand$dissim_final)
    
    
    p_snare_side_gene_list[[i]] <-
      snareUMAPPlot(dissim_sideseq, title = "SIDEREF (random) 100 Cells " %p% 
                      g %p% " Top Genes", 
                    celltypes = snare_sample@meta.data$celltype,
                    size = size_default)
    
  } else{
    ## Original 
    load(here("output/snare_sideseq/run1_100rand_ref_set_snare_1500.RData"))
    
    dissim_sideseq <- stats::as.dist(side_seq_ref_rand$dissim_final)
    
    p_snare_side_gene_list[[i]] <-
      snareUMAPPlot(dissim_sideseq, title = "SIDEREF (random) 100 Cells 150 Top Genes", 
                    celltypes = snare_sample@meta.data$celltype,
                    size = size_default)
  }
  
  ### Average silhouette width
  umap_res <-
    uwot::umap(dissim_sideseq,
               n_neighbors = n_neighbors,
               min_dist = min_dist) 
  
  sil_scores <- cluster::silhouette(as.integer(factor(snare_sample@meta.data$celltype)),
                                    dist = stats::dist(umap_res, 
                                                       method = "euclidean"))[,3]
  
  avg_sils[i] <- mean(sil_scores)
  
  i <- i + 1
}

ggarrange(p_snare_side_gene_list[[1]], p_snare_side_gene_list[[4]],
          p_snare_side_gene_list[[2]], p_snare_side_gene_list[[5]],
          p_snare_side_gene_list[[3]], p_snare_side_gene_list[[6]],
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/snare_different_top_genes.png"),
       width = 9, height = 12)

avg_sils




### ---------------------------------------------------------------------------
### ---------------------------------------------------------------------------
### ALREADY COVERED -----------------------------------------------------------
### ---------------------------------------------------------------------------
### ---------------------------------------------------------------------------

## TODO, put in loop
plot_list_num_peaks <- list()

load(here("output/embeds/snare_full_pca_1500"))
plot_list_num_peaks[[1]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "All 244,544 Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_150Kpeaks_pca_1500"))
plot_list_num_peaks[[2]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "150,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_100Kpeaks_pca_1500"))
plot_list_num_peaks[[3]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "100,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_50Kpeaks_pca_1500"))
plot_list_num_peaks[[4]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "50,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_5Kpeaks_pca_1500"))
plot_list_num_peaks[[5]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "5,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

ggarrange(plotlist = plot_list_num_peaks, nrow = 5,
          common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/snare_atac_num_peaks.png"), 
       width = 5, height = 25)

save(plot_list_num_peaks, 
     file = here("output/figures/plot_objects/atac_num_peaks_pca.RData"))


### For RNA-seq ---------------------------------------------------------------
### TODO: clean pipeline

DefaultAssay(snare) <- "RNA"

gene_counts <- c(dim(snare@assays$RNA@counts)[1], 
                 25000, 10000, 2500, 500)


## TODO: test if we really need to resample each time.
for(g in gene_counts) {
  print(g)
  set.seed(209827)
  snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))
  DefaultAssay(snare) <- "RNA"
  ### Preprocess
  snare_sample <- FindVariableFeatures(snare_sample, nfeatures = g)
  snare_sample <- NormalizeData(snare_sample)
  snare_sample <- ScaleData(snare_sample)
  
  snare_sample <- RunPCA(snare_sample, npcs = 50, 
                         reduction.key = "pcs_", reduction.name = "pcs")
  
  
  snare_samp_pca <- snare_sample@reductions$pcs@cell.embeddings %>% as.data.frame()
  save(snare_samp_pca, file = here(paste0("output/embeds/snare_", g, "_genes_pca_1500.RData")))
}


plot_list_n_genes <- list()

## TODO loop
load(here("output/embeds/snare_33160_genes_pca_1500.RData"))
plot_list_n_genes[[1]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "All 33,160 Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_25000_genes_pca_1500.RData"))
plot_list_n_genes[[2]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "25,000 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_10000_genes_pca_1500.RData"))
plot_list_n_genes[[3]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "10,000 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_2500_genes_pca_1500.RData"))
plot_list_n_genes[[4]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "2,500 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_500_genes_pca_1500.RData"))
plot_list_n_genes[[5]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "500 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

ggarrange(plotlist = list(plot_list_n_genes[[1]], plot_list_num_peaks[[1]], 
                          plot_list_n_genes[[2]], plot_list_num_peaks[[2]], 
                          plot_list_n_genes[[3]], plot_list_num_peaks[[3]], 
                          plot_list_n_genes[[4]], plot_list_num_peaks[[4]], 
                          plot_list_n_genes[[5]], plot_list_num_peaks[[5]]),
          nrow = 5, ncol = 2,
          common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/snare_atac_gene_sub_vs_atac_sub.png"), 
       width = 12, height = 18)



ggarrange(plotlist = list( plot_list_num_genes[[1]], 
                            plot_list_num_genes[[2]], 
                            plot_list_num_genes[[3]], 
                           plot_list_num_genes[[4]], 
                            plot_list_num_genes[[5]]),
          nrow = 5, ncol = 1,
          common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/snare_atac_gene_sub.png"), 
       width = 5, height = 25)



###############################################################################
### Weighted Nearest Neighbors Plot
###############################################################################

## SIDEseq RNA data


## PCA => UMAP ATAC Peak Counts Matrix
load(here("output/embeds/snare_full_pca_1500"))


###############################################################################
### Testing hyperparameters of UMAP
###############################################################################

neighbors_vec <- c(5, 15, 25, 50)
min_dist_vec <- c(1e-3, 0.01, 0.1, 0.5)
plot_list <- list()

avg_sil_vec <- c()

i <- 1
for(n in neighbors_vec) {
  n_neighbors = n
  for(m in min_dist_vec) {
    print(i)
    min_dist = m
    umap_snare_samp_atac_pca <-
      uwot::umap(snare_samp_pca,
                 n_neighbors = n_neighbors,
                 min_dist = min_dist)
    
    
    ## compute Silhouette Scores
    sil_scores <- cluster::silhouette(as.integer(factor(snare_sample@meta.data$celltype)),
                                      dist = stats::dist(umap_snare_samp_atac_pca, 
                                                         method = "euclidean"))[,3]
    avg_sil_vec[i] <- mean(sil_scores)
    
    p <-
      data.frame(umap_snare_samp_atac_pca) %>% 
      mutate(celltype = snare_sample@meta.data$celltype) %>%
      ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
      labs(x="UMAP1", y="UMAP2", col = "Cell Type") + 
      theme_bw() + 
      ggtitle(paste(n, "Neighbors,", m, "Minimum Dist,", 
                    round(avg_sil_vec[i], 2), "Sil. Score"))
    
    
    
    plot_list[[i]] <- p
    
    i <- i+1
    
  }
}


ggarrange(plotlist = plot_list, ncol = 4, nrow = 4,
          common.legend = TRUE, legend = "right")


ggsave(here("output/figures/snare_atac_pca_hypers.png"), 
       width = 20, height = 14)


save(plot_list, 
     file = here("output/figures/plot_objects/atac_pca_hyper_grid.RData"))











