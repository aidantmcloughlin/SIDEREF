###############################################################################
### SNARE RNA Results
###############################################################################


## constants
set.seed(727452)
n_neighbors = 15
min_dist = 0.01
size_default = 0.9
umap_reruns <- 5

## Load Training and Test Splits Data

load(here("output/snare_sample_1500.RData"))
load(here("output/snare_test.RData"))


DefaultAssay(snare_sample) <- "RNA"
DefaultAssay(snare_test)   <- "RNA"

snare_sample <- NormalizeData(snare_sample)
## This effectively just extracts normalized counts for variable genes
snare_sample <- ScaleData(snare_sample,
                          do.center = FALSE,
                          do.scale  = TRUE)

snare_test <- NormalizeData(snare_test)
## This effectively just extracts normalized counts for variable genes
snare_test <- ScaleData(snare_test,
                        do.center = FALSE,
                        do.scale  = TRUE)

snare_expr_mat <- snare_sample@assays$RNA@scale.data

snare_test_expr_mat <- snare_test@assays$RNA@scale.data


plot_list_snare_cell_embed <- list()


i <- 1

for(r in c(25, 50, 75, 100, 150, 200, 300, 400)) {
  load(here("output/snare_sideseq/run" %p% sample(1:5, 1) %p%
              "_" %p% r %p% "cell_embed_ref_set_snare_1500.RData"))
  
  dissim_sideref <- stats::as.dist(side_seq_ref_cell_embed$dissim_final)
  
  plot_list_snare_cell_embed[[i]] <- 
    snareUMAPPlot(data  = as.dist(dissim_sideref), 
                  title = r %p% " Cells of 25 Subgroups Reference Set",
                  celltypes = snare_sample@meta.data$celltype,
                  plot_sil_score = FALSE) +
    theme(legend.position = "right")

  
  i <- i + 1
}





plot_list_snare_rand <- list()


i <- 1

for(r in c(25, 50, 75, 100, 150, 200, 300, 400)) {
  load(here("output/snare_sideseq/run" %p% sample(1:5, 1) %p%
              "_" %p% r %p% "rand_ref_set_snare_1500.RData"))
  
  dissim_sideref <- stats::as.dist(side_seq_ref_rand$dissim_final)
  
  plot_list_snare_rand[[i]] <- 
    snareUMAPPlot(data  = as.dist(dissim_sideref), 
                  title = r %p% " Random Cells Reference Set",
                  celltypes = snare_sample@meta.data$celltype,
                  plot_sil_score = FALSE) + 
    theme(legend.position = "right")
  
  
  i <- i+1
}


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


### ---------------------------------------------------------------------------
### Full SIDEseq results
### ---------------------------------------------------------------------------

load(here("output/snare_sideseq/side_seq_snare_1500.RData"))

snareUMAPPlot(as.dist(full_side_seq),
              title = "SIDEseq",
              celltypes = snare_sample@meta.data$celltype,
              size = size_default,
              plot_sil_score = FALSE) +
  theme(legend.position = "right")


ggsave(here("output/figures/full_side_seq_snare_1500.png"),
       width = 5, height = 3.5)


###############################################################################
### More or less genes? #######################################################
###############################################################################


gene_list_sizes <- c(50, 150, 300, 450, 600, 900) 


p_snare_side_gene_list <- list()

i <- 1
for(g in gene_list_sizes) {
  print(g)
  
  load(here("output/snare_sideseq/side_seq_rand_" %p% g %p% "_genes.RData"))  
  
  dissim_sideseq <- stats::as.dist(side_seq_ref_rand$dissim_final)
  
  
  p_snare_side_gene_list[[i]] <-
    snareUMAPPlot(dissim_sideseq, title = "SIDEREF (random) 150 Cells " %p% 
                    g %p% " Top Genes", 
                  celltypes = snare_sample@meta.data$celltype,
                  size = size_default)
  
  
  i <- i + 1
}

ggarrange(p_snare_side_gene_list[[1]], p_snare_side_gene_list[[4]],
          p_snare_side_gene_list[[2]], p_snare_side_gene_list[[5]],
          p_snare_side_gene_list[[3]], p_snare_side_gene_list[[6]],
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/snare_different_top_genes.png"),
       width = 9, height = 12)





###############################################################################
###
### TEST SET RESULTS
###
###############################################################################


snare_test_cell_types <- snare_test@meta.data$celltype

test_set_point_size <- 0.7

### 1. SIDEREF ----------------------------------------------------------------



load(here("output/snare_sideseq/side_ref_test_150_cells_150_genes.RData"))
dissim_sideref <- stats::as.dist(side_seq_ref_rand$dissim_final)


p_snare_sideref <-
  snareUMAPPlot(dissim_sideref, title = "SIDEREF (random) 150 Cells, 150 Top Genes", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)




### 2. PCA => UMAP ------------------------------------------------------------

p_snare_25_pca_umap <-
  snareUMAPPlot(t(as.matrix(snare_test_expr_mat)),
                title = "Euclidean Applied to 25 PCs",
                celltypes = snare_test_cell_types,
                size = test_set_point_size,
                pcs = 25) +
  theme(legend.position = "right")


### 3. Euclid UMAP ------------------------------------------------------------

p_snare_euclid <-
  snareUMAPPlot(t(as.matrix(snare_test_expr_mat)), 
                title = "Euclidean Distance", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)


### 4. Pearson UMAP -----------------------------------------------------------

pear_dist <- 1 - abs(cor(snare_test_expr_mat))

p_snare_pear <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson Distance", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)


### 5. Spearman UMAP ----------------------------------------------------------

spearman_dist <- 1 - abs(cor(snare_test_expr_mat, method = "spearman"))

p_snare_spear <-
  snareUMAPPlot(as.dist(spearman_dist), 
                title = "Spearman Distance", 
                celltypes = snare_test_cell_types,
                size = test_set_point_size)


### 6. PCA --------------------------------------------------------------------


p_snare_pca <-
  prcomp(t(snare_test_expr_mat))$x[,1:2] %>% data.frame()


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
          p_snare_25_pca_umap,
          p_snare_spear,
          nrow = 3, ncol = 2,
          common.legend = TRUE, 
          legend = "bottom")


ggsave(file = "output/figures/snare_comparison.png", 
       width = 12.5, height = 12)





