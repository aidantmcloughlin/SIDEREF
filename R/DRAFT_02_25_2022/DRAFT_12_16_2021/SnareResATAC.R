###############################################################################
### SNARE ATAC Results
###############################################################################


## set constants load data
set.seed(786231)
n_neighbors = 15
min_dist = 0.01

load(here("output/snare_sample_1500.RData"))
DefaultAssay(snare_sample) <- "ATAC"



### ---------------------------------------------------------------------------
### LSI Applied to Peak Counts Matrix -----------------------------------------
### ---------------------------------------------------------------------------


### TFIDF, SVD
snare_sample <- FindTopFeatures(snare_sample, min.cutoff = 10)
snare_sample <- RunTFIDF(snare_sample)
snare_sample <- RunSVD(snare_sample, reduction.name = "newlsi", 
                       reduction.key = "newlsi_")


reduct <- snare_sample@reductions$newlsi@cell.embeddings[,2:50]


sub_peaks_vec <-  c(150, 75, 25, 5)


## Tracking what I've saved so far:
reduct_sub_peaks_list <- list()
i <- 1
for(sub_peaks in sub_peaks_vec) {
  print(sub_peaks)
  
  snare_sample_sub_peaks <- 
    FindVariableFeatures(snare_sample, 
                         assay = "ATAC",
                         nfeatures = sub_peaks * 1e3)
  
  ## Renormalize.
  snare_sample_sub_peaks <- 
    RunTFIDF(snare_sample_sub_peaks[snare_sample_sub_peaks@assays$ATAC@var.features, ])
  
  snare_sample_sub_peaks <- RunSVD(snare_sample_sub_peaks, 
                                   reduction.name = "newlsi", 
                                   reduction.key = "newlsi_")
  
  
  reduct_sub_peaks_list[[i]] <- 
    snare_sample_sub_peaks@reductions$newlsi@cell.embeddings[,2:50]
  
  i <- i + 1
  
}


p_peak_counts_lsi_full <-
  snareUMAPPlot(reduct, 
                title = "TFIDF and SVD, All Peak Counts", 
                celltypes = snare_sample@meta.data$celltype)

p_sub_peaks_lsi_list <- list()
for(i in 1:4) {
  p_sub_peaks_lsi_list[[i]] <-
    snareUMAPPlot(reduct_sub_peaks_list[[i]][,2:30],
                  title = "TF-IDF and SVD, " %p% sub_peaks_vec[i] %p% "K Peak Counts",
                  celltypes = snare_sample@meta.data$celltype)
}


ggarrange(p_sub_peaks_lsi_list[[1]],
          p_sub_peaks_lsi_list[[2]],
          p_sub_peaks_lsi_list[[3]],
          p_sub_peaks_lsi_list[[4]], nrow = 4,
          common.legend = TRUE, legend = "right")

ggsave(here("output/figures/snare_atac_num_peaks.png"),
       width = 5, height = 14)





###### ------------------------------------------------------------------------
### Gene Activity Matrix-Based Embeddings -------------------------------------
###### ------------------------------------------------------------------------



###### ------------------------------------------------------------------------
### On subset of variable Genes -----------------------------------------------

set.seed(702918)

DefaultAssay(snare_sample) <- "ATAC_gene"

snare_sample <- NormalizeData(snare_sample)
snare_sample <- ScaleData(snare_sample,
                          do.center = FALSE,
                          do.scale  = TRUE)

gene_act_mat_var_genes <- 
  snare_sample@assays$ATAC_gene@scale.data


## Centered Data for PCA
snare_sample <- ScaleData(snare_sample,
                          do.center = TRUE,
                          do.scale  = TRUE)


gene_act_mat_var_genes_centered <- 
  snare_sample@assays$ATAC_gene@scale.data


sideref_top_genes_vec <- c(50, 150, 300)

cat("Variable Genes SIDEREF computations...")
for(g in sideref_top_genes_vec) {
  sideref_gene_act <- 
    SIDEseqRefSet(expr_matrix = gene_act_mat_var_genes, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = g,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = 150,
                  n_clust = 25,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = 23,
                  ## other
                  verbose = TRUE)
  
  
  save(sideref_gene_act, 
       file = here("output/snare_atac/" %p% 
                     "gene_act_variable_ref_rand_150_cells_" %p%
                     g %p% "_genes.RData"))  
}


### plot SIDEREF
load(here("output/snare_atac/gene_act_variable_ref_rand_150_cells_150_genes.RData"))

p_sideref_gene_act_subset <-
  snareUMAPPlot(as.dist(sideref_gene_act$dissim_final),
                title = "SIDEREF, 3K Variable ATAC Genes\n150 Random Cells, Top 150 Genes",
                celltypes = snare_sample@meta.data$celltype)



### non-SIDEREF methods -------------------------------------------------------

### PCA is on centered data:
pc_embed_gene_act <-
  prcomp(t(gene_act_mat_var_genes_centered))$x[,1:25] %>% data.frame()

p_pca_gene_act_subset <-
  snareUMAPPlot(pc_embed_gene_act[,2:5], 
                title = "5 PCs, 3K Variable ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)



### Pearson on full!
pear_dist <- 1 - abs(cor(gene_act_mat_var_genes))

p_snare_pear_subset <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson, 3K Variable ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)


### SVD, LSI
snare_var_gene_subset <- 
  subset(snare_sample, features = VariableFeatures(snare_sample))

snare_lsi <- RunTFIDF(snare_var_gene_subset)

snare_lsi <- RunSVD(snare_lsi,
                    reduction.name = "svdact", 
                    reduction.key = "svdact_")@reductions$svdact@cell.embeddings[,2:5]

p_snare_lsi_gene_subset <-
  snareUMAPPlot(snare_lsi, 
                title = "5 Col. LSI, 3K Variable ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)




### ---------------------------------------------------------------------------
### Computations run on all Genes ---------------------------------------------
### ---------------------------------------------------------------------------

snare_all_genes <- snare_sample
DefaultAssay(snare_all_genes) <- "ATAC_gene"
VariableFeatures(snare_all_genes) <- row.names(snare_sample@assays$ATAC_gene@counts)
snare_all_genes <- NormalizeData(snare_all_genes)
snare_all_genes <- ScaleData(snare_all_genes,
                             do.center = FALSE,
                             do.scale  = TRUE)


gene_act_mat_all_genes <- 
  snare_all_genes@assays$ATAC_gene@scale.data


## Centered Data for PCA
snare_all_genes <- ScaleData(snare_all_genes,
                             do.center = TRUE,
                             do.scale  = TRUE)


gene_act_mat_all_genes_centered <- 
  snare_all_genes@assays$ATAC_gene@scale.data



cat("All Genes SIDEREF computations...")
for(g in sideref_top_genes_vec) {
  sideref_gene_act <- 
    SIDEseqRefSet(expr_matrix = gene_act_mat_all_genes, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = g,
                  ## cell reference set selection parameters
                  selection_method = "random",
                  size_ref_set = 150,
                  n_clust = 25,
                  B = 0,
                  D = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = 23,
                  ## other
                  verbose = TRUE)
  
  
  save(sideref_gene_act, 
       file = here("output/snare_atac/" %p% 
                     "gene_act_allgenes_ref_rand_150_cells_" %p%
                     g %p% "_genes.RData"))  
}


### plot SIDEREF
load(here("output/snare_atac/gene_act_allgenes_ref_rand_150_cells_150_genes.RData"))

p_sideref_gene_act_full <-
  snareUMAPPlot(as.dist(sideref_gene_act$dissim_final),
                title = "SIDEREF, All ATAC Genes\n150 Random Cells, Top 150 Genes",
                celltypes = snare_sample@meta.data$celltype)



### non-SIDEREF methods -------------------------------------------------------


### PCA:
pc_embed_gene_act_all <-
  prcomp(t(gene_act_mat_all_genes_centered))$x[,1:25] %>% 
  data.frame()

p_pca_gene_act_full <-
  snareUMAPPlot(pc_embed_gene_act_all[,2:5], 
                title = "5 PCs, All ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)



### Pearson on full!
pear_dist <- 1 - abs(cor(gene_act_mat_all_genes))

p_snare_pear_full <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson, All ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)


snare_lsi <- RunTFIDF(snare_all_genes)

snare_lsi <- 
  RunSVD(snare_lsi,
         reduction.name = "svdact", 
         reduction.key = "svdact_")@reductions$svdact@cell.embeddings[,2:5]

p_snare_lsi_gene_all <-
  snareUMAPPlot(snare_lsi, 
                title = "5 Col LSI, All ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)



##### -------------------------------------------------------------------------
### Combine Individual Plots into Final Arrangements --------------------------
##### -------------------------------------------------------------------------


ggarrange(p_sideref_gene_act_subset, p_sideref_gene_act_full, 
          ncol = 2, common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/snare_sideref_gene_act_mat.png"),
       width = 9, height = 5)

ggarrange(p_pca_gene_act_subset, p_pca_gene_act_full,
          p_snare_lsi_gene_subset, p_snare_lsi_gene_all,
          p_snare_pear_subset, p_snare_pear_full,
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/gene_activity_mat_embeddings.png"),
       width = 7.5, height = 10)

