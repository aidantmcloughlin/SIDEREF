## snare ATAC plots

load(here("output/snareseq_final.RData"))
DefaultAssay(snare) <- "ATAC"

## TODO: this belongs in a master script.
n_neighbors = 15
min_dist = 0.01

set.seed(209827)
snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))


snareUMAPPlot <- function(data, title, celltypes,
                          pcs = NULL) {
  
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
  
  sil_scores <- cluster::silhouette(as.integer(factor(celltypes)),
                                    dist = stats::dist(umap_res, 
                                                       method = "euclidean"))[,3]
  avg_sils <- sapply(1:length(unique(celltypes)), function(i){mean(sil_scores[as.integer(factor(celltypes)) == i])})
  
  
  levels(celltypes) <- paste0(levels(celltypes), " (", round(avg_sils, 2), ")") 
  
  p <-
    data.frame(umap_res) %>% 
    mutate(celltype = celltypes) %>%
    ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
    labs(x="UMAP1", y="UMAP2", col = "Cell Type (Sil. Score)") + 
    theme_bw() + 
    ggtitle(title)
  return(p)
}



### ---------------------------------------------------------------------------
### Look at gene act counts compare to RNA ------------------------------------
### ---------------------------------------------------------------------------

gene_act_counts_raw <- snare_sample@assays$ATAC_gene@counts

summary(colMeans(gene_act_counts_raw))

summary(colMeans(snare_sample@assays$RNA@counts))



### ---------------------------------------------------------------------------
### ---------------------------------------------------------------------------
### Alternative scATAC Approaches ---------------------------------------------
### ---------------------------------------------------------------------------
### ---------------------------------------------------------------------------


DefaultAssay(snare_sample) <- "ATAC"
set.seed(209827)
snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))


### TFIDF, SVD
snare_sample <- FindTopFeatures(snare_sample, min.cutoff = 10)
snare_sample <- RunTFIDF(snare_sample)
snare_sample <- RunSVD(snare_sample, reduction.name = "newlsi", 
                       reduction.key = "newlsi_")


reduct <- snare_sample@reductions$newlsi@cell.embeddings[,2:50]

sum(gene_act_counts_raw)



#load(here("output/embeds/snare_full_pca_1500"))

## Run directly on TFIDF data to check
#atac_normalized <- snare_sample@assays$ATAC@data


## PCA on Full Data: 
#snare_samp_pca <- prcomp(t(atac_normalized))$x[,1:50]
#save(snare_samp_pca, file = here("output/embeds/snare_full_pca_1500_test.RData"))


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
  
  ## PCA
  #atac_normalized_sub_peaks <- snare_sample_sub_peaks@assays$ATAC@data
  
  #snare_samp_pca <- prcomp(t(atac_normalized_sub_peaks))$x[,1:100]
  
  #save(snare_samp_pca, file = here("output/embeds/snare_" %p% 
   #                                  sub_peaks %p% "Kpeaks_pca_1500.RData"))
  
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



load(here("output/embeds/snare_full_pca_1500_test.RData"))

p_peak_counts_pca <- 
  snareUMAPPlot(data  = snare_samp_pca[,2:50], 
                title = "Euclidean on 50 PCs, All Peak Counts",
                celltypes = snare_sample@meta.data$celltype)


pear_dist <- 1 - abs(cor(as.matrix(snare_sample@assays$ATAC@data)))

p_peak_counts_pear <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson Distance, All Peak Counts", 
                celltypes = snare_sample@meta.data$celltype)


ggarrange(p_peak_counts_lsi, p_peak_counts_pca, p_peak_counts_pear, nrow = 3)


ggsave(here("output/figures/umap_peak_counts_matrix.png"),
       width = 6, height = 11)


##############################################################################
### Weighted Nearest Neighbors
##############################################################################


### RNA Result
load(here("output/snare_sideseq/300cell_embed_ref_set_snare_1500.RData"))


umap_res_rna <-
  uwot::umap(as.dist(side_seq_ref_cell_embed$dissim_final),
             n_neighbors = n_neighbors,
             min_dist = min_dist) 

### ATAC Result
load(here("output/embeds/snare_full_pca_1500"))

umap_res_atac <-
  uwot::umap(snare_samp_pca,
             n_neighbors = n_neighbors,
             min_dist = min_dist) 


snare_sample@reductions[['umapside']] <- 
  CreateDimReducObject(embeddings = umap_res_rna, key = 'umapside_',
                       assay = "RNA")

snare_sample@reductions[['umappcaatac']] <- 
  CreateDimReducObject(embeddings = umap_res_atac, key = 'umappcaatac_',
                       assay = "ATAC")

snare_sample <- FindMultiModalNeighbors(
  snare_sample, reduction.list = list('umapside', 'umappcaatac'),
  dims.list = list(1:2, 1:2), modality.weight.name = "RNA.weight",
  k.nn = 10
)


### --------------------------------------------------------------------------
### --------------------------------------------------------------------------


p_peak_counts_pca <- 
  snareUMAPPlot(data  = snare_samp_pca, 
                title = "Euclidean on 50 PCs, All Peak Counts",
                celltypes = snare_sample@meta.data$celltype)




### Gene Activity Matrix

### On subset of variable Genes -----------------------------------------------

DefaultAssay(snare_sample) <- "ATAC_gene"

snare_sample <- NormalizeData(snare_sample)
snare_sample <- ScaleData(snare_sample,
                          do.center = FALSE,
                          do.scale  = FALSE)

gene_act_mat <- snare_sample@assays$ATAC_gene@scale.data


sideseq_gene_act <- 
  SIDEseqRefSet(expr_matrix = gene_act_mat, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = "random",
                size_ref_set = 35,#150,
                n_clust = 25,
                B = 0,
                D = 1,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = 5,
                ## other
                verbose = TRUE)


save(sideseq_gene_act, file = here("output/snare_atac/gene_act_variable_ref_rand_150.RData"))

sideseq_gene_act <-
  SIDEseqSimult(expr_matrix = gene_act_mat, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                parallelize = TRUE,
                n_cores = 6)

save(sideseq_gene_act, file = here("output/snare_atac/gene_act_variable_ref_full.RData"))


### PCA:
pc_embed_gene_act <-
  prcomp(t(gene_act_mat))$x[,1:25] %>% data.frame()

p_pca_gene_act_subset <-
  snareUMAPPlot(pc_embed_gene_act[,2:5], 
                title = "5 PCs, 3000 Variable ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)



### Pearson on full!
pear_dist <- 1 - abs(cor(gene_act_mat))

p_snare_pear_subset <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson, 3K Variable ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)


### SVD, LSI
test <- subset(snare_sample, features = VariableFeatures(snare_sample))
snare_lsi <- RunTFIDF(test)

snare_lsi <- RunSVD(snare_lsi,
                    reduction.name = "svdact", 
                    reduction.key = "svdact_")@reductions$svdact@cell.embeddings[,2:5]

p_snare_lsi_gene_subset <-
  snareUMAPPlot(snare_lsi, 
                title = "5 Col. LSI, 3K Variable ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)




### ---------------------------------------------------------------------------
### ---------------------------------------------------------------------------
### Computations run on all Genes ---------------------------------------------
### ---------------------------------------------------------------------------
### ---------------------------------------------------------------------------

snare_all_genes <- snare_sample
DefaultAssay(snare_all_genes) <- "ATAC_gene"
VariableFeatures(snare_all_genes) <- row.names(snare_sample@assays$ATAC_gene@counts)
snare_all_genes <- NormalizeData(snare_all_genes)
snare_all_genes <- ScaleData(snare_all_genes,
                             do.center = FALSE,
                             do.scale  = FALSE)

gene_act_mat_all <- snare_all_genes@assays$ATAC_gene@scale.data


### SIDESEQ
sideref_gene_act <- 
  SIDEseqRefSet(expr_matrix = gene_act_mat_all, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = "random",
                size_ref_set = 25, #150,
                n_clust = 25,
                B = 0,
                D = 1,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = 6,
                ## other
                verbose = TRUE)


save(sideseq_gene_act, file = here("output/snare_atac/gene_act_allgenes_ref_rand_150.RData"))

sideseq_gene_act <-
  SIDEseqSimult(expr_matrix = gene_act_mat, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                parallelize = TRUE,
                n_cores = 6)

save(sideseq_gene_act, file = here("output/snare_atac/gene_act_allgenes_variable_ref_full.RData"))


### Plotting SIDEseq ATAC Results:
load(here("output/snare_atac/gene_act_allgenes_variable_ref_full.RData"))

p_side_gene_act_full <-
  snareUMAPPlot(as.dist(sideseq_gene_act$dissim_final),
                title = "SIDEseq, All ATAC Genes",
                celltypes = snare_sample@meta.data$celltype)


load(here("output/snare_atac/gene_act_allgenes_ref_rand_150.RData"))

p_sideref_gene_act_full <-
  snareUMAPPlot(as.dist(sideseq_gene_act$dissim_final),
                title = "SIDEREF, All ATAC Genes\n150 Random Cells, Top 150 Genes",
                celltypes = snare_sample@meta.data$celltype)


load(here("output/snare_atac/gene_act_variable_ref_full.RData"))

p_side_gene_act_variable <-
  snareUMAPPlot(as.dist(sideseq_gene_act$dissim_final),
                title = "SIDEseq, 3,000 Variable ATAC Genes",
                celltypes = snare_sample@meta.data$celltype)

load(here("output/snare_atac/gene_act_variable_ref_rand_150.RData"))

p_sideref_gene_act_variable <-
  snareUMAPPlot(as.dist(sideseq_gene_act$dissim_final),
                title = "SIDEREF, 3,000 Variable ATAC Genes\n150 Random Cells, Top 150 Genes",
                celltypes = snare_sample@meta.data$celltype)


ggarrange(p_sideref_gene_act_full,
          p_sideref_gene_act_variable,
          ncol = 2, common.legend=TRUE, legend = "bottom")


ggsave(here("output/figures/snare_sideref_gene_act_mat.png"),
       width = 8, height = 4.2)


### PCA:
pc_embed_gene_act_all <-
  prcomp(t(gene_act_mat_all))$x[,1:25] %>% data.frame()

p_pca_gene_act_full <-
  snareUMAPPlot(pc_embed_gene_act_all[,2:5], 
                title = "5 PCs, All ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)



### Pearson on full!
pear_dist <- 1 - abs(cor(gene_act_mat_all))

p_snare_pear_full <-
  snareUMAPPlot(as.dist(pear_dist), 
                title = "Pearson, All ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)




snare_lsi <- RunTFIDF(snare_all_genes)

snare_lsi <- RunSVD(snare_lsi,
                    reduction.name = "svdact", 
                    reduction.key = "svdact_")@reductions$svdact@cell.embeddings[,2:5]

p_snare_lsi_gene_all <-
  snareUMAPPlot(snare_lsi, 
                title = "5 Col LSI, All ATAC Genes", 
                celltypes = snare_sample@meta.data$celltype)




ggarrange(p_pca_gene_act_subset, p_pca_gene_act_full,
          p_snare_lsi_gene_subset, p_snare_lsi_gene_all,
          p_snare_pear_subset, p_snare_pear_full,
          nrow = 3, ncol = 2, 
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/gene_activity_mat_embeddings.png"),
       width = 7.5, height = 10)

