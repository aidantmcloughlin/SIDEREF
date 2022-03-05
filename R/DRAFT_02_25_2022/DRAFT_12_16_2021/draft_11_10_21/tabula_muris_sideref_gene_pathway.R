#### SIDEREF Tabula Muris Data

## Load datasets 
liver_df <- 
  readRDS(here("data/tab_muris_bulk/Liver.rds"))

heart_aorta_df <- 
  readRDS(here("data/tab_muris_bulk/Heart_and_Aorta.rds"))

### clean data by size counts
liver_mat <- 
  as.matrix(t(t(liver_df$counts) / liver_df$size_factor))

heart_aorta_mat <- 
  as.matrix(t(t(heart_aorta_df$counts) / heart_aorta_df$size_factor))


## Load GO gene lists
## Liver genes
# GO_cholest_metab <- 
#   fread(here("data/tab_muris_bulk/" %p%
#                "GO_term_cholesterol_metabolism.csv"))
# 
# 
# ## heart, aorta genes
# 
# GO_card_musc_contract <- 
#   fread(here("data/tab_muris_bulk/" %p% 
#                "GO_term_cardiac_muscle_contraction.csv"))
# 
# GO_antigen_pres <- 
#   fread(here("data/tab_muris_bulk/" %p% 
#                "GO_term_antigen_presentation.csv"))

GO_fibro_prolif <- 
  fread(here("data/tab_muris_bulk/" %p% 
               "GO_term_fibroblast_proliferation.csv"))


bait_genes_fibro_prolif <-
  GO_fibro_prolif$Symbol[GO_fibro_prolif$Symbol %in% rownames(heart_aorta_mat)]



## To Seurat Object

scaleDataKeepGenes <- function(exp_mat, n_variable = 3000, genes_to_keep=NULL) {
  exp_mat_seurat = CreateSeuratObject(counts = exp_mat)
  ##heart
  exp_mat_seurat <- FindVariableFeatures(exp_mat_seurat, nfeatures = 3000)
  ## append Gene Pathway of Interest
  exp_mat_seurat@assays$RNA@var.features = c(exp_mat_seurat@assays$RNA@var.features,
                                             setdiff(genes_to_keep,
                                                     exp_mat_seurat@assays$RNA@var.features))
  exp_mat_seurat <- NormalizeData(exp_mat_seurat)
  exp_mat_seurat <- ScaleData(exp_mat_seurat)
  
  exp_mat <- exp_mat_seurat@assays$RNA@scale.data
  
  return(exp_mat)
  
}

heart_aorta_mat_clean <- scaleDataKeepGenes(heart_aorta_mat,
                                            n_variable = 500,
                                            genes_to_keep = bait_genes_fibro_prolif)



## sample genes to focus SIDEREF on the pathway.
n_pathway_genes = length(bait_genes_fibro_prolif)

genes_keep = c(bait_genes_fibro_prolif,
               base::sample(setdiff(rownames(heart_aorta_mat_clean), bait_genes_fibro_prolif),
                            n_pathway_genes))

top_genes = floor(length(genes_keep)/4)


## Run SIDEREF save results

print("fibro heart aorta SIDEREF computation")
fibro_heart_aorta_sideref_res <- 
  SIDEseqRefSet(expr_matrix = heart_aorta_mat_clean[genes_keep, ],
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = top_genes,
                selection_method = "cell_embed_sample",
                size_ref_set = 100,
                n_clust = 15, 
                B = 0, 
                D = 1,
                parallelize = TRUE,
                n_cores = 23,
                verbose = TRUE)

save(fibro_heart_aorta_sideref_res, file = "fibro_heart_aorta_sideref.RData")


