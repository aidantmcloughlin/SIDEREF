#### SIDEREF Tabula Muris Data

## Load datasets 
print(here())

liver_df <- 
  readRDS(here("data/tab_muris_bulk/Liver.rds"))

heart_aorta_df <- 
  readRDS(here("data/tab_muris_bulk/Heart_and_Aorta.rds"))

### clean data by size counts
liver_mat <- 
  as.matrix(t(t(liver_df$counts) / liver_df$size_factor))

heart_aorta_mat <- 
  as.matrix(t(t(heart_aorta_df$counts) / heart_aorta_df$size_factor))


## 1924 cells:: Liver

## 619 cells: Heart-Aorta

## To Seurat Object
liver_seurat <- 
  CreateSeuratObject(counts = liver_mat)
heart_aorta_seurat <- 
  CreateSeuratObject(counts = heart_aorta_mat)

## Select 3000 Variable Genes
liver_seurat <- FindVariableFeatures(liver_seurat, nfeatures = 3000)
liver_seurat <- NormalizeData(liver_seurat)
liver_seurat <- ScaleData(liver_seurat)

liver_mat <- liver_seurat@assays$RNA@scale.data

##heart
heart_aorta_seurat <- FindVariableFeatures(heart_aorta_seurat, nfeatures = 3000)
heart_aorta_seurat <- NormalizeData(heart_aorta_seurat)
heart_aorta_seurat <- ScaleData(heart_aorta_seurat)

heart_aorta_mat <- heart_aorta_seurat@assays$RNA@scale.data

save(heart_aorta_mat, "test.RData")

## Run SIDEREF save results

print("heart aorta SIDEREF computation")
heart_aorta_sideref_res <- 
  SIDEseqRefSet(expr_matrix = heart_aorta_mat,
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 250,
                selection_method = "cell_embed_sample",
                size_ref_set = 150,
                n_clust = 150, 
                B = 0, 
                D = 1,
                parallelize = TRUE,
                n_cores = 23,
                verbose = TRUE)

save(heart_aorta_sideref_res, "heart_aorta_sideref.RData")

print("liver SIDEREF computation")
liver_sideref_res <- 
  SIDEseqRefSet(expr_matrix = liver_mat,
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 250,
                selection_method = "cell_embed_sample",
                size_ref_set = 150,
                n_clust = 150, 
                B = 0, 
                D = 1,
                parallelize = TRUE,
                n_cores = 23,
                verbose = TRUE)


save(liver_sideref_res, "liver_sideref.RData")



