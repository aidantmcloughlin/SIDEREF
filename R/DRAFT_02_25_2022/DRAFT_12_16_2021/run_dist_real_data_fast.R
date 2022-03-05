N_CORES = 6
load(here("output/real_data/pbmc_seurat_variable_sample.RData"))
source(here("R/subgroup_composition_gene_contrib.R"))

## get cell types
load(here("output/real_data/pbmc_seurat.RData"))

group_labels = pbmc_seurat@meta.data$celltype

group_labels = group_labels[rownames(pbmc_seurat@meta.data) %in% colnames(pbmc_var_features)]


full_res <- 
  runDimReds(pbmc_var_features, group_labels = group_labels, 
             r = 100, g=c(1, 1, 1), 
             n_clust = 10, n_cores = N_CORES) 


## save
save(full_res, file = here("output/real_data/pbmc_dist_res_fast.RData"))