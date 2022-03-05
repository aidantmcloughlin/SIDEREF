### Running the distance computations for real data sets

## Computation over PBMC data sample.
N_CORES = 18
load(here("output/real_data/pbmc_seurat_variable_sample.RData"))
source(here("R/subgroup_composition_gene_contrib.R"))

## get cell types
load(here("output/real_data/pbmc_seurat.RData"))

group_labels = pbmc_seurat@meta.data$celltype



full_res <- 
  runDimReds(pbmc_var_features, group_labels = group_labels, 
             r = 100, g=c(150,300,450), 
             n_clust = 10, n_cores = N_CORES) 


## save
save(full_res, file = here("output/real_data/pbmc_dist_res.RData")) 



## Test for no genes 
N_CORES = 4

load(here("output/real_data/pbmc_seurat_variable_sample.RData"))
source(here("R/subgroup_composition_gene_contrib.R"))

## get cell types
load(here("output/real_data/pbmc_seurat.RData"))

group_labels = pbmc_seurat@meta.data$celltype
group_labels = group_labels[rownames(pbmc_seurat@meta.data) %in% colnames(pbmc_var_features)]



full_res <- 
  runDimReds(pbmc_var_features, group_labels = group_labels, 
             r = 100, g=c(2,2,2), 
             n_clust = 10, n_cores = N_CORES) 


## Computation over Pollen data sample.





