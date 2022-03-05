###############################################################################
### Low Level Cluster Analysis
###############################################################################
library(here)
source(here('R/libraries.R'))

## Constants
SPECTRAL_DIMS <- 15
file_name <- "tab_muris_full_0.05pct"

### load Tabula Muris data
load(here('output/tab_muris_sc/' %p% file_name %p% '_scale_var_genes_3000.RData'))
load(here('output/tab_muris_sc/' %p% file_name %p% '.RData'))

### load distance measures
load(here('output/tab_muris_sc/dist_res/' %p% 
            'full_tab_muris_dist_res_0.05pct_3000g.RData'))


for(dist_name in names(dist_res_tab_muris_0.05pct_3000g$dist_list[
  !grepl("pca", names(dist_res_tab_muris_0.05pct_3000g$dist_list))])) {
  print(dist_name)
  n <- ncol(dist_res_tab_muris_0.05pct_3000g$dist_list[[1]])
  
  l <- length(dist_res_tab_muris_0.05pct_3000g$dist_list)
  
  dissim <- dist_res_tab_muris_0.05pct_3000g$dist_list[[dist_name]]
  
  dissim <- dist(spectralEmbed(dissim, ndim = SPECTRAL_DIMS))
  
  mat <- matrix(0, nrow = n, ncol = n)
  
  mat[lower.tri(mat)]  <- c(dissim)
  
  mat[upper.tri(mat)] <-
    t(mat)[upper.tri(mat)]
  
  dist_res_tab_muris_0.05pct_3000g$dist_list[[l+1]] <- mat
  names(dist_res_tab_muris_0.05pct_3000g$dist_list)[l+1] <- 
    dist_name %p% "_spectral_d" %p% SPECTRAL_DIMS %p% "_dist"
}



### Low-level clustering pipeline:
data_seurat <-
  CreateSeuratObject(data_full, project = "SeuratProject", 
                     assay = "RNA")

data_seurat@assays$RNA@scale.data <- 
  data_full_variable

clusts <- list()

for(x in names(dist_res_tab_muris_0.05pct_3000g$dist_list)) {
  print(x)
  dist = dist_res_tab_muris_0.05pct_3000g$dist_list[[x]]
  
  rownames(dist) = colnames(data_seurat)
  colnames(dist) = colnames(data_seurat)
  
  
  nn_graphs <- FindNeighbors(dist, distance.matrix=TRUE)
  
  data_seurat@graphs$RNA_nn <- nn_graphs$nn
  data_seurat@graphs$RNA_snn <- nn_graphs$snn
  
  data_seurat <- FindClusters(data_seurat)
  
  clusts[[x]] <- data_seurat@meta.data$seurat_clusters
  
}

save(clusts, file = here("output/tab_muris_sc/low_level_clust_res.RData"))


