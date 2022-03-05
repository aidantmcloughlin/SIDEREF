### Hclust Res Figures


library(here)
source(here("R/libraries.R"))



## load results data sets
tm_droplet_hclust_res_list <- 
  loadRData(here("output/tab_muris_sc/hclust_results_list_0.05pct_droplet_and_facs.RData"))

tm_droplet_hclust_res_list$hclust_ari_celltypes
tm_droplet_hclust_res_list$hclust_ari_celllevel
tm_droplet_hclust_res_list$hclust_ari_pca10
tm_droplet_hclust_res_list$hclust_ari_baseline

