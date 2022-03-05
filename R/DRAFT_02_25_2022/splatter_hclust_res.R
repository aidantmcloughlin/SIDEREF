### Splatter HClust computations and Results.
set.seed(1)
REPS <- 1
REP_SEEDS <- 50:(50+REPS-1)

SPECTRAL_DIMS <- c(5, 10, 25)

PCs <- c(25)

library(here)
source(here('R/libraries.R'))
source(here('R/hcclust.R'))
source(here('R/subgroup_composition_gene_contrib.R'))
source(here('R/dist_and_clust_pipelines.R'))

F_PATH <- "output/splatter_sim/"
SIM_NAME <- "splatter_sim_skin_eq_grps_" %p% 
  "lowvar"

for(f in seq_len(REPS)) {
  set.seed(REP_SEEDS[f])
  ## load splatter simulation
  load(here(F_PATH %p% SIM_NAME %p% "_rep" %p% f %p% ".RData"),
       verbose=TRUE)
  
  
  
  sim <- logNormCounts(splat_res$sim)
  sim_counts <- as.matrix(sim@assays@data$logcounts)
  ## remove genes with no expression
  sim_counts <- sim_counts[rowSums(sim_counts) > 0, ]
  
  splat_group_labels_low <- factor(sim@colData@listData$Group)
  
  splat_group_labels_high <- rep(0, length(splat_group_labels_low))
  
  subgroup_vars <- 
    names(rowData(sim))[grepl("shared_genes", names(rowData(sim)))]
  
  subgroup_vars <- 
    lapply(str_split(gsub("shared_genes", "", subgroup_vars), "_"),
           function(x) "Group" %p% x)
  
  for(i in seq_len(length(subgroup_vars))) {
    splat_group_labels_high[splat_group_labels_low %in% subgroup_vars[[i]]] <- 
      "Subgroup" %p% i
  }
  
  splat_group_labels_high[splat_group_labels_high == "0"] <- 
    as.character(splat_group_labels_low)[splat_group_labels_high == "0"]
  
  
  splat_group_labels_high <- 
    splat_group_labels_high %>% factor()
  
  ## distance computations
  splat_dist_res <- 
    splatDistResPipeline(sim_counts, splat_group_labels_low,
                         var_features=VAR_FEATURES,
                         G = G,
                         PCs = PCs,
                         spectral_dims = SPECTRAL_DIMS,
                         n_cores = N_CORES)
  
  
  ## hclust computations
  results_list <- list()
  
  ## From cell types
  splat_hclust_res_celltypes <- 
    lapply(names(splat_dist_res$dist_list),
           function(n) {
             return(
               RunHClustPipeline(
                 clusters = as.numeric(splat_group_labels_low), 
                 subg_inds = as.numeric(splat_group_labels_high),
                 hclust_dissim_mat = splat_dist_res$dist_list[[n]], 
                 p_dissim_mat      = splat_dist_res$dist_list$pca_25_dist,
                 cell_ids = colnames(sim_counts)))
           })
  
  names(splat_hclust_res_celltypes) <- 
    names(splat_dist_res$dist_list)
  
  results_list[["hclust_res_celltypes"]] <-
    splat_hclust_res_celltypes
  
  
  results_list[["hclust_ari_celltypes"]] <-
    lapply(splat_hclust_res_celltypes, 
           function(x) x$adj_rand_index_high)
  
  
  ### From cell level
  splat_hclust_res_celllevel <- 
    lapply(names(splat_dist_res$dist_list),
           function(n) {
             return(
               RunHClustPipeline(
                 clusters = seq_len(length(splat_group_labels_low)), 
                 subg_inds = as.numeric(splat_group_labels_high),
                 hclust_dissim_mat = splat_dist_res$dist_list[[n]], 
                 p_dissim_mat      = splat_dist_res$dist_list$pca_25_dist,
                 cell_ids = colnames(sim_counts)))
           })
  
  names(splat_hclust_res_celllevel) <- 
    names(splat_dist_res$dist_list)
  
  
  results_list[["hclust_res_celllevel"]] <-
    splat_hclust_res_celllevel
  
  results_list[["hclust_ari_celllevel"]] <-
    lapply(splat_hclust_res_celllevel, 
           function(x) x$adj_rand_index_high)
  
  results_list[["meta_data"]] <- 
    colData(splat_res$sim) %>% 
    data.frame() %>% 
    mutate(SubGroup = splat_group_labels_high)
  
  
  if(f == 1) {
    results_list[["dist_res"]] <- splat_dist_res}
  
  
  save(results_list,
       file = here(F_PATH %p% "hclust_res/" %p% SIM_NAME %p% "_rep" %p% 
                     f %p% "_hclustres.RData"))
  
}



