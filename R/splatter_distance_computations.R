### Splatter Dist computations and Results.
set.seed(1)

library(here)

source(here('R/libraries.R'))
source(here('R/relative_group_dist_comps.R'))
source(here("R/distance_compute_helpers.R"))


F_PATH <- "output/splatter_sim/"
DIST_REL_PATH <- "dist_res/"
SIM_NAME <- "splatter_sim"


runSplatterDistCompute <- 
  function(sim_name,
           f_path = F_PATH,
           d_rel_path = DIST_REL_PATH) {
  ## load splatter simulation
  load(here(f_path %p% sim_name %p% ".RData"),
       verbose=TRUE)
  
  sim <- logNormCounts(splat_res$sim)
  sim_counts <- as.matrix(sim@assays@data$logcounts)
  ## remove genes with no expression
  sim_counts <- sim_counts[rowSums(sim_counts) > 0, ]
  
  splat_group_labels_low <- factor(sim@colData@listData$Group)
  
  splat_group_labels_high <- rep(0, length(splat_group_labels_low))
  
  global_group_vars <- 
    names(rowData(sim))[grepl("shared_genes", names(rowData(sim)))]
  
  global_group_vars <- 
    lapply(str_split(gsub("shared_genes", "", global_group_vars), "_"),
           function(x) "Group" %p% x)
  
  for(i in seq_len(length(global_group_vars))) {
    splat_group_labels_high[splat_group_labels_low %in% global_group_vars[[i]]] <- 
      "Global_Group" %p% i
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
                         PCs = SPLAT_PCS,
                         spectral_dims = SPLAT_SPECTRAL_DIMS,
                         n_cores = N_CORES,
                         n_clust = N_CLUST,
                         n_pcs = N_PCS,
                         min_dist = MIN_DIST,
                         n_neighbors = N_NEIGHBORS,
                         save_loc = here(f_path %p% d_rel_path %p% sim_name))
  
  
  meta_data <- 
    colData(splat_res$sim) %>% 
    data.frame() %>% 
    mutate(Global_Group = splat_group_labels_high)
  
  
  save(meta_data,
       file = here(f_path %p% d_rel_path %p% sim_name %p% 
                     "_meta_data.RData"))
  
}

runSplatterDistCompute(sim_name = SIM_NAME,
                       d_rel_path = DIST_REL_PATH %p% "analysis/")


if(RUN_SPLATTER_PERFORMANCE_METRICS) {
  for(n in seq_len(N_SPLATTER_SIMS_METRICS)) {
    runSplatterDistCompute(sim_name = SIM_NAME %p% n,
                           d_rel_path = DIST_REL_PATH %p% "rep" %p% n %p% "/")
  }
}







