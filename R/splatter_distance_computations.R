### Splatter Dist computations and Results.
set.seed(1)

SPECTRAL_DIMS <- c(5, 10, 25)

PCs <- c(25)

library(here)
source(here('R/libraries.R'))
source(here('R/relative_group_dist_comps.R'))
source(here("R/distance_compute_helpers.R"))


F_PATH <- "output/splatter_sim/"
SIM_NAME <- "splatter_sim"


## load splatter simulation
load(here(F_PATH %p% SIM_NAME %p% ".RData"),
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
                       n_cores = N_CORES,
                       n_clust = N_CLUST,
                       n_pcs = N_PCS,
                       save_loc = here(F_PATH %p% "dist_res/" %p% SIM_NAME))


meta_data <- 
  colData(splat_res$sim) %>% 
  data.frame() %>% 
  mutate(SubGroup = splat_group_labels_high)


save(meta_data,
     file = here(F_PATH %p% "dist_res/" %p% SIM_NAME %p% 
                   "_meta_data.RData"))





