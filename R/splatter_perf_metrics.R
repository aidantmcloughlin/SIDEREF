library(here)
source(here("R/libraries.R"))
source(here('R/relative_group_dist_comps.R'))
source(here('R/perf_metrics_funcs.R'))

set.seed(1)


##### DATA LOAD, CLEAN, EVALUATE IN A LOOP:
splatter_metrics_res_list <- vector(mode="list", length = N_SPLATTER_SIMS_METRICS)
  
for(s in seq_len(N_SPLATTER_SIMS_METRICS)) {
  cat("Rep: " %p% s)
  
  ### Output file paths:
  
  splat_out_path <-
    here("output/splatter_sim/dist_res/rep" %p% s %p% "/")
  
  splat_out_files <- 
    list.files(splat_out_path)
  
  
  ## Load results
  for(f in splat_out_files) {
    load(splat_out_path %p% f,
         verbose = TRUE)  
  }
  
  
  ### Note: this is manually created dependent on the experiment format:
  splatter_global_group_map <-
    data.frame(Global_Group = c("Group1", "Group2",
                                "Global_Group" %p% 1:6)) %>% 
    mutate(cell_group_high = row_number()) %>% 
    mutate(global_group_label = c(
      "A [1] High Ind. D.E. Prob",
      "B [2] Low Ind. D.E. Prob",
      "C [3,4,5] High Ind. D.E., Low Shared D.E.",
      "D [6,7,8] High Ind. D.E., High Shared D.E.",
      "E [9,10,11] Low Ind D.E., High Shared D.E.",
      "F [12,13,14] Low Ind D.E., Low Shared D.E.",
      "G [15,16,17] Low Ind D.E., Low Shared D.E., High D.E. Variance",
      "H [18,19,20] High Ind D.E., Low Shared D.E., High D.E. Variance"
    )) %>% 
    
    mutate(hm_global_group_label = c(
      "A", "B", "C", "D", 
      "E", "F", "G", "H"
    ))
  
  ### Heatmaps ==================================================================
  meta_data <-
    meta_data %>% 
    mutate(cell_group_low = gsub("Group", "", 
                                 as.character(meta_data$Group)) %>% 
             as.numeric()) %>% 
    left_join(splatter_global_group_map)
  
  group_labels <- as.character(meta_data$cell_group_low)
  
  
  group_global_group_map <- 
    meta_data %>%
    dplyr::distinct(cell_group_low, cell_group_high) %>% 
    left_join(splatter_global_group_map) %>%
    dplyr::rename(cell_group_high1 = cell_group_high)
  
  grpwise_max_vio   <- sapply(
    append(main_dist_output$dist_list, othr_dist_output$dist_list),
    function(dist) {countGlobalGrpMaxViolationsMulti(
        dist, 
        group_labels = group_labels,
        group_global_group_map = group_global_group_map)
      })
  
  grpwise_min_vio   <- sapply(
    append(main_dist_output$dist_list, othr_dist_output$dist_list),
    function(dist) {countGlobalGrpMinViolationsMulti(
      dist, 
      group_labels = group_labels,
      group_global_group_map = group_global_group_map)
    })
  
  
  splatter_metrics_res_list[[s]] = 
    data.frame(method = names(grpwise_max_vio),
      rep = s,
      grpwise_max_vio = grpwise_max_vio,
      grpwise_min_vio = grpwise_min_vio
    )
  
}


### Reduce, summarize, table.
METHOD_LIST <- c("side_ref_g300_dist",
                 "pca_wtd25_dist",
                 "euclid_dist", "pear_dist", "spearman_dist",
                 "RAFSIL", "SIMLR_5_dims")

splatter_metrics_res_df <- 
  purrr::reduce(splatter_metrics_res_list, rbind) %>% 
  dplyr::filter(method %in% METHOD_LIST) %>% 
  group_by(method) %>% 
  summarize(grpwise_max_vio  = mean(grpwise_max_vio),
            grpwise_min_vio  = mean(grpwise_min_vio))

## reshuffle
idx=c()
for(m in METHOD_LIST) {
  idx <- c(idx, which(splatter_metrics_res_df$method==m))
}

splatter_metrics_res_df <- splatter_metrics_res_df[idx,]

splatter_metrics_res_df %>% 
  mutate(grpwise_max_vio = grpwise_max_vio*100,
         grpwise_min_vio = grpwise_min_vio*100)

