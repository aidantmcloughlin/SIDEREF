library(here)
source(here('R/libraries.R'))
source(here("R/converge_analysis_helpers.R"))
source(here("R/distance_compute_helpers.R"))
source(here('R/relative_group_dist_comps.R'))

set.seed(1)

heatmapSimilarity <- function(dist1, dist2,
                              meta_data1, meta_data2,
                              cell_id_vec = NULL) {
  
  dist_df1 <-
    getSubDistDf(dist = dist1,
                 meta_data = meta_data1,
                 cell_id_vec = cell_id_vec)
  
  dist_df2 <-
    getSubDistDf(dist = dist2,
                 meta_data = meta_data2,
                 cell_id_vec = cell_id_vec)
  
  diffs <- dist_df1$dist - dist_df2$dist
  
  norm_diff <- norm(diffs , 
                    type = "2")
  
  
  quantile_diffs <- quantile(abs(diffs), c(.90, .95, .975, .99, 1))
  
  return(list(norm_diff = norm_diff, quantile_diffs = quantile_diffs))
  
}


SIDEREFHeatmapConvergenceDf <- 
  function(
    assay_data,
    meta_data,
    min_sample = 10,
    max_sample = 20,
    steps = 3,
    sample_sizes = NULL,
    reps = 1,
    seeds = NULL,
    g = 300,
    r = 100,
    n_pcs = NULL,
    n_clust = NULL,
    min_dist = 0.01,
    n_neighbors = 15,
    n_cores = 1,
    save_intermed_file = NULL) {
    
    if(is.null(seeds)) {
      seeds = seq_len(reps)
    }
    
    if(is.null(sample_sizes)) {
      sample_sizes <- 
        round(seq(min_sample, max_sample, length.out = steps))  
    } else{
      steps <- length(sample_sizes)
    }
    
    
    heatmap_convergence_res_df <- 
      data.frame(n = sample_sizes[2:steps],
                 n_prev = sample_sizes[1:(steps-1)],
                 group_dist_cells = length(unique(meta_data$group_label))**2,
                 norm_diff = rep(0, steps-1),
                 quantile_90 = rep(0, steps-1),
                 quantile_95 = rep(0, steps-1),
                 quantile_97.5 = rep(0, steps-1),
                 quantile_99 = rep(0, steps-1),
                 max_diff = rep(0, steps-1)
      )
    
    res_df_list <- lapply(seq_len(reps), function(x) heatmap_convergence_res_df)
    
    for(j in seq_len(reps)) {
      print("Rep " %p% j %p% " started.")
      set.seed(seeds[j])
      
      for(i in seq_len(steps)) {
        print(i)
        print("Sample Size: " %p% sample_sizes[i])
        ## Sample from Data
        samp_list <- sampleNPerClust(
          assay_data, 
          meta_data, 
          n=sample_sizes[i])
        
        print(dim(samp_list$assay_data_sample))
        ## Run SIDEREF:
        sideref_res <-
          SIDEREF(expr_matrix = samp_list$assay_data_sample,
                  n_top_genes = g,
                  R = r,
                  n_pcs = n_pcs,
                  n_clust = n_clust,
                  min_dist = min_dist,
                  n_neighbors = n_neighbors,
                  parallelize = TRUE,
                  n_cores = n_cores,
                  verbose = TRUE)
        
        ## Compute group dist df
        if(i == 1) {
          sideref_res_old <- sideref_res
          meta_data_old <- samp_list$meta_data_sample
        } else{
          heatmap_sim_res <-
            heatmapSimilarity(dist1 = sideref_res, 
                              dist2 = sideref_res_old,
                              meta_data1 = samp_list$meta_data_sample,
                              meta_data2 = meta_data_old)
          
          sideref_res_old <- sideref_res
          meta_data_old <- samp_list$meta_data_sample
          
          ## Compare and store to results df
          res_df_list[[j]]$norm_diff[i-1] <- heatmap_sim_res$norm_diff
          res_df_list[[j]][i-1, 5:ncol(res_df_list[[j]])] <- 
            heatmap_sim_res$quantile_diffs
          
          if(!is.null(save_intermed_file)) {
            save(res_df_list, file = save_intermed_file)
          }
        }
        ## end of i loop
      }
      ## end of j loop
    }
    
    return(res_df_list)
  }



###############################################################################
#### MAIN
###############################################################################

### Scripts.
source(here("R/tm_helper_functions.R"))

### Load data.

## Get full metadata
load(here("output/tab_muris_sc/tab_muris_full.RData"),
     verbose=TRUE)
rm(data_full)


## Get preprocessed assay data
load(here("output/tab_muris_sc/tab_muris_full_scale_var_genes_" %p% 
            VAR_FEATURES %p% ifelse(USE_DROPLET_ONLY, "_droplet_only", "") %p% 
            ".RData"),
     verbose=TRUE)

assay_data <- tm_preproc[[1]]
rownames(meta_data_full) <- meta_data_full$cell_id

meta_data <- meta_data_full[colnames(assay_data), ]
rm(meta_data_full, tm_preproc)

meta_data <- prepTMGroupLabels(meta_data, include_tissue = TRUE)

heatmap_convergence_res_df_list <-
  SIDEREFHeatmapConvergenceDf(
    assay_data,
    meta_data,
    min_sample = 2,
    max_sample = 50,
    steps = 25,
    sample_sizes = c(2, 5, 10, 20, 40, 60, 80, 100),
    reps = CONVERGE_RES_REPS,
    seeds = NULL,
    g = G,
    r = R,
    n_pcs = N_PCS,
    n_clust = N_CLUST,
    min_dist = MIN_DIST,
    n_neighbors = N_NEIGHBORS,
    n_cores = N_CORES,
    save_intermed_file = here("output/computations/sideref_converge_df_list.RData"))

