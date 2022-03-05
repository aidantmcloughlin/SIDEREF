library(here)
source(here('R/libraries.R'))
source(here("R/distance_compute_helpers.R"))
source(here('R/relative_group_dist_comps.R'))

set.seed(100)

getSubDistDf <- 
  function(dist,
           meta_data,
           cell_id_vec = NULL) {
    
    dist <-
      selectFullDistIdx(
        cell_id_vec = cell_id_vec,
        dist = dist,
        meta_data)
    group_labels <- getGroupLabels(meta_data = meta_data, 
                                   cell_id_vec = cell_id_vec)
    
    dist_df <- getDistDf(group_labels, 
                         dist, 
                         gather = TRUE)
    
    return(dist_df)
  }




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
  
  
  pct_diffs <- quantile(
    abs(diffs), c(.90, .95, .975, .99, 1))
  
  return(list(norm_diff = norm_diff, 
              pct_diffs = pct_diffs))
  
}


sampleNPerClust <- 
  function(
    assay_data,
    meta_data, 
    n
  ) {
    group_based_idx <- 
      lapply(unique(meta_data$group_label),
             function(x) which(meta_data$group_label == x))
    
    
    group_idx_shuf <-
      lapply(group_based_idx, function(x){
        sample(x, size = min(n, length(x)), replace=FALSE)})
    
    cell_idx <- c()
    
    for(g in seq_len(length(group_idx_shuf))) {
      cell_idx <- c(cell_idx, 
                    group_idx_shuf[[g]])
    }
    
    assay_data_sample <- assay_data[, cell_idx]
    meta_data_sample <- meta_data[cell_idx, ]
    
    return(list(assay_data_sample = assay_data_sample,
                meta_data_sample  = meta_data_sample))
    
  }


SIDEREFHeatmapConvergenceDf <- 
  function(
    assay_data,
    meta_data,
    min_sample = 10,
    max_sample = 20,
    steps = 3,
    reps = 1,
    g = 300,
    r = 100,
    n_clust = 20,
    n_cores = 1,
    save_intermed_file = NULL) {
    
    sample_sizes <- 
      round(seq(min_sample, max_sample, length.out = steps))
    
    heatmap_convergence_res_df <- 
      data.frame(n = sample_sizes[1:steps],
                 group_dist_cells = length(unique(meta_data$group_label))**2,
                 norm_diff = rep(0, steps),
                 pctile_90 = rep(0, steps),
                 pctile_95 = rep(0, steps),
                 pctile_97.5 = rep(0, steps),
                 pctile_99 = rep(0, steps),
                 max_diff = rep(0, steps)
      )
    
    res_df_list <- lapply(seq_len(reps), function(x) heatmap_convergence_res_df)
    
    for(j in seq_len(reps)) {
      print("Rep " %p% j %p% " started.")
      for(i in seq_len(steps)) {
        print(i)
        print("Sample Size: " %p% sample_sizes[i])
        ## Sample from Data
        samp_list1 <- sampleNPerClust(
          assay_data, 
          meta_data, 
          n=sample_sizes[i])
        
        
        print(dim(samp_list1$assay_data_sample))
        ## Run SIDEREF:
        sideref_res1 <-
          SIDEREF(expr_matrix = samp_list1$assay_data_sample,
                  n_top_genes = g,
                  R = r,
                  n_clust = n_clust,
                  parallelize = TRUE,
                  n_cores = n_cores,
                  verbose = TRUE)
        
        sideref_res2 <-
          SIDEREF(expr_matrix = samp_list1$assay_data_sample,
                  n_top_genes = g,
                  R = r,
                  n_clust = n_clust,
                  parallelize = TRUE,
                  n_cores = n_cores,
                  verbose = TRUE)
        
        ## Compute group dist df
        heatmap_sim_res <-
          heatmapSimilarity(dist1 = sideref_res1, 
                            dist2 = sideref_res2,
                            meta_data1 = samp_list1$meta_data_sample,
                            meta_data2 = samp_list1$meta_data_sample,)
        
        ## Compare and store to results df
        res_df_list[[j]]$norm_diff[i] <- heatmap_sim_res$norm_diff
        res_df_list[[j]][i, 4:ncol(res_df_list[[j]])] <- 
          heatmap_sim_res$pct_diffs
        
        if(!is.null(save_intermed_file)) {
          save(res_df_list, file = save_intermed_file)
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
load(here("output/tab_muris_sc/preproc_data/tab_muris_full.RData"),
     verbose=TRUE)
rm(data_full)


## Get preprocessed assay data
load(here("output/tab_muris_sc/preproc_data/tab_muris_full_scale_var_genes_" %p% 
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
    steps = 20,
    reps = 1, 
    g = G,
    r = R,
    n_clust = N_CLUST,
    n_cores = N_CORES,
    save_intermed_file = here("output/sideref_converge_df_list.RData"))


# 
# ###############################################################################
# ### Plot Results
# 
# load(here("output/sideref_converge_df_list.RData"))
# 
# ### average over runs
# AVG_RUNS <- FALSE
# 
# if(AVG_RUNS) {
#   nonzero_reps <- sum(sapply(res_df_list, function(x) x$norm_diff[nrow(x)] > 0))
#   
#   avg_res_df <- Reduce("+", res_df_list) / nonzero_reps  
# } else{
#   avg_res_df <- res_df_list[[1]]
# }
# 
# last_entry <- max(which(avg_res_df$norm_diff > 0))
# 
# ### Plot Res
# avg_res_df[1:last_entry, ] %>%
#   mutate(mse = norm_diff**2 / group_dist_cells) %>% 
#   gather(cat, val, pctile_95:mse) %>%
#   ggplot(aes(x=n, y=val, color=cat)) +
#   geom_line() + 
#   geom_point(size=0.5) +
#   theme_bw() + 
#   ggtitle("Groupwise Heatmap Differences in Tabula Muris")
# 


## TODO: remove bootstrap (and keep random averaging option to sideref?)

## TODO: name changes to SIDEREF functions. 

## TODO: SIDEREF progress bar?

## TODO: SIDEREF purity scores irrelevant now? [run through code and remove list calls of sideref.]

## TODO: cleanup TM preproc and rename scripts to more intuitive.

## TODO: find and replace all instances of "SIDEseqRefSet" with "SIDEREF"

## TODO: default argument in SIDEREF for "R" once complete.

## TODO: simplify "setSimilarityMeasure" function.

## TODO: dist_and_clust_pipelines should be at least broken up into spectral embedding script and other helpers (and go through for no longer relevant functionalities

