set.seed(1)
library(here)
source(here("R/libraries.R"))
source(here("R/converge_analysis_helpers.R"))


compute_time_SIDEREF_SIDEseq <- 
  function(assay_data, 
           meta_data,
           n_cells_sizes,
           n_reps = 1,
           g = G,
           r = R,
           n_pcs = N_PCS,
           n_clust = N_CLUST,
           min_dist = MIN_DIST,
           n_neighbors = N_NEIGHBORS,
           n_cores = COMPUTE_TIME_N_CORES,
           save_intermed_file = NULL) {
    
    compute_time_df <- 
      data_frame(n_cells = n_cells_sizes,
                 SIDEREF_comp_time = 0,
                 SIDEseq_comp_time = 0)
    
    for(n in seq_len(length(n_cells_sizes))) {
      cat("computing cell size: " %p% n_cells_sizes[n])
      ## initialize
      sideseq_difftime <- 0
      sideref_difftime <- 0
      for(s in seq_len(n_reps)) {
        
        ## get data sample
        samp_inds <- sample(ncol(assay_data),
                            n_cells_sizes[n])
        
        assay_sample <- assay_data[, samp_inds]
        
        ## COMPUTE SIDEseq
        start_time_sideseq <- Sys.time()
        SIDEseq_res <- 
          SIDEseq(assay_sample,
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = g,
                  parallelize = TRUE,
                  n_cores = n_cores)
        end_time_sideseq <- Sys.time()
        
        sideseq_difftime <- sideseq_difftime + 
          difftime(end_time_sideseq, start_time_sideseq, units='mins')
        
        start_time_sideref <- Sys.time()
        SIDEREF_res_strat <-
          SIDEREF(assay_sample,
                  n_top_genes = g,
                  R = r,
                  selection_method = "cell_embed_sample",
                  n_pcs = n_pcs,
                  n_clust = n_clust,
                  min_dist = MIN_DIST,
                  n_neighbors = N_NEIGHBORS,
                  parallelize = TRUE,
                  n_cores = n_cores,
                  verbose = TRUE)
        end_time_sideref <- Sys.time()
        
        sideref_difftime <- sideref_difftime + 
          difftime(end_time_sideref, start_time_sideref, units='mins')
        
      } ## end n_reps loop  
      ## average and append
      sideseq_difftime <- sideseq_difftime / n_reps
      sideref_difftime <- sideref_difftime / n_reps
      
      compute_time_df$SIDEseq_comp_time[n] <- sideseq_difftime
      compute_time_df$SIDEREF_comp_time[n] <- sideref_difftime
      
      if(!is.null(save_intermed_file)) {
        save(compute_time_df, file = save_intermed_file)
      }
      
    } ## end n_cells_sizes loop
    
    return(compute_time_df)
    
  }






###############################################################################
### MAIN: DATA LOAD AND RUN.
###############################################################################

###########################################
### Load and prep TM data.
source(here("R/tm_helper_functions.R"))


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

### Run analysis function with loaded data:

compute_time_SIDEREF_SIDEseq(
  assay_data, 
  meta_data,
  n_cells_sizes = N_CELLS_SIZES,
  n_reps = COMPUTE_TIME_REPS,
  g = G,
  r = R,
  n_pcs = N_PCS,
  n_clust = N_CLUST,
  min_dist = MIN_DIST,
  n_neighbors = N_NEIGHBORS,
  n_cores = COMPUTE_TIME_N_CORES,
  save_intermed_file = here("output/computations/" %p%
                              "compute_time_SIDEREF_SIDEseq.RData"))




