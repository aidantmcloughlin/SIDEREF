###############################################################################
### Compare full side seq to SIDEREF runs.
###############################################################################

library(here)
source(here("R/libraries.R"))
source(here("R/converge_analysis_helpers.R"))
set.seed(1)




snn_SIDEREF_SIDEseq <- 
  function(assay_data, 
           meta_data,
           ref_set_sizes,
           sample_size,
           n_reps,
           g = G,
           n_pcs = N_PCS,
           n_clust = N_CLUST,
           min_dist = MIN_DIST,
           n_neighbors = N_NEIGHBORS,
           snn_k = SNN_K,
           n_cores = N_CORES,
           save_intermed_file = NULL) {
  
    
    ## Sample the data:
    ## get data sample
    samp_inds <- sample(ncol(assay_data),
                        sample_size)
    
    assay_sample <- assay_data[, samp_inds]
    
    ## COMPUTE SIDEseq
    SIDEseq_res <- 
      SIDEseq(assay_sample, 
              diff_expr_method = "diff_expr_norm",
              similarity_method = "n_intersect",
              n_top_genes = g,
              parallelize = TRUE,
              n_cores = n_cores)
    
    knn_SIDEseq <- getKNNList(SIDEseq_res, k = snn_k)    
    
    embed_res <- list()
    rand_res <- list()
    
    embed_table     <- list()
    rand_table      <- list()
    
    methods <- c("Cell Type Stratified Sample",
                 "Random Sample")
    
    nn_consistency <- 
      as.data.frame(
        expand.grid(method = methods, 
                    ref_set_size = ref_set_sizes),
        stringsAsFactors = FALSE)
    
    for(n in seq_len(n_reps)) {
      nn_consistency['avg_snn_run' %p% n] <- numeric()
    } 
    for(n in seq_len(n_reps)) {
      nn_consistency['compute_time_run' %p% n] <- numeric()
    }
    
    ## SIDEREF COMPUTATION LOOP
    i <- 1
    for(r in ref_set_sizes) {
      if(!is.null(save_intermed_file)) {
        save(nn_consistency, file = save_intermed_file)
      }
      
      cat("Ref Set Size " %p% r)
      for(s in 1:n_reps) {
        cat("Rep " %p%s)
        ## stratified sample
        start_time_strat <- Sys.time()
        SIDEREF_res_strat <-
          SIDEREF(expr_matrix = assay_sample,
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
        end_time_strat <- Sys.time()
        
        ## random sample
        start_time_rand <- Sys.time()
        SIDEREF_res_rand <-
          SIDEREF(expr_matrix = assay_sample,
                  n_top_genes = g,
                  R = r,
                  selection_method = "random",
                  parallelize = TRUE,
                  n_cores = n_cores,
                  verbose = TRUE)
        end_time_rand <- Sys.time()
        
        knn_SIDEREF_embed <- 
          getKNNList(SIDEREF_res_strat, k = snn_k)
        knn_SIDEREF_rand <- 
          getKNNList(SIDEREF_res_rand,  k = snn_k)
        
        snn_embed <- 
          sapply(seq_len(nrow(SIDEseq_res)), 
                 function(i){length(intersect(knn_SIDEREF_embed[[i]],
                                              knn_SIDEseq[[i]]))})
        snn_rand <- 
          sapply(seq_len(nrow(SIDEseq_res)), 
                 function(i){length(intersect(knn_SIDEREF_rand[[i]],
                                              knn_SIDEseq[[i]]))})
        
        avg_embed     <- mean(snn_embed)
        avg_rand      <- mean(snn_rand)
        
        nn_consistency['avg_snn_run' %p% s][i,] <- 
          avg_embed
        nn_consistency['compute_time_run' %p% s][i,] <- 
          difftime(end_time_strat, start_time_strat, units='secs')
        
        
        nn_consistency['avg_snn_run' %p% s][i+1,] <- 
          avg_rand
        nn_consistency['compute_time_run' %p% s][i+1,] <- 
          difftime(end_time_rand, start_time_rand, units='secs')
        
      }
      i <- i + 2
    }
    
    ## summary statistics
    nn_consistency$avg_snn <- rep(0, nrow(nn_consistency))
    nn_consistency$max_snn <- rep(0, nrow(nn_consistency))
    nn_consistency$min_snn <- rep(snn_k, nrow(nn_consistency))
    nn_consistency$avg_compute_time <- rep(0, nrow(nn_consistency))
    
    for(s in seq_len(n_reps)) {
      nn_consistency$avg_snn <- nn_consistency$avg_snn + 
        nn_consistency['avg_snn_run' %p% s][,1]
      
      nn_consistency$max_snn <- 
        pmax(nn_consistency$max_snn, 
             nn_consistency['avg_snn_run' %p% s][,1])
      
      nn_consistency$min_snn <- 
        pmin(nn_consistency$min_snn, 
             nn_consistency['avg_snn_run' %p% s][,1])
      
      nn_consistency$avg_compute_time <- 
        nn_consistency$avg_compute_time + 
        nn_consistency['compute_time_run' %p% s][,1]
    }
    nn_consistency$avg_snn <- nn_consistency$avg_snn / n_reps
    nn_consistency$avg_compute_time <- nn_consistency$avg_compute_time / n_reps
    
      
    if(!is.null(save_intermed_file)) {
      save(nn_consistency, file = save_intermed_file)
    }
    return(nn_consistency)
    
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

snn_SIDEREF_SIDEseq(
  assay_data, 
  meta_data,
  ref_set_sizes = REF_SET_SIZES,
  sample_size = SNN_SAMP_SIZE,
  n_reps = SNN_REPS,
  save_intermed_file = here("output/computations/" %p% 
                              "snn_SIDEREF_SIDEseq.RData"))


