## TODO: change file name

evaluateSIDEseqRun <- 
  function(
    ## Parameters specific to SIDEseq evaluation:
    cell_types = NULL, ## vector of the cell labels.
    n_runs = 1,
    n_spect_clusts = 50,
    ## Parameters for SIDEseqRefSet run:
    expr_matrix, 
    diff_expr_method,
    similarity_method = "top_n_gene_intersect",
    n_top_genes = 300,
    selection_method,
    size_ref_set,
    n_clust = 10,
    B = 1,
    prior_mat_method,
    parallelize = TRUE,
    n_cores = detectCores()-1,
    verbose = TRUE) {
    
    ## TODO: no need to specify n_clust if cell_types arg is provided.
    
    results_df <- data.frame(diff_expr_method = diff_expr_method,
                             similarity_method = similarity_method,
                             n_top_genes = n_top_genes,
                             selection_method = selection_method,
                             size_ref_set = size_ref_set,
                             n_clust = n_clust,
                             B = B,
                             prior_mat_method,
                             parallelize = parallelize,
                             n_cores = n_cores,
                             ## results variables
                             run_time = rep(NA, n_runs),
                             avg_sil_width = rep(NA, n_runs),
                             rand_index = rep(NA, n_runs))
    
    
    
    for(r in seq_len(n_runs)) {
      
      ## tracking full run time.
      ptm <- proc.time()
      side_seq_res <-
        SIDEseqRefSet(expr_matrix = expr_matrix, 
                      diff_expr_method = diff_expr_method,
                      similarity_method = similarity_method,
                      n_top_genes = n_top_genes,
                      ## cell reference set selection parameters
                      selection_method = selection_method,
                      size_ref_set = size_ref_set,
                      n_clust = n_clust,
                      B = B,
                      prior_mat_method = prior_mat_method, ## c("new", "average")
                      ## parallelization parameters
                      parallelize = parallelize,
                      n_cores = n_cores,
                      ## other
                      verbose = verbose)
      
      ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
      
      results_df$run_time[r] <- ct_elap
      
      ## compute average silhouette width
      
      if(!is.null(cell_types)) {
        n_clust <- length(unique(cell_types))
      }
      
      
      avg_sil_width_vec <- rep(NA, n_spect_clusts)
      adj_rand_ind_vec  <- rep(NA, n_spect_clusts)
      
      for(s in seq_len(n_spect_clusts)) {
        spect_clust_res <- spectralCluster(side_seq_res$sim_matrix,
                                           n_clust = n_clust)
        
        paired_dists   <- stats::dist(spect_clust_res$G)
        
        ## silhouette scores
        sil_scores <- 
          cluster::silhouette(spect_clust_res$cluster, 
                              paired_dists)[,3]
        
        avg_sil_width_vec[s] <- mean(sil_scores)
        
        ## adjust rand index if test clusters are given.
        if(!is.null(cell_types)) {
          adj_rand_ind_vec[s] <- 
            mclust::adjustedRandIndex(cell_types, spect_clust_res$clusters)
        }
        
      }
      
      ## avg results across cluster attempts:
      results_df$avg_sil_width[r] <- mean(avg_sil_width_vec)
      if(!is.null(cell_types)) {
        results_df$rand_index[r] <- mean(adj_rand_ind_vec)
      }
      
    }
    
    eval_vars <- c("run_time", "avg_sil_width", "rand_index")
    
    ## average results 
    results_df <-
      results_df %>%
      group_by_at(setdiff(names(results_df), eval_vars)) %>%
      summarise_at(all_of(eval_vars), mean) %>%
      ungroup()
      
    
    return(results_df)
    
    
}


