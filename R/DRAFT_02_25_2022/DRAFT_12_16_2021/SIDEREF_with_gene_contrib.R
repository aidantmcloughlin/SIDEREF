SIDEseqRefSetGeneContrib <- function(expr_matrix, 
                                     diff_expr_method = "diff_expr_norm",
                                     similarity_method = "n_intersect",
                                     rbo_p = NULL,
                                     n_top_genes = 300,
                                     ## cell reference set selection parameters
                                     selection_method = "random",
                                     gamma = NULL,
                                     size_ref_set,
                                     n_clust = 10,
                                     B = 0,
                                     D = 1,
                                     purity_tol = 0.85, ## minimum purity score for algorithm to complete.
                                     n_neighbors = 15, 
                                     min_dist = 0.01,
                                     clust_repeats = 5,
                                     ## parallelization parameters
                                     parallelize = TRUE,
                                     n_cores = detectCores()-1,
                                     ## other
                                     n_top_contribs = NULL, ## set to the number of desired gene contributions, if any.
                                     verbose = TRUE) {
  
  
  if(!is.null(n_top_contribs) & D > 1) {stop("Gene Contribution Scores Only Supported for Single Dissim Matrix Computation")}
  
  ## set constants, parallel backend
  J <- dim(expr_matrix)[1]
  N <- dim(expr_matrix)[2]
  
  SIDESeq_method <- "ref_set"
  
  cell_pairs <- combn(N, 2)
  
  if(parallelize == TRUE) {
    cl <- makeCluster(n_cores) 
  }
  
  
  ## Set Differential Expression measure
  if(diff_expr_method == "diff_expr_norm") {
    diffExprMeasure <- function(i, j, expr_matrix) {
      return(abs(expr_matrix[, i] - expr_matrix[, j]) / 
               sqrt(expr_matrix[, i] + expr_matrix[, j]))
    }
  }
  
  ## Set Similarity measure
  similarityMeasure <- 
    setSimilarityMeasure(similarity_method = similarity_method, 
                         side_seq_method = SIDESeq_method,
                         rbo_p = rbo_p)
  
  
  
  ## Iterate procedure B times or until convergence criteria is reached
  b <- 1
  prior_mat <- NULL
  dissim_matrix <- NULL
  
  ## storing the convergence criterion
  purity_scores <- rep(0, B)
  
  ### Dissimilarity Matrix Computation Workflow:
  computeDissim <- function(contribs_val) {
    
    if(class(ref_cells) == "logical") {ref_cells <- which(ref_cells)}
    
    ## create ref set combn list
    ref_set_pairs <-
      matrix(c(rep(seq_len(N), length(ref_cells)),
               sapply(ref_cells, rep, N)), 
             byrow = TRUE,
             nrow  = 2)
    
    ## Compute Differential Expression Vector for each pair of cells 
    diff_expr_matrix <-
      computeDiffExprMat(cell_pairs = ref_set_pairs,
                         expr_matrix = expr_matrix, 
                         diff_expr_measure = diffExprMeasure,
                         n_top_genes = n_top_genes,
                         parallelize = parallelize,
                         cl = cl) 
    
    
    ## Compute similarity scores
    if(parallelize == TRUE) {
      clusterExport(cl, varlist=c("diff_expr_matrix", "N", "J",
                                  "ref_cells",
                                  "similarityMeasure", "rbo",
                                  "contribs_val"),
                    envir=environment())
      similarity_vec <-
        parallel::parApply(cl = cl, X = cell_pairs, 2, 
                           FUN = function(x) {
                             return(
                               similarityMeasure(i = x[1], j = x[2], N = N, J = J,
                                                 ref_set = ref_cells,
                                                 d_e_mat = diff_expr_matrix,
                                                 contribs_val = contribs_val)
                             )})
    } else{
      
      similarity_vec <-
        apply(cell_pairs, 2, 
              function(x) {
                similarityMeasure(i = x[1], j = x[2], N = N, J = J,
                                  ref_set = ref_cells,
                                  d_e_mat = diff_expr_matrix,
                                  contribs_val = contribs_val)})
      
    }
    
    ## separate gene contribution scores and similarity scores if needed.
    if(!is.null(contribs_val)) {
      gene_contrib_list <- lapply(similarity_vec, function(x) x[[2]])
      gene_contrib_vals <- lapply(similarity_vec, function(x) x[[3]])
      
      similarity_vec <- sapply(similarity_vec, function(x) x[[1]])
    }
    
    ## dissimilarity measure is max score - similarity_score
    dissimilarity_vec <- max(similarity_vec) - similarity_vec
    
    ## convert vector to cell-by-cell matrix
    dissim_matrix <- matrix(0, nrow = N, ncol = N)
    
    dissim_matrix[lower.tri(dissim_matrix)] <- dissimilarity_vec
    
    
    dissim_matrix[upper.tri(dissim_matrix)] <- 
      t(dissim_matrix)[upper.tri(dissim_matrix)]
    
    if(!is.null(contribs_val)) {
      return(list(dissim_matrix = dissim_matrix,
                  gene_contrib_list = gene_contrib_list,
                  gene_contrib_vals = gene_contrib_vals))
    } else{return(list(dissim_matrix = dissim_matrix))}
    
  }
  
  ### Run cluster convergence loop
  while(b <= B & purity_scores[max(1, b-1)] < purity_tol) {
    if(verbose == TRUE) {cat("Running convergence iteration", b, "of", B, "...")}
    ## extract reference set
    ref_cells <-
      selectRefSet(expr_matrix,
                   selection_method = selection_method,
                   size_ref_set = size_ref_set,
                   gamma = gamma,
                   dissim_matrix = dissim_matrix,
                   n_clust,
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
    
    ## compute dissim matrix
    dissim_matrix <- computeDissim(contribs_val = NULL ## don't need gene contributions during internal runs. 
                                   )$dissim_matrix
    
    ## compare to prior matrix
    if(!is.null(prior_mat)) {
      
      internal_purity_vec <- c()
      for(r in seq_len(clust_repeats)) {
        ## UMAP embedding + Kmeans
        umap_embed_current <-
          uwot::umap(stats::as.dist(dissim_matrix),
                     n_neighbors = n_neighbors,
                     min_dist = min_dist)
        
        umap_embed_prior <- 
          uwot::umap(stats::as.dist(prior_mat),
                     n_neighbors = n_neighbors,
                     min_dist = min_dist)
        
        clust_res_current <- kmeanspp(umap_embed_current, k = n_clust)$cluster
        clust_res_prior   <- kmeanspp(umap_embed_prior,   k = n_clust)$cluster
        
        ## storing purity of current cluster with previous clustering:
        internal_purity_vec[r] <-
          clusterPurity(clust_res_current, 
                        clust_res_prior) 
      }
      purity_scores[b] <- mean(internal_purity_vec)
      if(verbose == TRUE){cat("Current clustering purity score", purity_scores[b])}
    }
    
    ## update prior
    prior_mat <- dissim_matrix
    
    ## end of while loop
    b <- b + 1
  }
  
  
  dissim_final <- NULL
  
  ## Run final averaging loop
  ref_cell_list <- 
    lapply(seq_len(D),
           function(d) {
             return(selectRefSet(expr_matrix,
                                 selection_method = selection_method,
                                 size_ref_set = size_ref_set,
                                 gamma = gamma,
                                 dissim_matrix = dissim_matrix,
                                 n_clust,
                                 n_neighbors = n_neighbors,
                                 min_dist = min_dist))
           })
  
  
  
  ## clear space
  rm(prior_mat,dissim_matrix)
  
  for(d in seq_len(D)) {
    if(verbose == TRUE) {cat("Running averaging iteration", d, "of", D, "...")}
    
    ref_cells <- ref_cell_list[[d]]
    
    dissim_curr <- computeDissim(contribs_val = n_top_contribs)
    
    
    if(is.null(dissim_final)){dissim_final <- dissim_curr$dissim_matrix}
    ## incorporate new matrix into running average.
    dissim_final <-
      1/(d+1) * (d * dissim_final + dissim_curr$dissim_matrix)
    
  }
  
  
  ## close workers if needed
  if(parallelize == TRUE) {
    stopCluster(cl)
  }
  
  if(!is.null(n_top_contribs)) {
    return(list(dissim_final      = dissim_final, 
                purity_scores     = purity_scores,
                gene_contrib_list = dissim_curr$gene_contrib_list,
                gene_contrib_vals = dissim_curr$gene_contrib_vals))
  } else{return(list(dissim_final      = dissim_final, 
                     purity_scores     = purity_scores))}
  
}

