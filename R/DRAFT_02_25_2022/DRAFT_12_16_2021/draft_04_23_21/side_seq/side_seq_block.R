
## Get number of cells per block
computeBlockSize <- function(ncols,
                             nrows,
                             block_memory,
                             entry_size = 4.001) {
  ## Block memory is given in MegaBytes
  ## matrix of positive int at 4.001 bytes per int
  entries = ncols*nrows
  
  ## test is full block fits in memory limit
  if(entries*entry_size*1e-6 <= block_memory) {
    block_cols = ncols
    n_blocks = 1
  } else{## find max ncols per memory limit
    max_block_cols = (block_memory * 1e6) / (nrows * entry_size)
    
    ## number of possible dissim measures per block
    k_poss = floor(0.5*(2*N - 1 - sqrt((2*N-1)^2 - 8 * max_block_cols)))
    
    ## smooth block size across number of blocks
    n_blocks = ceiling(N / k_poss)
    
    ## TODO: check whether this needed
    k_per_block = min(k_poss, ceiling(N / n_blocks))
  }
  
  return(k_per_block)
}


## Compute as many dissimilarity scores at a time as memory constraints allow 
SIDEseqBlockSimult <- function(expr_matrix, 
                               diff_expr_method = "diff_expr_norm",
                               similarity_method = "top_n_gene_intersect",
                               n_top_genes = 1000,
                               block_memory = 3e3,
                               parallelize = TRUE,
                               ## parallelization options
                               max_seed = 1e9,
                               n_cores = detectCores()-1) {
  ## Block memory is given in MegaBytes
  ## set constants, parallel backend
  J <- dim(expr_matrix)[1]
  N <- dim(expr_matrix)[2]
  
  SIDESeq_method <- "block"
  
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
                         side_seq_method   = SIDESeq_method)
  
  ## Get number of cols in each calculation block
  n_pairs <- dim(cell_pairs)[2]
  k_per_block <- computeBlockSize(ncols = n_pairs,
                                  nrows = n_top_genes,
                                  block_memory = block_memory)
  
  ## initialize dissimilarity vector
  dissimilarity_vec <- rep(0, n_pairs)
  
  ## export needed parallelization objects
  if(parallelize == TRUE) {
    clusterExport(cl, varlist = c("diffExprMeasure", "expr_matrix",
                                  "n_top_genes", "block_list"))
  }
  ## initialize
  lower_ind <- 0
  upper_ind <- 0
  ## iterate over each block
  while(upper_ind < N) {
    cat("Cells", lower, "to", block_list[[2]])
    col_indices <- 
      (block_list[[1]] * (b-1)+1):
      (min(n_pairs, block_list[[1]] * b))
    
    cell_pairs_subset <- cell_pairs[, col_indices]
    ## Compute Differential Expression Vector for each pair of cells
    diff_expr_matrix <-
      computeDiffExprMat(cell_pairs_subset, expr_matrix, 
                         diff_expr_measure = diffExprMeasure,
                         n_top_genes = n_top_genes,
                         parallelize = parallelize,
                         cl = cl)
    
    ## Compute dissimilarity scores
    
    ## First, get similarity scores
    if(parallelize == TRUE) {
      
      clusterExport(cl, varlist=c("diff_expr_matrix", "N",
                                  "similarityMeasure", "getCombnPairIndex"),
                    envir=environment())
      dissimilarity_vec[col_indices] <-
        parApply(cl = cl, X = cell_pairs_subset, 2, 
                 FUN = function(x) {
                   return(
                     similarityMeasure(i = x[1], j = x[2], N = N,
                                       block_cols = block_list[[1]],
                                       d_e_mat = diff_expr_matrix)
                   )})
    } else{
      dissimilarity_vec[col_indices] <-
        apply(cell_pairs_subset, 2, 
              function(x) {
                print(x)
                similarityMeasure(i = x[1], j = x[2], N = N, 
                                  block_cols = block_list[[1]],
                                  d_e_mat = diff_expr_matrix)})
      
    }
    
    
  }
  
  ## dissimilarity measure is max score - similarity_score
  dissimilarity_vec <- max(dissimilarity_vec) - dissimilarity_vec
  
  ## convert vector to cell-by-cell matrix
  D <- matrix(0, nrow = N, ncol = N)
  
  D[lower.tri(D)] <- dissimilarity_vec
  
  
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  
  ## close workers if needed
  if(parallelize == TRUE) {
    stopCluster(cl)
  }
  
  return(D)
  
}

