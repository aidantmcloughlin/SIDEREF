SIDEseqIterated <- function(expr_matrix, 
                            diff_expr_method = "diff_expr_norm",
                            similarity_method = "top_n_gene_intersect",
                            n_top_genes = 150,
                            parallelize = TRUE,
                            n_cores = detectCores()-1) {
  
  ## set constants, parallel backend
  J <- dim(expr_matrix)[1]
  N <- dim(expr_matrix)[2]
  
  SIDESeq_method <- "iterative"
  
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
                         side_seq_method = SIDESeq_method)
  
  ## get similarity val function:
  getSimilarityVal <- function(x) {
    ## Relevant cell pairings for each S_{i,j} calculation:
    cell_pairs_subset <- matrix(c(rep(x[1],(N-2)), rep(x[2],(N-2)), 
                                  rep((1:N)[c(-x[1],-x[2])],2)), 
                                nrow=2, byrow=TRUE)
    
    ## Compute Differential Expression Vector for each pair of cells 
    diff_expr_matrix <-
      computeDiffExprMat(cell_pairs_subset, expr_matrix, 
                         diff_expr_measure = diffExprMeasure,
                         n_top_genes = n_top_genes,
                         parallelize = FALSE)
    
    return(similarityMeasure(i = x[1], j = x[2], N = N, 
                             d_e_mat = diff_expr_matrix))
  }
  
  
  
  ## initialize dissimilarity vector
  dissimilarity_vec <- rep(0, dim(cell_pairs)[2])
  
  ## idea: run entire pipeline to dissimilarity vector in order to limit the 
  ##    memory requirements storing full diff-exp matrix.
  if(parallelize == TRUE) {
      clusterExport(cl, varlist = c("N", "diffExprMeasure", "expr_matrix", 
                                    "n_top_genes", "parallelize",
                                    "computeDiffExprMat",
                                    "similarityMeasure", "getCombnPairIndex",
                                    "getSimilarityVal",
                                    "clusterExport"),
                    envir = environment())
      
      dissimilarity_vec <-
        parApply(cl, cell_pairs, 2, getSimilarityVal)
      
  } else{
    dissimilarity_vec <-
      apply(cell_pairs, 2, getSimilarityVal)
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
  
  
  
  