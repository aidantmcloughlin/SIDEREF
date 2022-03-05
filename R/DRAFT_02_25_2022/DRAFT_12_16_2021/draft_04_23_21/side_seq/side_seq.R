
 
SIDEseqSimult <- function(expr_matrix, 
                          diff_expr_method = "diff_expr_norm",
                          similarity_method = "n_intersect",
                          n_top_genes = 150,
                          parallelize = TRUE,
                          n_cores = detectCores()-1) {
  
  ## set constants, parallel backend
  J <- dim(expr_matrix)[1]
  N <- dim(expr_matrix)[2]
  
  SIDESeq_method <- "simultaneous"
  
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
  
  
  
  ## Compute Differential Expression Vector for each pair of cells 
  diff_expr_matrix <-
    computeDiffExprMat(cell_pairs, expr_matrix, 
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
    dissimilarity_vec <-
      parApply(cl = cl, X = cell_pairs, 2, 
               FUN = function(x) {
                 return(
                   similarityMeasure(i = x[1], j = x[2], N = N,
                                     d_e_mat = diff_expr_matrix)
                 )})
  } else{
    dissimilarity_vec <-
      apply(cell_pairs, 2, 
            function(x) {
              similarityMeasure(i = x[1], j = x[2], N = N, 
                                d_e_mat = diff_expr_matrix)})
    
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






