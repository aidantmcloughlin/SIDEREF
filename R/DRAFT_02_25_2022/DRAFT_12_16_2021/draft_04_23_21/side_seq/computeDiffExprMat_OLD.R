computeDiffExprMat <- function(cell_pairs, expr_matrix, diff_expr_measure,
                               n_top_genes,
                               parallelize = TRUE,
                               cl = NULL,
                               max_seed = 1e9) {
  ## cell pairs: 2xP matrix of cell pairings to compute diff expr between
  
  if(parallelize == TRUE) {
    clusterExport(cl, varlist=c("cell_pairs", "expr_matrix",
                                "diff_expr_measure", "n_top_genes"),
                  envir=environment())}
  
  
  ## computations:
  if(parallelize == TRUE) {
    ## need to establish seeds for reproducibility
    cell_pairs <- rbind(
      cell_pairs,
      sample(max_seed, size = dim(cell_pairs)[2], replace = TRUE))
    
    ## Compute differential expression measure for each gene, then rank
    diff_expr_matrix <- 
      parApply(cl = cl, X = cell_pairs, 2,
               FUN = function(x) {
                 ## set seed for reproducibility
                 set.seed(x[3])
                 diff_expr = diff_expr_measure(x[1], x[2], expr_matrix)
                 ## Take top n diff expr genes for each column
                 top_genes <-
                   which(data.table::frank(-diff_expr, na.last = TRUE,
                                           ties.method = "random") <= 
                           n_top_genes)
                 return(top_genes)
               })
    
  } else{## no parallel compute:
    ## Compute differential expression measure for each gene
    diff_expr_matrix <-
      apply(cell_pairs, 2, 
            function(x) {
              diff_expr = diff_expr_measure(x[1], x[2], expr_matrix)
              ## Take top n diff expr genes for each column
              top_genes <-
                which(data.table::frank(-diff_expr, na.last = TRUE,
                                        ties.method = "random") <= 
                        n_top_genes)
              return(top_genes)
            })
    
  }
  return(diff_expr_matrix)
}

