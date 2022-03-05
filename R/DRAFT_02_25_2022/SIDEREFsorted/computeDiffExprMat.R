computeDiffExprMat <- function(cell_pairs, expr_matrix, diff_expr_measure,
                               n_top_genes,
                               parallelize = TRUE,
                               cl = NULL,
                               max_seed = 1e9,
                               do_sort = FALSE) {
  ## cell pairs: 2xP matrix of cell pairings to compute diff expr between
  
  if(parallelize == TRUE) {
    clusterExport(cl, varlist=c("cell_pairs", "expr_matrix",
                                "diff_expr_measure", "n_top_genes",
                                "do_sort"),
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
                 if(do_sort) {
                   top_genes <-
                     sort(
                       order(-diff_expr, sample(seq_len(length(diff_expr)))
                       )[1:n_top_genes])
                 } else{
                   top_genes <- 
                     order(-diff_expr, sample(seq_len(length(diff_expr)))
                     )[1:n_top_genes]
                 }
                 return(top_genes)
               })
    
  } else{## no parallel compute:
    ## Compute differential expression measure for each gene
    diff_expr_matrix <-
      apply(cell_pairs, 2, 
            function(x) {
              diff_expr = diff_expr_measure(x[1], x[2], expr_matrix)
              ## Take top n diff expr genes for each column
              if(do_sort) {
                top_genes <-
                  sort(
                  order(-diff_expr, sample(seq_len(length(diff_expr)))
                  )[1:n_top_genes])
              } else{
                top_genes <- 
                  order(-diff_expr, sample(seq_len(length(diff_expr)))
                  )[1:n_top_genes]
              }
              return(top_genes)
            })
  }
  
  return(diff_expr_matrix)
}


