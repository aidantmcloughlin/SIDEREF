### SIDEseq simulation examples -----------------------------------------------


### examine distribution of normalized counts for highly expressed genes

mean_top_expr <-
  apply(umi_counts_norm, 2, function(x) -mean(sort(-x)[1:10]))




## Running DE analyses for subset of cell types
table(celltype_df$celltype)

sub_list <- c("A1.E1", "EN06", "IN03", "MX", "P1")


celltype_df_sub <- celltype_df %>% 
  dplyr::filter(celltype %in% sub_list)


## number of expressed genes per cell type
n_expr_genes <- 
  apply(umi_counts_norm, 2, function(x) {
      return(sum(x > log(2)))
    })

### top n expressed genes of each cell type
n_top_genes <- 50

top_n_genes <- 
  apply(umi_counts_norm, 2, function(x) {
    which(data.table::frank(-x, na.last = TRUE,
                            ties.method = "random") <= n_top_genes)
  })

  
expr_matrix_subset <- umi_counts_norm[,celltype_df_sub$barcode]

cell_pairs <- combn(ncol(expr_matrix_subset), 2)


## subset to highly variable genes
var_genes <- getHighlyVariableGenes(expr_matrix_subset,
                                    p_val_thres = 1e-10)

## TODO: GiniClust
expr_matrix_subset <- 
  expr_matrix_subset[var_genes,]


diffExprMeasure <- function(i, j, expr_matrix) {
  return(abs(expr_matrix[, i] - expr_matrix[, j]))
}

cl <- makeCluster(6)
clusterExport(cl, varlist=c("expr_matrix_subset",
                            "diffExprMeasure"),
              envir=environment())
## Compute top 20 diff expr measures for each genes


diff_expr_mat <- 
  parApply(cl = cl, X = cell_pairs, 2,
           FUN = function(x) {
             ## set seed for reproducibility
             diff_expr = diffExprMeasure(x[1], x[2], expr_matrix_subset)
             
             return(-sort(-diff_expr)[1:150])
           })



stopCluster(cl)

## TODO: sorting implementation of differential expression matrix computation comparison
## TODO: utilizing gene ontologies to develop separations.

SIDEseq_subset_res <-
  SIDEseqSimult(expr_matrix = expr_matrix_subset, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "top_n_gene_intersect",
                n_top_genes = 150,
                parallelize = TRUE)


## group by cell type compute DE intersection
cell_pairs_indices <- vector(mode = "list", length = length(sub_list))
diffcell_pairs_indices <- vector(mode = "list", length = length(sub_list))

names(cell_pairs_indices) <- sub_list
names(diffcell_pairs_indices) <- sub_list

for(cg in sub_list) {
  print(cg)
  cell_indices <- which(celltype_df_sub$celltype == cg)
  other_indices <- setdiff(seq_len(nrow(celltype_df_sub)), cell_indices)
  
  cell_pairs_indices[[cg]] <- 
    data.frame(cg    = paste0(cg, "same"),
               index = which(apply(cell_pairs, 2, function(x) {
    x[1] %in% cell_indices & x[2] %in% cell_indices})))
  
  diffcell_pairs_indices[[cg]] <- 
    data.frame(cg    = paste0(cg, "diff"),
               index = which(apply(cell_pairs, 2, function(x) {
    (x[1] %in% cell_indices & x[2] %in% other_indices) | 
      (x[1] %in% other_indices & x[2] %in% cell_indices)})))
  
}





cell_diff_expr_analysis <-
  celltype_df %>% 
  mutate(mean_top_expr = mean_top_expr,
         n_expr_genes  = n_expr_genes,
         top_n_genes   = top_n_genes)

## TODO: Frequency histogram
cell_diff_expr_analysis %>% 
  ggplot(aes(x=mean_top_expr)) + geom_histogram() + facet_wrap(~celltype)


cell_diff_expr_analysis %>% 
  ggplot(aes(x=mean_top_expr)) + geom_density() + facet_wrap(~celltype)


cell_diff_expr_analysis %>% 
  ggplot(aes(x=n_expr_genes)) + geom_density() + facet_wrap(~celltype)


cell_diff_expr_analysis %>% 
  ggplot(aes(x=n_expr_genes, color=celltype)) + geom_density()





### SIDEseq OG simulation -----------------------------------------------------

#small_scrna <- 








     