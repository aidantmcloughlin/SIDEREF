### prelim gene contrib Splatter results


## common paste helper function:
`%p%` <- function(x,y) {paste0(x,y)}


n_neighbors <- 15
min_dist <- 0.01

## Load Functions:
source(here('R/computeDiffExprMat.R'))
source(here('R/setSimilarityMeasure_new.R'))
source(here('R/SIDEseq.R'))
source(here('R/SIDEREF_with_gene_contrib.R'))
source(here('R/PlotFunctions.R'))

## load simulation data
load(here("output/splatter_sim2.RData"))

sim2 <- logNormCounts(sim2)
sim2_counts <- as.matrix(sim2@assays@data$logcounts)
## remove genes with no expression
sim2_counts <- sim2_counts[rowSums(sim2_counts) > 0, ]

splat_group_labels <- factor(sim2@colData@listData$Group)


levels(splat_group_labels) <- 
  c("1: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 2, 3)", 
    "2: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 1, 3)",
    "3: D.E. Prob 3.0%\nShared D.E. Prob 2.0% (with Grp. 1, 2)",
    "4: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 5)",
    "5: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 4)",
    "6: D.E. Prob 4.0%\nNo Shared D.E.")


## run gene contrib SIDEREF.

side_seq_ref_rand <- 
  SIDEseqRefSetGeneContrib(expr_matrix = sim2_counts, 
                           diff_expr_method = "diff_expr_norm",
                           similarity_method = "n_intersect",
                           n_top_genes = 50,
                           ## cell reference set selection parameters
                           selection_method = "random",
                           size_ref_set = r,
                           n_clust = 5,
                           B = 0,
                           D = 1,
                           ## parallelization parameters
                           parallelize = TRUE,
                           n_cores = 5,
                           ## other
                           n_top_contribs = 15,
                           verbose = TRUE)

## TODO: convert this to a function.
group_2_ind <- which(factor(sim2@colData@listData$Group) == "Group2")

N <- dim(sim2_counts)[2]

index_mat <- matrix(0, nrow = N, ncol=N)

index_mat[lower.tri(index_mat)] <- 1:length(side_seq_ref_rand$gene_contrib_list)

group_2_list_inds <- c(index_mat[group_2_ind, ])

group_2_list_inds <- group_2_list_inds[group_2_list_inds != 0]

group_2_gene_lists <- side_seq_ref_rand$gene_contrib_list[group_2_list_inds]
                             
group_2_gene_scores <- side_seq_ref_rand$gene_contrib_vals[group_2_list_inds]


full_gene_list <- rep(0, dim(sim2_counts)[1])

n_comparisons <- length(group_2_gene_lists)

for(i in seq_len(n_comparisons)) {
  genes <- group_2_gene_lists[[i]]
  
  full_gene_list[genes] <- full_gene_list[genes] + group_2_gene_scores[[i]] / n_comparisons
  
}

group2_gene_list <- order(-full_gene_list, sample(seq_len(length(full_gene_list))))[1:15]

group2_gene_scores <- full_gene_list[group2_gene_list]




