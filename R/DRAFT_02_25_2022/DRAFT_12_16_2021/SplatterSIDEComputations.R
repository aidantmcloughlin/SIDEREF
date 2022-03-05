###############################################################################
### Primary Splatter Computations
###############################################################################

set.seed(4797827)


def_n_top_genes <- 150



## load simulation data
load(here("output/splatter_sideseq/splatter_sim1.RData"))
## Extract log counts
sim1 <- logNormCounts(sim1)
sim1_counts <- as.matrix(sim1@assays@data$logcounts)

## Different Reference Set Sizes of Interest
ref_set_sizes <- c(5, 10, 25, 50, 75, 100, 150) 


## Loop of reference set simulations  
for(s in 1:5) {
  for(r in ref_set_sizes) {
    print(r)
    side_seq_ref_rand <- 
      SIDEseqRefSet(expr_matrix = sim1_counts, 
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = def_n_top_genes,
                    ## cell reference set selection parameters
                    selection_method = "random",
                    size_ref_set = r,
                    n_clust = 5,
                    B = 0,
                    D = 1,
                    ## parallelization parameters
                    parallelize = TRUE,
                    n_cores = detectCores()-2,
                    ## other
                    verbose = TRUE)
    
    
    save(side_seq_ref_rand, 
         file = here(paste0("output/splatter_sideseq/run", s, 
                            "_", r, "rand_ref_set_splatter.RData")))
    
    side_seq_ref_cell_embed <- 
      SIDEseqRefSet(expr_matrix = sim1_counts, 
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = def_n_top_genes,
                    ## cell reference set selection parameters
                    selection_method = "cell_embed_sample",
                    size_ref_set = r,
                    n_clust = 5,
                    B = 3,
                    D = 1,
                    ## parallelization parameters
                    parallelize = TRUE,
                    n_cores = detectCores()-2,
                    ## other
                    verbose = TRUE)
    
    
    save(side_seq_ref_cell_embed, 
         file = here(paste0("output/splatter_sideseq/run", s,
                            "_", r, "cell_embed_ref_set_splatter.RData")))
    
    if(r >= 10) {
      side_seq_ref_cell_embed <- 
        SIDEseqRefSet(expr_matrix = sim1_counts, 
                      diff_expr_method = "diff_expr_norm",
                      similarity_method = "n_intersect",
                      n_top_genes = def_n_top_genes,
                      ## cell reference set selection parameters
                      selection_method = "cell_embed_sample",
                      size_ref_set = r,
                      n_clust = 10,
                      B = 3,
                      D = 1,
                      ## parallelization parameters
                      parallelize = TRUE,
                      n_cores = detectCores()-2,
                      ## other
                      verbose = TRUE)
      
      
      save(side_seq_ref_cell_embed, 
           file = here(paste0("output/splatter_sideseq/run", s,
                              "_", r, "cell_embed_ref_set_splatter_10_clust.RData"))) 
    }
    
  }
  
  
}
