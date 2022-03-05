###############################################################################
###
### SIDEseq computations on SNARE RNA data.
###
###############################################################################


set.seed(5673290)

## Load Training and Test Splits Data

load(here("output/snare_sample_1500.RData"))
load(here("output/snare_test.RData"))

DefaultAssay(snare_sample) <- "RNA"
DefaultAssay(snare_test)   <- "RNA"

snare_sample <- NormalizeData(snare_sample)
## This effectively just extracts normalized counts for variable genes
snare_sample <- ScaleData(snare_sample,
                          do.center = FALSE,
                          do.scale  = TRUE)

snare_test <- NormalizeData(snare_test)
## This effectively just extracts normalized counts for variable genes
snare_test <- ScaleData(snare_test,
                        do.center = FALSE,
                        do.scale  = TRUE)

snare_expr_mat <- snare_sample@assays$RNA@scale.data

snare_test_expr_mat <- snare_test@assays$RNA@scale.data


# 
# ### Compute Full SIDESeq ------------------------------------------------------
# full_side_seq <-
#   SIDEseqSimult(expr_matrix = snare_expr_mat, 
#                 diff_expr_method = "diff_expr_norm",
#                 similarity_method = "n_intersect",
#                 n_top_genes = 150,
#                 parallelize = TRUE,
#                 n_cores = 23)
# 
# save(full_side_seq, file = here("output/snare_sideseq/side_seq_snare_1500.RData"))
# 
# 
# ### Run multiple iterations at different ref set sizes ------------------------
# 
# set.seed(448213)
# 
# ref_set_sizes <- c(25, 50, 75, 100, 150, 200, 300, 400) 
# 
# 
# stable_num <- 5
# 
# for(r in ref_set_sizes) {
#   print(r)
#   for(s in 1:stable_num) {
#     ptm <- proc.time()
#     side_seq_ref_rand <- 
#       SIDEseqRefSet(expr_matrix = snare_expr_mat,
#                     diff_expr_method = "diff_expr_norm",
#                     similarity_method = "n_intersect",
#                     n_top_genes = 150,
#                     selection_method = "random",
#                     size_ref_set = r,
#                     n_clust = 25, 
#                     B=0, 
#                     D=1,
#                     parallelize = TRUE,
#                     n_cores = 23,
#                     verbose = TRUE)
#     
#     ct_elap <- (proc.time() - ptm)[3]
#     
#     side_seq_ref_rand[[3]] <- ct_elap 
#     
#     save(side_seq_ref_rand, 
#          file = here(paste0("output/snare_sideseq/run", s, 
#                             "_", r, "rand_ref_set_snare_1500.RData")))
#     
#     
#     
#     ptm <- proc.time()
#     
#     side_seq_ref_cell_embed <- 
#       SIDEseqRefSet(expr_matrix = snare_expr_mat,
#                     diff_expr_method = "diff_expr_norm",
#                     similarity_method = "n_intersect",
#                     n_top_genes = 150,
#                     selection_method = "cell_embed_sample",
#                     size_ref_set = r,
#                     n_clust = 25, 
#                     B = 3, 
#                     D = 1,
#                     parallelize = TRUE,
#                     n_cores = 23,
#                     verbose = TRUE)
#     
#     ct_elap <- (proc.time() - ptm)[3]
#     
#     side_seq_ref_cell_embed[[3]] <- ct_elap 
#     
#     save(side_seq_ref_cell_embed, 
#          file = here(paste0("output/snare_sideseq/run", s, 
#                             "_", r, "cell_embed_ref_set_snare_1500.RData")))
#     
#     
#   }
# }
# 


### Testing stratified initialization, one run ---------------------------------

set.seed(981734)

ref_set_sizes <- c(25, 50, 75, 100, 150, 200, 300, 400) 


stable_num <- 5

for(r in ref_set_sizes) {
  print(r)
  for(s in 1:stable_num) {

    ptm <- proc.time()
    
    side_seq_ref_cell_embed_pca_umap <- 
      SIDEseqRefSet(expr_matrix = snare_expr_mat,
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = 150,
                    selection_method = "cell_embed_sample",
                    size_ref_set = r,
                    n_clust = 25, 
                    B = 0, 
                    D = 1,
                    parallelize = TRUE,
                    n_cores = 23,
                    verbose = TRUE)
    
    ct_elap <- (proc.time() - ptm)[3]
    
    side_seq_ref_cell_embed_pca_umap[[3]] <- ct_elap 
    
    save(side_seq_ref_cell_embed_pca_umap, 
         file = here(paste0("output/snare_sideseq/run", s, 
                            "_", r, "cell_embed_pca_umap_ref_set_snare_1500.RData")))
    
    
  }
}



# 
# ### Compute analyses for different levels of top genes ------------------------
# set.seed(1739287)
# 
# gene_nums <- c(50, 150, 300, 450, 600, 900)
# 
# for(g in gene_nums) {
#   side_seq_ref_rand <- 
#     SIDEseqRefSet(expr_matrix = snare_expr_mat,
#                   diff_expr_method = "diff_expr_norm",
#                   similarity_method = "n_intersect",
#                   n_top_genes = g,
#                   selection_method = "random",
#                   size_ref_set = 100,
#                   n_clust = 25, 
#                   B=0, 
#                   D=1,
#                   parallelize = TRUE,
#                   n_cores = 23,
#                   verbose = TRUE)
#   
#   save(side_seq_ref_rand, 
#        file = here(paste0("output/snare_sideseq/side_seq_rand_", g, 
#                           "_genes.RData")))
# }
# 
# 
# 
# 
# ### Run for the full test set data --------------------------------------------
# set.seed(7489992)
# 
# ptm <- proc.time()
# 
# side_seq_ref_rand <- 
#   SIDEseqRefSet(expr_matrix = snare_test_expr_mat,
#                 diff_expr_method = "diff_expr_norm",
#                 similarity_method = "n_intersect",
#                 n_top_genes = 150,
#                 selection_method = "random",
#                 size_ref_set = 100,
#                 n_clust = 25, 
#                 B=0, 
#                 D=1,
#                 parallelize = TRUE,
#                 n_cores = 23,
#                 verbose = TRUE)
# 
# ct_elap <- (proc.time() - ptm)[3]
# 
# side_seq_ref_rand[[3]] <- ct_elap 
# 
# 
# save(side_seq_ref_rand, 
#      file = here("output/snare_sideseq/side_seq_test_100_cells_300_genes.RData"))
# 
