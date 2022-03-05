
###############################################################################
### Differential Expression Test with conjoined clusters
###############################################################################





###############################################################################
### Comparing Similarity Scores from SIDEREF
###############################################################################

print("running full distance computation")

set.seed(data_store_list[[1]][[4]])

data_full <- data_store_list[[1]][[2]]
meta_data_full <- data_store_list[[1]][[3]]

sub_res <- TabMurisProcSubset(full_data = data_full, 
                              meta_data = meta_data_full,
                              var_features = VAR_FEATURES,
                              return_seurat = TRUE)

full_data_seurat <- sub_res[[1]]

full_data_subset_variable <- full_data_seurat@assays$RNA@scale.data
cell_inds <- sub_res[[2]]

## reset n_clust for small data
n_clust = min(floor(ncol(full_data_subset_variable) / 3), N_CLUST)

cells_to_keep <- which(meta_data_full$tissue == "Marrow")

## Overall-distance matrices:
sideref_res <- 
  SIDEseqRefSet(expr_matrix = full_data_subset_variable[, cells_to_keep], 
                n_top_genes = 150,
                size_ref_set = R,
                selection_method = "cell_embed_sample",
                return_sim = TRUE,
                n_clust = N_CLUST,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = N_CORES,
                ## other
                verbose = TRUE)

sideseq_res <- SIDEseqSimult(full_data_subset_variable[, cells_to_keep], 
                             n_top_genes = 150,
                             parallelize = TRUE,
                             n_cores = N_CORES, 
                             return_sim = TRUE)



compareSimScore <- function(sim_mat, cell_idx) {
  sim_as_vec <- sim_mat[lower.tri(sim_mat)]
  
  
  cell_group_sim <- sim_mat[cell_idx, cell_idx]
  cell_group_vs_sim <- 
    sim_mat[cell_idx, setdiff(seq_len(ncol(sim_mat)), cell_idx)]
  
  ## remove zeroes along diagonal 
  cell_group_sim_as_vec <- cell_group_sim[lower.tri(cell_group_sim)]
  
  within_sim_avg = mean(cell_group_sim_as_vec)
  vs_sim_avg = mean(cell_group_vs_sim)
  
  return(list(within_sim_avg = within_sim_avg,
              vs_sim_avg = vs_sim_avg))
}



