###############################################################################
### Low Level Cluster Analysis
###############################################################################

set.seed(1)

library(here)
print(here())
source(here('R/libraries.R'))
source(here('R/dist_and_clust_pipelines.R'))
source(here('R/subgroup_composition_gene_contrib.R'))
source(here('R/hcclust.R'))


############################# CONSTANTS ####################################
## Data subset Constants
KEEP_DROPLET <- FALSE
KEEP_FACS <- TRUE

data_pct <- ""
file_name <- "tab_muris_full" %p% ifelse(data_pct == "", "", "_" %p% data_pct)
print(file_name)
REPS <- 10
FOLDS <- 20
FOLD_SEEDS <- 100:(100+REPS-1)
## sampling cell types by size 20 per type:
SAMPLE_CELL_TYPES = TRUE
SAMPLE_PCT = FALSE
SAMPLE_SIZE = 20

## Distance setting constants:
SPECTRAL_DIMS <- 25
G <- c(50, 150, 300)
R <- 100
PCs <- c(3, 10, 25)
N_CLUST <- 20 ## k-means clusters during SIDEREF
VAR_FEATURES <- 3000


## Other:
N_CORES <- 16
 

## Whether it's needed to run the preproc pipeline:
RUN_DIST <- TRUE
RUN_RES <- FALSE


tm_data_list <- 
  tmLoadCut(file_name = file_name, data_pct = data_pct, 
            sample_size = SAMPLE_SIZE,
            sample_cell_types = SAMPLE_CELL_TYPES,
            keep_droplet = KEEP_DROPLET, keep_facs = KEEP_FACS,
            sample_pct = SAMPLE_PCT)

## extract objects
for(n in names(tm_data_list)) {
  assign(n, tm_data_list[[n]])
}

rm(tm_data_list)


for(x in data_store_list[1:length(data_store_list)]) {

    cat("running fold " %p% x[[1]])
    set.seed(data_store_list[[1]][[4]])


    tm_clust_by_file <- list()

    data_full <- x[[2]]
    print("outer loop object track: " %p% dim(data_full))
    meta_data_full <- x[[3]]

    if(RUN_DIST) {
      print("running full distance computation")
      ## Overall-distance matrices:
      full_dist_res <- 
        TabMurisDistResPipeline(data_full, 
                                meta_data_full, 
                                filename = NULL, 
                                spectral_dims = SPECTRAL_DIMS, 
                                n_clust = N_CLUST,
                                n_cores = N_CORES,
                                R = R,
                                G = G,
                                PCs = PCs,
                                var_features = VAR_FEATURES)
      
      print(mem_used())
      ### save data 
      if(x[[1]] == 1){
        save(full_dist_res, meta_data_full,
             file = here("output/tab_muris_sc/dist_res/full_dist_res_" %p%
               data_pct %p% "_" %p% file_tag %p% "_fold" %p% x[[1]] %p% 
                 "_samp_ct" %p% SAMPLE_CELL_TYPES %p% 
                 "_samp_pct" %p% SAMPLE_PCT %p%
                 "_samp_size" %p% SAMPLE_SIZE %p% ".RData")) }

    }

   
    
    ###############################################################################
    ### Results assessment
    ###############################################################################
    
    
    if(RUN_RES) {
      
      print("beginning results assessment..")
      
      results_list <- list()
      
      
      ## cluster results
      #	load(here("output/tab_muris_sc/clust_res/tm_clust_by_file_" %p% 
      #			   data_pct %p% "_" %p% file_tag %p% "_fold" %p% x[[1]]  %p% ".RData"), 
      #	     verbose=TRUE)
      #
      ## full distance matrix (and full clustering results, if produced)
      #	load(here("output/tab_muris_sc/dist_res/full_dist_res_" %p% data_pct %p% "_" %p% file_tag %p% "_fold" %p% x[[1]] %p% ".RData"))
      
      
      
      #	print("creating cluster results dataframe")
      
      
      #	## gather cluster res from all distance measures into organized dataframe.
      #	cell_name_vec <- 
      #	  unlist(lapply(tm_clust_by_file, function(x) names(x$clusts[[1]])))
      #
      #	cell_filename_vec <- 
      #	  unlist(lapply(names(tm_clust_by_file), function(x) {
      #	    rep(x, length(tm_clust_by_file[[x]]$clusts[[1]]))
      #	  }))
      #
      #	clust_res_df <- 
      #	  data.frame(file_source = cell_filename_vec, 
      #		     cell_id = cell_name_vec) %>% 
      #	  left_join(meta_data_full %>% select(cell_id, cell_type))
      #
      #	## cluster res by distance measure:
      #	for(d in names(tm_clust_by_file[[1]]$clusts)) {
      #
      #	  clust_res_vec <-
      #	    unlist(lapply(tm_clust_by_file, function(x) x$clusts[[d]]))
      #
      #
      #	  clust_res_df['clust_res_' %p% d] <- 
      #	    as.numeric(factor(clust_res_df$file_source %p% as.character(clust_res_vec))) 
      #
      #	}
      #
      #	results_list[["clust_res_df"]] <- clust_res_df
      #
      #
      #	clust_res_cols <- names(clust_res_df)[4:ncol(clust_res_df)]
      #
      #
      #	## compute ARI wrt original clusters
      #	adj_rand_index_celltype <- 
      #	  clust_res_df %>% 
      #	  summarise_at(clust_res_cols, ~ mclust::adjustedRandIndex(cell_type, .x))
      #
      #	adj_rand_index_sideref50 <- 
      #	  clust_res_df %>% 
      #	  summarise_at(setdiff(clust_res_cols, "clust_res_side_ref_g50_dist"), 
      #		       ~ mclust::adjustedRandIndex(clust_res_side_ref_g50_dist, .x))
      #
      #	adj_rand_index_pca10 <- 
      #	  clust_res_df %>% 
      #	  summarise_at(setdiff(clust_res_cols, "clust_res_pca_10_dist"), 
      #		       ~ mclust::adjustedRandIndex(clust_res_pca_10_dist, .x))
      #
      #	results_list[["ari_wrt_celltype"]] <-
      #	  adj_rand_index_celltype %>%
      #	    gather(dist_measure_ari)
      #
      #
      #	results_list[["ari_wrt_sideref50"]] <-
      #	  adj_rand_index_sideref50 %>%
      #	    gather(dist_measure_ari)
      #
      #
      #	results_list[["ari_wrt_pca10"]] <-
      #	  adj_rand_index_pca10 %>%
      #	    gather(dist_measure_ari)
      #
      #
      #	results_list[["num_clusts"]] <-
      #	  clust_res_df %>% 
      #	  summarise_at(clust_res_cols, max) %>% 
      #	  gather(dist_measure, num_clust)
      #
      #
      #	print("computing average modularity scores across files")
      #
      #
      #	### modularity scores
      #	n_dists = length(tm_clust_by_file[[1]]$clusts)
      #	dist_names = names(tm_clust_by_file[[1]]$clusts)
      #
      #	repped_file_names = sapply(names(tm_clust_by_file), function(x) rep(x, n_dists))
      #
      #
      #	for(i in 1:2) {
      #	  ## 1 corresponding to NN
      #	  ## 2 corresponding to SNN
      #
      #	  nn_module = lapply(tm_clust_by_file, function(x) {
      #	    out = sapply(x$modularities, function(y) {
      #	      return(y[i])
      #	    })
      #	    names(out) = dist_names
      #	    return(out)
      #
      #	  })
      #
      #	  nn_nums = unlist(nn_module)
      #	  names(nn_nums) = NULL
      #
      #	  nn_modularity_df <-
      #	    data.frame(file_name = c(repped_file_names),
      #		       nn = nn_nums,
      #		       dist = rep(names(nn_module[[1]]), length(tm_clust_by_file)))
      #
      #	  results_list[[ifelse(i==1,"nn_res","snn_res")]] <-
      #	    nn_modularity_df %>%
      #	    mutate(dist = factor(dist, levels = unique(dist))) %>% 
      #	    group_by(dist) %>% 
      #	    summarise(mean_nn = mean(nn), 
      #		      median_nn = median(nn))
      #
      #	}
      #
      #	print("computing avg sil scores")
      #
      #	## Average Silhouette Scores of UMAP embeddings
      #
      #	avg_sils = lapply(tm_clust_by_file, function(x) {
      #	  out = sapply(seq_len(n_dists), function(i) {
      #	    c = x$clusts[[i]]
      #	    umap_res = x$umaps[[i]]
      #	    if(length(unique(c)) == 1) {
      #	      avg_sil = NA
      #	    } else{
      #	      sil_scores <- cluster::silhouette(as.integer(factor(c)),
      #						dist = stats::dist(umap_res, 
      #								   method = "euclidean"))[,3]
      #	      avg_sil = mean(sil_scores)  
      #	    }
      #
      #	    return(avg_sil)
      #	  })
      #	  names(out) = dist_names
      #	  return(out)
      #
      #	})
      #
      #	avg_sils_nums = unlist(avg_sils)
      #	names(avg_sils_nums) = NULL
      #
      #	avg_sils_df <-
      #	  data.frame(file_name = c(repped_file_names),
      #		     avg_sil = avg_sils_nums,
      #		     dist = rep(names(nn_module[[1]]), length(tm_clust_by_file)))
      #
      #	results_list[["avg_sil_scores_UMAP"]] <-
      #	  avg_sils_df %>% 
      #	  mutate(dist = factor(dist, levels = unique(dist))) %>% 
      #	  group_by(dist) %>% 
      #	  summarise(mean_avg_sil = mean(avg_sil, na.rm=TRUE))
      #
      #
      #	avg_sil_full_dist <-
      #	  data.frame(dist = names(tm_clust_by_file[[1]]$avg_sils),
      #		     avg_sil = sapply(names(tm_clust_by_file[[1]]$avg_sils), 
      #				      function(x) { mean(
      #					sapply(tm_clust_by_file, function(t) {
      #					  return(t$avg_sils[[x]])
      #					}))
      #				      }))
      #
      #	rownames(avg_sil_full_dist) <- NULL
      #
      #	results_list[["avg_sil_scores_fulldist"]] <-
      #	  avg_sil_full_dist
      #
      #
      #
      #
      #
      #	## Use clusters to generate tissue-level dendrograms:
      
      
      ##############################################################################
      ###   Gather the needed idx lists for hclust
      
      
      annotated_cells <- which(!is.na(meta_data_full$cell_type))
      
      tissue_subg_inds_all <- 
        as.numeric(
          factor(meta_data_full$tissue,
                 levels = unique(meta_data_full$tissue)))
      
      tissue_subg_inds_annotated <- 
        as.numeric(
          factor(meta_data_full[annotated_cells,]$tissue,
                 levels = unique(meta_data_full[annotated_cells,]$tissue)))
      
      celltype_g_inds <- 
        as.numeric(
          factor(meta_data_full[annotated_cells,]$cell_type,
                 levels = unique(meta_data_full[annotated_cells,]$cell_type))
        )
      
      # celltype_pca10_inds <- 
      #   as.numeric(factor(clust_res_df$clust_res_pca_10_dist),
      #              levels = unique(clust_res_df$clust_res_pca_10_dist))
      
      #
      #	print("hierarchical clustering results")
      #
      #	## HClust Results (wrt automatic low-level clusters)
      #
      #	print("wrt baseline")
      #	### 1: With Respect to their Baseline Distance Measurement:
      #	tab_muris_hclust_res_to_same_dist_clust <-
      #	  lapply(names(full_dist_res$dist_list),
      #		 function(n) {
      #		   return(
      #		     RunHClustPipeline(
      #		       clusters = as.numeric(factor(clust_res_df["clust_res_" %p% n])), 
      #		       subg_inds = tissue_subg_inds_all,
      #		       hclust_dissim_mat = full_dist_res$dist_list[[n]], 
      #		       p_dissim_mat      = 
      #			 full_dist_res$dist_list$pca_25_dist))
      #		 })
      #
      #
      #	names(tab_muris_hclust_res_to_same_dist_clust) <- 
      #	  names(full_dist_res$dist_list)
      #
      #	results_list[["tab_muris_hclust_res_to_same_dist_clust"]] <-
      #	  tab_muris_hclust_res_to_same_dist_clust
      #
      #	results_list[["hclust_ari_baseline"]] <-
      #	  lapply(tab_muris_hclust_res_to_same_dist_clust, 
      #		 function(x) x$adj_rand_index_high)
      #
      #
      #	print("wrt pca10")
      #	### 2: With Respect to PCA 10:
      #	tab_muris_hclust_res_to_pca_10_dist_clust <-
      #	  lapply(names(full_dist_res$dist_list),
      #		 function(n) {
      #		   return(
      #		     RunHClustPipeline(
      #		       clusters = celltype_pca10_inds, 
      #		       subg_inds = tissue_subg_inds_all,
      #		       hclust_dissim_mat = full_dist_res$dist_list[[n]], 
      #		       p_dissim_mat      = 
      #			 full_dist_res$dist_list$pca_25_dist))
      #		 })
      #
      #
      #	names(tab_muris_hclust_res_to_pca_10_dist_clust) <- 
      #	  names(full_dist_res$dist_list)
      #
      #	results_list[["hclust_res_pca10"]] <-
      #	  tab_muris_hclust_res_to_pca_10_dist_clust
      #
      #	results_list[["hclust_ari_pca10"]] <-
      #	  lapply(tab_muris_hclust_res_to_pca_10_dist_clust, 
      #		 function(x) x$adj_rand_index_high)
      #
      #
      #
      #
      print("wrt cell level")
      ## 3. from cell-level to tissue-level dendrograms.
      tab_muris_hclust_res_from_cell_level_res_df <-
        lapply(full_dist_res$dist_list,
               function(x) {
                 res <- 
                   RunHClustPipeline(
                     clusters = seq_len(nrow(meta_data_full)), 
                     subg_inds = tissue_subg_inds_all,
                     hclust_dissim_mat = x, 
                     p_dissim_mat      = full_dist_res$dist_list$pca_25_dist,
                     cell_ids = meta_data_full$cell_id)
               })
      
      
      names(tab_muris_hclust_res_from_cell_level_res_df) <- 
        names(full_dist_res$dist_list)
      
      results_list[["hclust_res_celllevel"]] <-
        tab_muris_hclust_res_from_cell_level_res_df
      
      results_list[["hclust_ari_celllevel"]] <-
        lapply(tab_muris_hclust_res_from_cell_level_res_df, 
               function(x) x$adj_rand_index_high)
      
      
      
      print("wrt cell types")
      ## 4. With pre-computed "cell types"
      
      tab_muris_hclust_res_with_celltype <-
        lapply(full_dist_res$dist_list,
               function(x) {
                 res <- 
                   RunHClustPipeline(
                     clusters = celltype_g_inds, 
                     subg_inds = tissue_subg_inds_annotated,
                     hclust_dissim_mat = x[annotated_cells, annotated_cells], 
                     p_dissim_mat      = full_dist_res$dist_list$pca_25_dist[annotated_cells, 
                                                                             annotated_cells],
                     cell_ids = meta_data_full[annotated_cells,]$cell_id)
               })
      
      
      names(tab_muris_hclust_res_with_celltype) <- 
        names(full_dist_res$dist_list)
      
      results_list[["hclust_res_celltypes"]] <-
        tab_muris_hclust_res_with_celltype
      
      results_list[["hclust_ari_celltypes"]] <-
        lapply(tab_muris_hclust_res_with_celltype, 
               function(x) x$adj_rand_index_high)
      
      
      
      ## 5. Starting from the tissue level.
      tab_muris_hclust_res_at_tissue <-
        lapply(full_dist_res$dist_list,
               function(x) {
                 res <- 
                   HClust_with_start_clusters(
                     clusters = tissue_subg_inds_annotated, 
                     subg_inds = rep(1, length(annotated_cells)),
                     dissim_mat = x[annotated_cells, annotated_cells], 
                     cell_ids = meta_data_full[annotated_cells,]$cell_id)
               })
      
      names(tab_muris_hclust_res_at_tissue) <- 
        names(full_dist_res$dist_list)
      
      results_list[["hclust_res_tissue"]] <-
        tab_muris_hclust_res_at_tissue
      
      results_list[["meta_data"]] <- 
        meta_data_full
      
      save(results_list, file = here("output/tab_muris_sc/no_low_level_hclust_results_list_" %p%
                                       data_pct %p% "_" %p% file_tag %p% "_fold" %p% x[[1]] %p% "of" %p% FOLDS %p%  ".RData"))
      
      
    }
    
    
}
