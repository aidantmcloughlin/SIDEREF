###############################################################################
### Low Level Cluster Analysis
###############################################################################

set.seed(1)

library(here)
source(here('R/libraries.R'))
source(here('R/tab_muris_clust_pipeline.R'))
source(here('R/hcclust.R'))
source(here('R/subgroup_composition_gene_contrib.R'))

############################# CONSTANTS ####################################
## Data subset Constants
KEEP_DROPLET <- TRUE
KEEP_FACS <- TRUE
data_pct <- "0.05pct"
file_name <- "tab_muris_full_" %p% data_pct

## Distance setting constants:
SPECTRAL_DIMS <- 10
G <- c(50, 150, 300)
R <- 100
PCs <- c(3, 10, 25)
N_CLUST <- 20 ## k-means clusters during SIDEREF
VAR_FEATURES <- 3000

## Other:
N_CORES <- 1

## Whether it's needed to run the preproc pipeline:
RUN_DIST_AND_CLUST <- FALSE




### load Tabula Muris data
load(here('output/tab_muris_sc/preproc_data/' %p% file_name %p% '_scale_var_genes_3000.RData'),
     verbose = TRUE)
load(here('output/tab_muris_sc/preproc_data/' %p% file_name %p% '.RData'),
     verbose = TRUE)

### full data distance results:


droplet_cells <- which(grepl("droplet", meta_data_full$file_source))
facs_cells <- which(grepl("facs", meta_data_full$file_source))

if(!KEEP_FACS & KEEP_DROPLET) {
  data_full <- data_full[,droplet_cells]
  meta_data_full <- meta_data_full[droplet_cells,]
  ## TODO: compute
  #load(DROPLET_FULL_DISTANCE_FILE)
  
  file_tag = "droplet_only"
  
} else if(KEEP_FACS & !KEEP_DROPLET) {
  data_full <- data_full[,facs_cells]
  meta_data_full <- meta_data_full[facs_cells,]
  ## TODO: compute
 # load(FACS_FULL_DISTANCE_FILE)
  
  file_tag = "facs_only"
} else{
  full_dist_res <- 
    loadRData(here("output/tab_muris_sc/dist_res/full_tab_muris_dist_res_" %p%
                     data_pct %p% "_" %p% file_tag %p% ".RData"))
  
  file_tag = "droplet_and_facs"
}



tm_filenames <- unique(meta_data_full$file_source)

tm_clust_by_file <- list()


if(RUN_DIST_AND_CLUST) {
  ## File-specific clustering results:
  j=1
  for(i in seq_len(length(tm_filenames))) {
    print(tm_filenames[i])
    tm_clust_by_file[[j]] <- 
      TabMurisClustResPipeline(data_full, meta_data_full, 
                               filename = tm_filenames[i],
                               spectral_dims = SPECTRAL_DIMS, 
                               n_cores = N_CORES, 
                               R = R, 
                               G = G,
                               PCs = PCs,
                               var_features = VAR_FEATURES)
    
    names(tm_clust_by_file)[j] <- tm_filenames[i]
    j = j+1    			     
  }
  
  save(tm_clust_by_file, 
       file = here("output/tab_muris_sc/clust_res/tm_clust_by_file_" %p% 
                     data_pct %p% "_" %p% file_tag %p% ".RData"))
  
  
#  ## Overall-distance matrices:
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

 save(full_dist_res,
      file = here("output/tab_muris_sc/dist_res/full_dist_res_" %p%
                    data_pct %p% "_" %p% file_tag %p% ".RData"))
  
}






###############################################################################
### Results assessment
###############################################################################


load(here("output/tab_muris_sc/clust_res/tm_clust_by_file_" %p% 
            data_pct %p% "_" %p% file_tag %p% ".RData"), 
     verbose=TRUE)


## append cell barcodes to cluster results
## TODO: this should be deprecated in a rerun

# for(t in tm_filenames) {
#   print(t)
#   sub_res <- TabMurisProcSubset(data_full, meta_data_full, t,
#                                 var_features = VAR_FEATURES) 
#   
#   full_data_subset_variable <- sub_res[[1]]
#   cell_inds <- sub_res[[2]]
#   cell_names <- colnames(full_data_subset_variable)
#   tm_clust_by_file[[t]] = 
#     lapply(tm_clust_by_file[[t]], function(x) {
#       names(x) = cell_names
#       return(x)})
# }


## gather cluster res from all distance measures into organized dataframe.
cell_name_vec <- 
  unlist(lapply(tm_clust_by_file, function(x) names(x$clusts[[1]])))

cell_filename_vec <- 
  unlist(lapply(names(tm_clust_by_file), function(x) {
    rep(x, length(tm_clust_by_file[[x]]$clusts[[1]]))
  }))

clust_res_df <- 
  data.frame(file_source = cell_filename_vec, 
             cell_id = cell_name_vec) %>% 
  left_join(meta_data_full %>% select(cell_id, cell_type))

## cluster res by distance measure:
for(d in names(tm_clust_by_file[[1]]$clusts)) {
  
  clust_res_vec <-
    unlist(lapply(tm_clust_by_file, function(x) x$clusts[[d]]))
  
  
  clust_res_df['clust_res_' %p% d] <- 
    as.numeric(factor(clust_res_df$file_source %p% as.character(clust_res_vec))) 
  
}




clust_res_cols <- names(clust_res_df)[4:ncol(clust_res_df)]


## compute ARI wrt OG clusters
adj_rand_index_og <- 
  clust_res_df %>% 
  summarise_at(clust_res_cols, ~ mclust::adjustedRandIndex(cell_type, .x))

adj_rand_index_sideref50 <- 
  clust_res_df %>% 
  summarise_at(setdiff(clust_res_cols, "clust_res_side_ref_g50_dist"), 
               ~ mclust::adjustedRandIndex(clust_res_side_ref_g50_dist, .x))

adj_rand_index_pca10 <- 
  clust_res_df %>% 
  summarise_at(setdiff(clust_res_cols, "clust_res_pca_10_dist"), 
               ~ mclust::adjustedRandIndex(clust_res_pca_10_dist, .x))


adj_rand_index_sideref50 %>% 
  gather(dist_measure, ari) #%>% bb
# ggplot(aes(x=dist_measure,y=ari)) + 
# geom_bar(stat="identity")

adj_rand_index_pca10 %>% 
  gather(dist_measure, ari) 


clust_res_df %>% 
  summarise_at(clust_res_cols, max) %>% 
  gather(dist_measure, num_clust)

### modularity scores
n_dists = length(tm_clust_by_file[[1]]$clusts)
dist_names = names(tm_clust_by_file[[1]]$clusts)

repped_file_names = sapply(names(tm_clust_by_file), function(x) rep(x, n_dists))

nn_module = lapply(tm_clust_by_file, function(x) {
  out = sapply(x$modularities, function(y) {
    return(y[1])
  })
  names(out) = dist_names
  return(out)
  
})

nn_nums = unlist(nn_module)
names(nn_nums) = NULL

nn_modularity_df <-
  data.frame(file_name = c(repped_file_names),
             nn = nn_nums,
             dist = rep(names(nn_module[[1]]), length(tm_clust_by_file)))

nn_summary<-
  nn_modularity_df %>% 
  mutate(dist = factor(dist, levels = unique(dist))) %>% 
  group_by(dist) %>% 
  summarise(mean_nn = mean(nn), 
            median_nn = median(nn)) %>% head(24)

## Average Silhouette Scores of UMAP embeddings

avg_sils = lapply(tm_clust_by_file, function(x) {
  out = sapply(seq_len(n_dists), function(i) {
    c = x$clusts[[i]]
    umap_res = x$umaps[[i]]
    if(length(unique(c)) == 1) {
      avg_sil = NA
    } else{
      sil_scores <- cluster::silhouette(as.integer(factor(c)),
                                        dist = stats::dist(umap_res, 
                                                           method = "euclidean"))[,3]
      avg_sil = mean(sil_scores)  
    }
    
    return(avg_sil)
  })
  names(out) = dist_names
  return(out)
  
})

avg_sils_nums = unlist(avg_sils)
names(avg_sils_nums) = NULL

avg_sils_df <-
  data.frame(file_name = c(repped_file_names),
             avg_sil = avg_sils_nums,
             dist = rep(names(nn_module[[1]]), length(tm_clust_by_file)))

avg_sils_sum <- 
  avg_sils_df %>% 
  mutate(dist = factor(dist, levels = unique(dist))) %>% 
  group_by(dist) %>% 
  summarise(mean_avg_sil = mean(avg_sil, na.rm=TRUE))





## Use clusters to generate tissue-level dendrograms:

## HClust Results (wrt automatic low-level clusters)

#load(here("output/tab_muris_sc/full_tab_muris_dist_res_0.05pct_3000g.RData"), verbose=TRUE)

if(SORT_INCONSISTENCIES <- TRUE) {
  extra_dists <- 
    setdiff(names(full_dist_res$dist_list), names(tm_clust_by_file[[1]]$clusts))
  
  extra_dists_inds <- which(names(full_dist_res$dist_list) %in% extra_dists)
  
  if(!is_empty(extra_dists_inds)) {
    full_dist_res$dist_list <- full_dist_res$dist_list[-extra_dists_inds]  
  }
  
  
  missing_dists <- 
    setdiff(names(tm_clust_by_file[[1]]$clusts), names(full_dist_res$dist_list))
  
  missing_dists_inds <- which(names(tm_clust_by_file[[1]]$clusts) %in% missing_dists)
  
  tm_clust_by_file_sub <- 
    lapply(tm_clust_by_file, function(x) {
      x$clusts       <- x$clusts[-missing_dists_inds]
      x$modularities <- x$modularities[-missing_dists_inds]
      x$umaps        <- x$umaps[-missing_dists_inds]
      
      return(x)
    }) 
  
  if(ncol(full_dist_res$dist_list[[1]]) > nrow(meta_data_full)) {
    if(file_tag=="droplet_only"){inds= droplet_cells}else{inds=facs_cells}
    for(n in names(full_dist_res$dist_list)) {
      full_dist_res$dist_list[[n]] <- full_dist_res$dist_list[[n]][inds, inds] 
    }
    
  }
}


### 1: With Respect to their Baseline Distance Measurement:
tab_muris_hclust_res_to_same_dist_clust <-
  lapply(names(full_dist_res$dist_list),
         function(n) {
           return(
             RunHClustPipeline(
               clusters = as.numeric(factor(clust_res_df["clust_res_" %p% n])), 
               subg_inds = as.numeric(factor(meta_data_full$tissue)),
               hclust_dissim_mat = full_dist_res$dist_list[[n]], 
               p_dissim_mat      = 
                 full_dist_res$dist_list$pca_25_dist))
         })


names(tab_muris_hclust_res_to_same_dist_clust) <- 
  names(full_dist_res$dist_list)

lapply(tab_muris_hclust_res_to_same_dist_clust, 
       function(x) x$adj_rand_index)


### 2: With Respect to PCA 3 (wtd):
tab_muris_hclust_res_to_pca_25_dist_clust <-
  lapply(names(full_dist_res$dist_list),
         function(n) {
           return(
             RunHClustPipeline(
               clusters = as.numeric(factor(clust_res_df["clust_res_pca_wtd25_dist"])), 
               subg_inds = as.numeric(factor(meta_data_full$tissue)),
               hclust_dissim_mat = full_dist_res$dist_list[[n]], 
               p_dissim_mat      = 
                 ## TODO: adjust this
                 full_dist_res$dist_list$pca_25_dist))
         })


names(tab_muris_hclust_res_to_pca_25_dist_clust) <- 
  names(full_dist_res$dist_list)

lapply(tab_muris_hclust_res_to_pca_25_dist_clust, 
       function(x) x$adj_rand_index)




### 4: To Cell Type
tab_muris_hclust_res_to_celltype_clust <-
  lapply(names(full_dist_res$dist_list),
         function(n) {
           return(
             RunHClustPipeline(
               clusters = as.numeric(factor(clust_res_df$file_source %p% clust_res_df$cell_type)), 
               subg_inds = as.numeric(factor(meta_data_full$tissue)),
               hclust_dissim_mat = full_dist_res$dist_list[[n]], 
               p_dissim_mat      = 
                 ## TODO: adjust this
                 full_dist_res$dist_list$pca_25_dist))
         })


names(tab_muris_hclust_res_to_celltype_clust) <- 
  names(full_dist_res$dist_list)

lapply(tab_muris_hclust_res_to_celltype_clust, 
       function(x) x$adj_rand_index)



### 3: With Respect to "Best" Low Level Clusters? (data adaptive way for this?)
# return(
#   RunHClustPipeline(
#     clusters = as.numeric(factor(clust_res_df["clust_res_pca_10_dist"])), 
#     subg_inds = as.numeric(factor(meta_data_full$tissue)),
#     hclust_dissim_mat = full_dist_res$dist_list[[n]], 
#     p_dissim_mat      = 
#       ## TODO: adjust this
#       full_dist_res$dist_list$pca_25_dist))
# })
# 
# 
# names(tab_muris_hclust_res_to_pca_10_dist_clust) <- 
#   names(full_dist_res$dist_list)
# 
# lapply(tab_muris_hclust_res_to_pca_10_dist_clust, 
#        function(x) x$adj_rand_index)
# 




## from cell-level to tissue-level dendrograms.
## Repeated analysis with no low level cluster assumptions

tab_muris_hclust_res_from_cell_level_res_df <-
  lapply(names(full_dist_res$dist_list),
         function(n) {
           return(
             RunHClustPipeline(
               clusters = seq_len(nrow(meta_data_full)), 
               subg_inds = as.numeric(factor(meta_data_full$tissue)),
               hclust_dissim_mat = full_dist_res$dist_list[[n]], 
               p_dissim_mat      = 
                 ## TODO: adjust this
                 full_dist_res$dist_list$pca_25_dist))
         })


names(tab_muris_hclust_res_from_cell_level_res_df) <- 
  names(full_dist_res$dist_list)

lapply(tab_muris_hclust_res_from_cell_level_res_df, 
       function(x) x$adj_rand_index)







## TODO: clustering algorithm applied to full data each distance measure 
## (separate script)


### ^^This is a compute job (Simply Run My Clustering Function on the full distance matrix.)


## Marker Gene Analysis




## install.packages("diceR")
##
# library(diceR)
# 
# 
# data(hgsc)
# dat <- hgsc[1:100, 1:50]
# x <- consensus_cluster(dat, nk = 4, reps = 1, algorithms = c("hc", "diana"),
#                        progress = FALSE)
# 
# CSPA(x, k = 4)




