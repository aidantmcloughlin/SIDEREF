#### Snare-seq SIDEREF runs

## Get Distance Matrix of Full SIDEREF SNAREseq Distance Matrix.
load(here("output/snare_sideseq/side_ref_test_150_cells_150_genes.RData"))

## Get actual test data for cell type labels.
load(here("output/snare_test.RData"))

snare_test_cell_types <- snare_test@meta.data$celltype
snare_test_expr_mat <- snare_test@assays$RNA@scale.data

snare_test_cols <- ncol(snare_test_expr_mat)


### (A) SIDEREF RES ------------------------------------------------------------

### Plotting
snare_group_p_hclust <-
  groupwiseDistanceHeatmap(
    group_labels = snare_test_cell_types, 
    dist_mat = side_seq_ref_rand$dissim_final, 
    title = "SNARE Data Cell Types\nSIDEREF (150C, 150G)\nClustered",
    hclust = TRUE)

snare_group_p_no_hclust <-
  groupwiseDistanceHeatmap(
    group_labels = snare_test_cell_types, 
    dist_mat = side_seq_ref_rand$dissim_final, 
    title = "SNARE Data Cell Types\nSIDEREF (150C, 150G)",
    hclust = FALSE)

snare_group_p_hclust
snare_group_p_no_hclust

snare_hclust_group_lvls <- levels(snare_group_p_hclust$data$group)

### Group Distance DF
snare_sideref_gp_dist_df <- 
  getDistDf(group_labels = snare_test_cell_types, 
            dist_mat = side_seq_ref_rand$dissim_final,
            gather = TRUE)


### (B) PCA RES ---------------------------------------------------------------

## Load PCA, Euclid, Spearman, Pearson distances for test set.
pca_res <- prcomp(t(snare_test_expr_mat))$x[,1:50]

pca25_dist <- dist(pca_res[,1:25])

pca_mat <- matrix(0, nrow = snare_test_cols,
                  ncol = snare_test_cols)

pca_mat[lower.tri(pca_mat)]  <- c(pca25_dist)

pca_mat[upper.tri(pca_mat)] <- 
  t(pca_mat)[upper.tri(pca_mat)]

pca25_dist <- pca_mat

### Plotting
pca25_group_p_hclust <-
  groupwiseDistanceHeatmap(
    group_labels = snare_test_cell_types, 
    dist_mat = pca25_dist, 
    title = "SNARE Data Cell Types\nPCA (25) => Euclidean Dist\nClustered (by SIDEREF)",
    hclust = FALSE,
    preset_levels = snare_hclust_group_lvls)


### Group Distance DF
snare_pca25_gp_dist_df <- 
  getDistDf(group_labels = snare_test_cell_types, 
            dist_mat = pca25_dist,
            gather = TRUE)

## remove dist mat
rm(pca_mat, pca25_dist)



### (C) Pearson RES -----------------------------------------------------------

pear_dist <- 1 - abs(cor(snare_test_expr_mat))

### Plotting
pear_group_p_hclust <-
  groupwiseDistanceHeatmap(
    group_labels = snare_test_cell_types, 
    dist_mat = pear_dist, 
    title = "SNARE Data Cell Types\nPearson Dist\nClustered (by SIDEREF)",
    hclust = FALSE,
    preset_levels = snare_hclust_group_lvls)


### Group Distance DF
snare_pear_gp_dist_df <- 
  getDistDf(group_labels = snare_test_cell_types, 
            dist_mat = pear_dist,
            gather = TRUE)

## remove dist mat
rm(pear_dist)

### (D) Spearman RES ----------------------------------------------------------


spearman_dist <- 1 - abs(cor(snare_test_expr_mat, method = "spearman"))

### Plotting
spearman_group_p_hclust <-
  groupwiseDistanceHeatmap(
    group_labels = snare_test_cell_types, 
    dist_mat = spearman_dist, 
    title = "SNARE Data Cell Types\nSpearman Dist\nClustered (by SIDEREF)",
    hclust = FALSE,
    preset_levels = snare_hclust_group_lvls)


### Group Distance DF
snare_spearman_gp_dist_df <- 
  getDistDf(group_labels = snare_test_cell_types, 
            dist_mat = spearman_dist,
            gather = TRUE)

## remove dist mat
rm(spearman_group_p_hclust)




###############################################################################
### (2) Exclusion of Vip/Lamp5 Astro on Vis -----------------------------------

## get Vip/Lamp5, Astro cell indices.
cell_indices_remove = 
  which(grepl("Vip|Astro", snare_test_cell_types))

cell_indices_keep =
  setdiff(seq_len(snare_test_cols), cell_indices_remove)

snare_group_p_hclust_no_astro_vip <-
  groupwiseDistanceHeatmap(
    group_labels = as.factor(as.character(snare_test_cell_types[cell_indices_keep])), 
    dist_mat = side_seq_ref_rand$dissim_final[cell_indices_keep, cell_indices_keep], 
    title = "SNARE Data Cell Types\nSIDEREF (150C, 150G)\nClustered",
    hclust = TRUE)


### some non-modular work to create reduced heatmap on full data group distance values:

no_atro_vip_lvls <- levels(snare_group_p_hclust_no_astro_vip$data$group)

n_groups = length(no_atro_vip_lvls)

snare_sideref_gp_dist_df_no_atro_vip <- 
  snare_sideref_gp_dist_df %>% 
  dplyr::filter(group %in% no_atro_vip_lvls & group2 %in% no_atro_vip_lvls)


snare_sideref_gp_dist_df_no_atro_vip$group = 
  factor(snare_sideref_gp_dist_df_no_atro_vip$group, 
         levels = no_atro_vip_lvls)
snare_sideref_gp_dist_df_no_atro_vip$group2 = 
  factor(snare_sideref_gp_dist_df_no_atro_vip$group2, 
         levels = no_atro_vip_lvls)

snare_group_p_full_trim_astro_vip <-
  snare_sideref_gp_dist_df_no_atro_vip %>% 
  ggplot(aes(x=group2, y=group, fill=dist)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "royalblue4", high = "orangered2", mid = "white", 
                       midpoint = 0.5, limit = c(0,1)) + 
  labs(x = "Source", y = "Target", fill = "Group-Wise\nDistance") + 
  ggtitle("SNARE Data Cell Types\nSIDEREF (150C, 150G)\nClustered")




###############################################################################
### (3) SNARESEQ Group Ranking Comparisons ------------------------------------



dist_compare_list <- 
  list("PCA25" = snare_pca25_gp_dist_df$dist,
       "Pearson" = snare_pear_gp_dist_df$dist,
       "Spearman" = snare_spearman_gp_dist_df$dist)

## get rank results df
rank_res <-
  lapply(seq_len(length(dist_compare_list)),
         function(x) {
           df =
             computeRankAccuracy(group_labels = snare_test_cell_types,
                                 dist_mat = side_seq_ref_rand$dissim_final, 
                                 ground_truth = dist_compare_list[[x]],
                                 dist_name = names(dist_compare_list)[x])
           df = data.frame(group = df$group_name_vec,
                           groupwise_rank_acc        = df$groupwise_rank_acc,
                           groupwise_mean_rank_error = df$groupwise_mean_rank_error,
                           overall_rank_acc          = df$overall_rank_acc,
                           overall_mean_rank_error   = df$overall_mean_rank_error,
                           method = df$dist_name)
         })

rank_res <- do.call("rbind", rank_res)

## bar plots


## groupwise
rank_res %>% 
  ggplot(aes(x=group, y=groupwise_rank_acc, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="", y="Rank Accuracy", fill = "Method") + 
  theme_bw() +
  ggtitle("Group Rank Consistency with SIDEREF")


rank_res %>% 
  ggplot(aes(x=group, y=groupwise_mean_rank_error, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="", y="Mean Absolute Rank Error", fill = "Method") + 
  theme_bw() +
  ggtitle("Group Rank Consistency with SIDEREF")


## overall
rank_res %>% 
  distinct(overall_rank_acc, overall_mean_rank_error, method) %>%
  mutate(group = "overall") %>%
  ggplot(aes(x=group, y=overall_rank_acc, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="", y="Mean Rank Accuracy", fill = "Method") + 
  theme_bw() +
  ggtitle("Overall Group Rank Consistency with SIDEREF")

rank_res %>% 
  distinct(overall_rank_acc, overall_mean_rank_error, method) %>%
  mutate(group = "overall") %>%
  ggplot(aes(x=group, y=overall_mean_rank_error, fill=method)) + 
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="", y="Mean Absolute Rank Error", fill = "Method") + 
  theme_bw() +
  ggtitle("Overall Group Rank Consistency with SIDEREF")

