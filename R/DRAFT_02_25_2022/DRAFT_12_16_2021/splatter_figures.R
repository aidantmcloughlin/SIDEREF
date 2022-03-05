###############################################################################
### Final Splatter Simulation Results Figures

library(here)
source(here("R/libraries.R"))


######### Constants for plotting purposes =====================================

### Output file paths:
splat_hclust_out_path <- 
  here("output/splatter_sim/hclust_res/")

splat_hclust_out_files <- 
  list.files(splat_hclust_out_path)


core_names <- 
  c("pca_3", "pca_10", "pca_25",
    "euclid", "pear", "spearman",
    "side_ref_g50", "side_ref_g150", "side_ref_g300")

plot_names <- 
  c("PCA (3)", "PCA (10)", "PCA (25)", 
    "Euclidean", "Pearson", "Spearman", 
    "SIDE REF (50)", "SIDE REF (150)", "SIDE REF (300)")

names_map <- 
  data.frame(core_name = core_names,
             plot_name = plot_names)


splatter_subgroup_map <-
  data.frame(group_num = levels(results_list$meta_data$SubGroup)) %>% 
  mutate(cell_group_high = row_number()) %>% 
  mutate(subgroup_label = c(
    "[1] High Ind. D.E. Prob",
    "[2] Low Ind. D.E. Prob",
    "[3,4,5] High Ind. D.E., Low Shared D.E.",
    "[6,7,8] High Ind. D.E., High Shared D.E.",
    "[9,10,11] Low Ind D.E., High Shared D.E.",
    "[12, 13, 14] Low Ind D.E., Low Shared D.E.",
    "[15, 16, 17] Low Ind D.E., Low Shared D.E., High D.E. Variance",
    "[18, 19, 20] High Ind D.E., Low Shared D.E., High D.E. Variance"
  ))

splatter_new_group_map <- 
  data.frame(
    cell_group_low = 1:20,
    new_cell_group_low = c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                           2, 
                           12, 13, 14, 15, 16, 17, 18, 19, 20)
  )


### Heatmaps ==================================================================

SPLAT_FILE_NAME <- 
  "splatter_sim_skin_eq_grps_lowvar_rep1_hclustres_3000g_5_10_25_spect.RData"

## Load results
load(here("output/splatter_sim/hclust_res/" %p% SPLAT_FILE_NAME),
     verbose = TRUE)

df_group_list <- 
  data.frame(cell_group_low = 
               as.numeric(gsub("Group", "", 
                               as.character(results_list$meta_data$Group)))) %>% 
  left_join(splatter_new_group_map)

group_subgroup_map <- 
  results_list$hclust_res_celltypes[[1]]$hclust_res$clust_res_df %>%
  dplyr::distinct(cell_group_low, cell_group_high) %>% 
  left_join(splatter_new_group_map) %>% 
  left_join(splatter_subgroup_map)

X_SIZE = 5
Y_SIZE = 6
TITLE_SIZE = 7

p1 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$side_ref_g300_dist,
  title = "(a) SIDEREF",
  hclust = TRUE,
  preset_levels = as.character(1:20)) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))

p2 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$pca_25_dist,
  title = "(b) PCA (25 Dims.)",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))

p3 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$side_ref_g300_dist_spectral_d5_dist,
  title = "(c) SIDEREF Spectral (5 Dims.)",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))


p4 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$pca_wtd25_dist_spectral_d5_dist,
  title = "(d) PCA (25 Dims.), Spectral (5 Dims.)",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))

p5 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$side_ref_g300_dist_spectral_d10_dist,
  title = "(e) SIDEREF Spectral (10 Dims.)",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))


p6 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$pca_wtd25_dist_spectral_d10_dist,
  title = "(f) PCA (25 Dims.), Spectral (10 Dims.)",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))



ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, 
          common.legend = TRUE, legend="bottom")



p7 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$euclid_dist,
  title = "(b) Euclidean",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))

p8 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$pear_dist,
  title = "(c) 1 - |Pearson Correlation|",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))

p9 <- groupwiseDistanceHeatmap(
  gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
  results_list$dist_res$dist_list$spearman_dist,
  title = "(d) 1 - |Spearman Correlation|",
  hclust = TRUE,
  preset_levels = as.character(1:20)) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE))



ggarrange(p1, p7, p8, p9, nrow = 2, ncol = 2,
          common.legend = TRUE, legend="bottom")


