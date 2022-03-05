################################################################################
### Title:      GeneFishing examination on SIDEREF clusters    
### Author:       Aidan McLoughlin
################################################################################

## some modules currently specific to this script.
library(GeneFishing)
source(here("R/subgroup_composition_gene_contrib.R"))



###############################################################################
###  Continued modules on splatter 2 distance objective functions
###############################################################################



### Simulation 2: Shared DE groups, load process
### Simulation 3: Noisier

## load simulation data
RERUN_SPLATTER_DIST_RES <- FALSE
ADD_SPECTRAL_EMBEDDING <- TRUE
SPECTRAL_DIMS <- 10
main_file_name = "splatter_sim_skin_eq_grps_lowvar"
sim_df_file = main_file_name %p% ".RData"
dist_output_file = main_file_name %p% "_dist_res.RData"



loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
sim_list <- loadRData(here("output/" %p% sim_df_file))

sim <- sim_list$sim

sim <- logNormCounts(sim)
sim_counts <- as.matrix(sim@assays@data$logcounts)
## remove genes with no expression
sim_counts <- sim_counts[rowSums(sim_counts) > 0, ]

splat_group_labels <- factor(sim@colData@listData$Group)





if(RERUN_SPLATTER_DIST_RES) {
  ## get full dist, umap, plot results
  full_res <- runDimReds(sim_counts   = sim_counts,
                         group_labels = splat_group_labels, 
                         r = 50)
  
  ## save
  save(full_res, file = here("output/" %p% dist_output_file)) 
} else{
  load(here("output/" %p% dist_output_file))
}


if(ADD_SPECTRAL_EMBEDDING) {
  
  for(dist_name in names(full_res$dist_list[!grepl("pca", names(full_res$dist_list))])) {
    print(dist_name)
    n <- ncol(full_res$dist_list[[1]])
    
    l <- length(full_res$dist_list)
    
    dissim <- full_res$dist_list[[dist_name]]
    
    dissim <- dist(spectralEmbed(dissim, ndim = SPECTRAL_DIMS))
    
    mat <- matrix(0, nrow = n, ncol = n)
    
    mat[lower.tri(mat)]  <- c(dissim)
    
    mat[upper.tri(mat)] <- 
      t(mat)[upper.tri(mat)]
    
    full_res$dist_list[[l+1]] <- mat
    names(full_res$dist_list)[l+1] <- dist_name %p% "_spectral_d" %p% SPECTRAL_DIMS %p% "_dist"
    }
  }




### Compute Group Differences.


row_data_sim <- data.frame(rowData(sim)) %>% 
  dplyr::select(-(Gene:GeneMean))

groups <- seq_len(length(unique(splat_group_labels)))
n_groups <- length(groups)

de_stats <- data.frame(expand.grid(to_group = groups, from_group = groups))

n_diff_genes <- 
  apply(de_stats, 1, function(x) {
    sum((abs(1 - row_data_sim[,x[1]]) > 0.001 & 
           abs(1 - row_data_sim[,x[2]]) <= 0.001) |
          (abs(1 - row_data_sim[,x[1]]) <= 0.001 & 
             abs(1 - row_data_sim[,x[2]]) > 0.001))
    })

de_stats <- de_stats %>% 
  mutate(total_gens   = nrow(row_data_sim),
         n_diff_genes = n_diff_genes)

de_stats_test <- de_stats %>% 
  group_by(from_group) %>% 
  mutate(dist_norm = (n_diff_genes - min(n_diff_genes)) / 
           (max(n_diff_genes) - min(n_diff_genes)))




ground_truth_heat_plot <- 
  de_stats_test %>%
  mutate(from_group = factor(as.character(from_group),
                             levels = as.character(groups)),
         to_group   = factor(as.character(to_group),
                             levels = as.character(groups))) %>% 
  ggplot(aes(x = from_group, y = to_group, fill = dist_norm)) +
  geom_tile() + 
  scale_fill_gradient2(low = "royalblue4", high = "orangered2", mid = "white", 
                       midpoint = 0.5, limit = c(0,1)) + 
  labs(x = "Source", y = "Target", fill = "Group-Wise\nDistance") + 
  ggtitle("Ground Truth") + 
  theme_bw()

### gather ground truth subgroups ---------------------------------------------

splat_subgroups <- 
  lapply(names(row_data_sim)[grepl("shared", names(row_data_sim))],
         function(x){
           x = gsub("shared_genes", "", x)
           x = strsplit(x, split='_')[[1]]
           return(as.numeric(x))
         })


head(splat_group_labels)

splat_subg_inds <- rep(NA, length(splat_group_labels))

for(i in seq_len(length(splat_subgroups))) {
  splat_subg_inds <- 
    ifelse(is.na(splat_subg_inds) & 
             as.numeric(splat_group_labels) %in% splat_subgroups[[i]], 
           i, 
           splat_subg_inds)
}


## include labels for the individual groups
splat_subgs_inds_w_singletons <- splat_subg_inds
ind_groups <- setdiff(seq_len(n_groups), do.call(c, splat_subgroups))
ind_subg_labs <- max(splat_subg_inds, na.rm=TRUE)+seq_len(length(ind_groups))

for(i in seq_len(length(ind_groups))) {
  ind_cells = which(as.numeric(splat_group_labels) == ind_groups[i])
  splat_subgs_inds_w_singletons[ind_cells] = ind_subg_labs[i]
}

color_mapping <- 
  seq_len(n_groups)

g = n_groups+1
for(x in splat_subgroups){
  color_mapping[x] <- g
  g = g+1
}
rm(g)

splat_color_mapping = as.numeric(as.factor(color_mapping))

###############################################################################
### Compute the relevant statistics and plot
###############################################################################


### Dynamic titling of distance results.
gg_names <-
  c("SIDEREF (50 Genes)",
    "SIDEREF (150 Genes)",
    "SIDEREF (300 Genes)",
    "PCA (3 PCs)",
    "PCA (25 PCs)",
    "PCA (75 PCs)",
    "Euclidean",
    "Pearson",
    "Spearman")

if(ADD_SPECTRAL_EMBEDDING) {
  gg_names <-
    c(gg_names,
      gg_names[!grepl("PCA", gg_names)] %p%
        " Spectral (" %p% SPECTRAL_DIMS %p% " Dims)")
}


method_levels <-
  c(gg_names[grepl("SIDEREF", gg_names)],
    gg_names[!grepl("SIDEREF", gg_names)])



name_xwalk <-
  data.frame(ggplot_name = gg_names,
             method = names(full_res$dist_list),
             stringsAsFactors = FALSE)


# ######### UMAP PLOTS ----------------------------------------------------------
# 
# umap_plot_list <- 
#   lapply(1:length(full_res$dist_list), 
#          function(i) {
#            splatUMAPPlot(
#              full_res$dist_list[[i]],
#              title = name_xwalk$ggplot_name[
#                which(name_xwalk$method == names(full_res$dist_list)[i])],
#              celltypes = splat_group_labels,
#              color_mapping = splat_color_mapping)
#          })
# 
# names(umap_plot_list) <- name_xwalk$method
# 
# 
# 
# 
# 
# do.call(ggarrange, 
#         c(umap_plot_list[!grepl("spectral", names(umap_plot_list), ignore.case = TRUE)],
#           nrow=2, ncol = 5,
#           legend = "bottom", common.legend = TRUE))
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% "_umap.png"), 
#        width = 18, height = 9)
# 
# 
# ### spectral measures
# do.call(ggarrange, 
#         c(umap_plot_list[grepl("spectral", names(umap_plot_list), ignore.case = TRUE)],
#           nrow=2, ncol = 4,
#           legend = "bottom", common.legend = TRUE))
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% "_spectral_umap.png"), 
#        width = 15, height = 9)
# 
# 
# 
# 
# ### Create Base Heatmaps for each Distance Measure
# 
# plot_list <- 
#   lapply(1:length(full_res$dist_list), 
#          function(i) {
#            groupwiseDistanceHeatmap(
#              group_labels = gsub("Group", "", splat_group_labels), 
#              dist_mat = full_res$dist_list[[i]], 
#              title = name_xwalk$ggplot_name[which(name_xwalk$method == 
#                                                     names(full_res$dist_list)[i])],
#              preset_levels = as.character(groups))
#          })
# 
# 
# ## Overlay without ground truth
# 
# ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
#           plot_list[[4]], plot_list[[5]], plot_list[[6]], 
#           plot_list[[7]], plot_list[[8]], plot_list[[9]], 
#           ncol=3, nrow=3)
# 
# # ggsave(here("output/figures/dist_heatmaps.png"), 
# #        width = 15.5, height = 9)
# 
# ### Overlay with ground truth
# 
# ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
#           plot_list[[7]], ground_truth_heat_plot, plot_list[[9]], 
#           plot_list[[4]], plot_list[[5]], plot_list[[8]], 
#           ncol=3, nrow=3)
# 
# 
# ggsave(here("output/figures/" %p% 
#               main_file_name %p% 
#               "_gt_dist_heatmaps.png"), 
#        width = 15.5, height = 9)
# 
# 
# ### Rank accuracy res --------------------------
# 
# 
# ## TODO: fix this for now group_labels arg.
# rank_res <-
#   lapply(seq_len(length(full_res$dist_list)),
#          function(x) {
#            rank_df =
#              computeRankAccuracy(group_labels = splat_group_labels,
#                                  dist_mat = full_res$dist_list[[x]],
#                                  ground_truth = de_stats_test$dist_norm,
#                                  dist_name = names(full_res$dist_list)[x])
#            rank_df = data.frame(group = factor(as.character(groups), 
#                                                levels = as.character(groups)),
#                            groupwise_spearman        = rank_df$groupwise_spearman,
#                            groupwise_rank_acc        = rank_df$groupwise_rank_acc,
#                            groupwise_mean_rank_error = rank_df$groupwise_mean_rank_error,
#                            overall_mean_spearman     = rank_df$overall_mean_spearman,
#                            overall_rank_acc          = rank_df$overall_rank_acc,
#                            overall_mean_rank_error   = rank_df$overall_mean_rank_error,
#                            method = rank_df$dist_name)
#          })
# 
# rank_res <- do.call("rbind", rank_res) %>% 
#   left_join(name_xwalk)
# 
# 
# 
# rank_res$ggplot_name <-
#   forcats::fct_relevel(rank_res$ggplot_name,
#                        levels = method_levels)
# 
# 
# 
# ### Spearman Results ----------------------------------------------------------
# 
# ### Groupwise results:
# rank_res %>% 
#   dplyr::filter(method != "pca_dist_3") %>% 
#   dplyr::filter(group != "all") %>%
#   dplyr::filter(grepl("spectral", method, ignore.case = TRUE)) %>%
#   ggplot(aes(x=group, y=groupwise_spearman, fill=ggplot_name)) + 
#   geom_bar(stat="identity", position = "dodge") + 
#   labs(x="", y="Spearman Cor.", fill = "Method") + 
#   theme_bw() + 
#   theme(legend.position = "bottom") + 
#   ylim(0, 1)
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% 
#               "spearman_spectral_" %p% ".png"),
#        width = 20, height = 5)
# 
# rank_res %>% 
#   dplyr::filter(method != "pca_dist_3") %>% 
#   dplyr::filter(group != "all") %>%
#   dplyr::filter(!grepl("spectral", method, ignore.case = TRUE)) %>%
#   ggplot(aes(x=group, y=groupwise_spearman, fill=ggplot_name)) + 
#   geom_bar(stat="identity", position = "dodge") + 
#   labs(x="", y="Spearman Cor.", fill = "Method") + 
#   theme_bw() + 
#   theme(legend.position = "bottom") + 
#   ylim(0, 1)
#   
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% 
#               "spearman_no_spectral_" %p% ".png"),
#        width = 20, height = 5)
# 
# 
# ### Overall results:
# rank_res %>% 
#   #dplyr::filter(method != "pca_dist_3") %>% 
#   distinct(overall_mean_spearman, overall_rank_acc, overall_mean_rank_error, 
#            method, ggplot_name) %>%
#   mutate(group = "overall") %>%
#   ggplot(aes(x=group, y=overall_mean_spearman, fill=ggplot_name)) + 
#   geom_bar(stat="identity", position = "dodge") + 
#   labs(x="", y="Mean Spearman Cor.", fill = "Method") + 
#   theme_bw() + 
#   theme(legend.position = "bottom") + 
#   ylim(0, 1)
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% 
#               "spearman_overall_" %p% ".png"),
#        width = 12, height = 6)
# 
# 
# ## TODO: remove the self-ranking from Spearman computation? (if time....)
# 
# 
# ### Subgroup Silhouette Widths ---------------------------------
# 
# subgroup_sil_width_res <- 
#   lapply(seq_len(length(full_res$dist_list)),
#          function(x) {
#            subgroupSilWidth(group_indices = as.numeric(splat_group_labels),
#                                subgroup_indices = splat_subg_inds,
#                                dist_mat = full_res$dist_list[[x]],
#                                dist_name = names(full_res$dist_list)[x],
#                                ind_stat = "mean",
#                             group_stat = "mean")
#          })
# 
# subgroup_sil_width_res <- 
#   do.call("rbind", subgroup_sil_width_res) %>% 
#   left_join(name_xwalk)
# 
# 
# subgroup_sil_width_res$ggplot_name <-
#   forcats::fct_relevel(subgroup_sil_width_res$ggplot_name,
#                        levels = method_levels)
# 
# 
# subgroup_sil_width_res <- 
#   subgroup_sil_width_res %>% 
#   mutate(group = gsub("group", "", group),
#          group = factor(group,
#                         levels = c(as.character(seq_len(n_groups)),
#                                    "all")))
# 
# subgroup_sil_width_res %>% 
#   dplyr::filter(method != "pca_dist_3") %>% 
#   dplyr::filter(grepl("spectral", method, ignore.case = TRUE)) %>% 
#   dplyr::filter(group != "all") %>%
#   ggplot(aes(x=group, y=sil_width, fill=ggplot_name)) + 
#   geom_bar(stat="identity", position = "dodge") + 
#   labs(x="", y="Mean Sil. Width", fill = "Method") + 
#   theme_bw() +
#   theme(legend.position = 'bottom') +
#   ylim(-.25, 0.8)
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% 
#               "silwidth_spectral_" %p% ".png"),
#        width = 20, height = 5)
# 
# subgroup_sil_width_res %>% 
#   dplyr::filter(method != "pca_dist_3") %>% 
#   dplyr::filter(!grepl("spectral", method, ignore.case = TRUE)) %>% 
#   dplyr::filter(group != "all") %>%
#   ggplot(aes(x=group, y=sil_width, fill=ggplot_name)) + 
#   geom_bar(stat="identity", position = "dodge") + 
#   labs(x="", y="Mean Sil. Width", fill = "Method") + 
#   theme_bw() + 
#   theme(legend.position = 'bottom') + 
#   ylim(-.25, 0.8)
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% 
#               "silwidth_no_spectral_" %p% ".png"),
#        width = 20, height = 5)
# 
# 
# subgroup_sil_width_res %>% 
#   dplyr::filter(group == "all" & method != "pca_dist_3") %>%
#   ggplot(aes(x=group, y=sil_width, fill=ggplot_name)) + 
#   geom_bar(stat="identity", position = "dodge") + 
#   labs(x="", y="Mean Sil. Width", fill = "Method") + 
#   theme_bw() +
#   theme(legend.position = "bottom")
# 
# ggsave(here("output/figures/subgs/" %p% main_file_name %p% 
#               "silwidth_overall_" %p% ".png"),
#        width = 12, height = 6)
# 
# 
# 
# 
# ###############################################################################
# ###  GeneFishing Module
# ###############################################################################
# 
# 
# ##    (4) Relevant G.O. bait list and include cell groups based on distance.
# 
# 
# brain_go_term_summary <- 
#   read_xlsx(path = here("data/gene_ontology/" %p% 
#                           "GO_term_summary_20210218_185754.xlsx")) %>%
#   clean_names()
# 
# ## Extract highly expressed oligo genes
# 
# ## Oligodenrocytes help maintain and generate the myelin sheath that surrounds axons
# myelin_genes <- which(grepl("myelin", brain_go_term_summary$annotated_term)|
#                                         grepl("myelin", brain_go_term_summary$context))
# 
# myelin_gene_names <- unique(brain_go_term_summary$symbol[myelin_genes])
# 
# myelin_go_data <- brain_go_term_summary[myelin_genes,]
# 
# 
# sum(myelin_gene_names %in% rownames(snare_test_expr_mat))
# 
# #Do we get more myelin genes with more variable features included?? 
# DefaultAssay(snare_test) <- "RNA"
# snare_test <- FindVariableFeatures(snare_test, nfeatures = 2500)
# snare_test <- NormalizeData(snare_test)
# snare_test <- ScaleData(snare_test)
# 
# 
# min(rowVars(snare_test@assays$RNA@scale.data))
# 
# sum(myelin_gene_names %in% rownames(snare_test@assays$RNA@scale.data))
# 
# min(rowVars(snare_test@assays$RNA@scale.data))
# 
# 
# 
# ### (1) Reasonably Variable Oligo Genes and Locations of Oligo Cells
# myelin_gene_names_filt = 
#   myelin_gene_names[which(myelin_gene_names %in% 
#                             rownames(snare_test@assays$RNA@scale.data))]
# 
# oligos <- which(snare_test_cell_types == "Oligo")
# 
# 
# ### (2) Get SIDEREF group distances of Oligo cells
# 
# ## "similar cells"
# snare_sideref_gp_dist_df %>% 
#   dplyr::filter(group2 == "Oligo") %>% 
#   arrange(dist) %>% pull(group)
# 
# 
# ## "distal cells"
# snare_sideref_gp_dist_df %>% 
#   dplyr::filter(group2 == "Oligo") %>% 
#   arrange(-dist) %>% pull(group)
# 
# 
# ## let's add the Astro cells
# astros <- which(snare_test_cell_types == "Astro")
# 
# ## Remove Zero Variance Genes
# 
# snare_gf_final = snare_test@assays$RNA@scale.data[, c(oligos, astros)]
# dim(snare_gf_final)
# snare_gf_final = snare_gf_final[which(rowVars(snare_gf_final)>1e-8),]
# dim(snare_gf_final)
# 
# ## did we retain all the myelin genes?
# myelin_gene_names_filt %in% rownames(snare_gf_final)
# summary(rowVars(snare_gf_final))
# 
# ## Load Modules
# #gf_modules <- list.files(here("../GeneFishing/R/"))
# 
# # for(i in seq_len(length(gf_modules))) {
# #   source(here("../GeneFishing/R/" %p% gf_modules[i]))
# # }
# # 
# # library(doParallel)
# # library(doSNOW)
# library(GeneFishing)
# 
# 
# chosen_bait_genes <- 
#   probeFishability(X = snare_gf_final, 
#                    potential_bait = myelin_gene_names_filt, 
#                    n_rounds = 50, 
#                    min_tightness = 0.5, 
#                    alpha = 5,
#                    umap = TRUE,
#                    ncores = 6)
# 
# 
# discovered_bait <- 
#   geneFishing(X = snare_gf_final, 
#               bait_genes = myelin_gene_names_filt, 
#               k = 2, 
#               umap = TRUE,
#               alpha = 5,
#               fishing_rounds = 500,
#               n_probing_rounds = 50,
#               ## The default here is 0.5 so probably a poor setting.
#               min_tightness = 2.5,
#               ncores = 6)
# 
# 
# 
# 
# 
# ###############################################################################
# ### TODO: Gene Contribution Module
# ###############################################################################
# 
# 
# 
# 
