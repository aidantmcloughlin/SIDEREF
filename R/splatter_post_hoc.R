###############################################################################
### Final Splatter Simulation Results Figures

library(here)
source(here("R/libraries.R"))
source(here('R/relative_group_dist_comps.R'))


### Output file paths:
splat_out_path <- 
  here("output/splatter_sim/dist_res/")

splat_out_files <- 
  list.files(splat_out_path)


## Load results
for(f in splat_out_files) {
  load(splat_out_path %p% f,
       verbose = TRUE)  
}



##### NAME MAPPING ====================================================

### Note: this is manually created dependent on the experiment format:

splatter_global_group_map <-
  data.frame(Global_Group = c("Group1", "Group2",
                              "Global_Group" %p% 1:6)) %>% 
  mutate(cell_group_high = row_number()) %>% 
  mutate(global_group_label = c(
    "A [1] High Ind. D.E. Prob",
    "B [2] Low Ind. D.E. Prob",
    "C [3,4,5] High Ind. D.E., Low Shared D.E.",
    "D [6,7,8] High Ind. D.E., High Shared D.E.",
    "E [9,10,11] Low Ind D.E., High Shared D.E.",
    "F [12,13,14] Low Ind D.E., Low Shared D.E.",
    "G [15,16,17] Low Ind D.E., Low Shared D.E., High D.E. Variance",
    "H [18,19,20] High Ind D.E., Low Shared D.E., High D.E. Variance"
  )) %>% 
  
  mutate(hm_global_group_label = c(
    "A", "B", "C", "D", 
    "E", "F", "G", "H"
  ))

### Heatmaps ==================================================================
meta_data <-
  meta_data %>% 
  mutate(cell_group_low = gsub("Group", "", 
                               as.character(meta_data$Group)) %>% 
           as.numeric()) %>% 
  left_join(splatter_global_group_map)

group_labels <- as.character(meta_data$cell_group_low)


group_global_group_map <- 
  meta_data %>%
  dplyr::distinct(cell_group_low, cell_group_high) %>% 
  left_join(splatter_global_group_map)

X_SIZE = 11
Y_SIZE = 12
LEGEND_SIZE = 8
LEGEND_TEXT_SIZE = 8
TITLE_SIZE = 9.5

SYMMETRIZE = FALSE




p1 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = main_dist_output$dist_list$side_ref_g300_dist,
  title = "(a) SIDEREF",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) +
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 10))



p2 <- groupwiseDistanceHeatmap(
  group_labels = group_labels,
  dist_mat = main_dist_output$dist_list$pca_wtd25_dist,
  title = "(b) PCA (25 Dims.)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))

p3 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = spectral_dist_res$dist_list$
    side_ref_g300_dist_spectral_d5_dist,
  title = "(c) SIDEREF Spectral (5 Dims.)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))

p4 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = spectral_dist_res$dist_list$
    pca_wtd25_dist_spectral_d5_dist,
  title = "(d) PCA (25 Dims.), Spectral (5 Dims.)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))

p5 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = spectral_dist_res$dist_list$
    side_ref_g300_dist_spectral_d10_dist,
  title = "(e) SIDEREF Spectral (10 Dims.)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))

p6 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = spectral_dist_res$dist_list$
    pca_wtd25_dist_spectral_d10_dist,
  title = "(f) PCA (25 Dims.), Spectral (10 Dims.)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))


p7 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = main_dist_output$dist_list$euclid_dist,
  title = "(b) Euclidean",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))

p8 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = main_dist_output$dist_list$pear_dist,
  title = "(c) 1 - |Pearson Correlation|",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE))+
  theme(plot.title = element_text(size = 12))

p9 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = main_dist_output$dist_list$spearman_dist,
  title = "(d) 1 - |Spearman Correlation|",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE)) +
  theme(plot.title = element_text(size = 12))


##### OTHER DISTANCE MEASURES ============================
 
p10 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = othr_dist_output$dist_list$RAFSIL,
  title = "(b) RAFSIL",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE)) +
  theme(plot.title = element_text(size = 12))

p11 <- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = othr_dist_output$dist_list$SIMLR_5_dims,
  title = "(c) SIMLR (5 Connected Components)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE)) +
  theme(plot.title = element_text(size = 12))

p12<- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = othr_dist_output$dist_list$SIMLR_10_dims,
  title = "(d) SIMLR (10 Connected Components)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE)) +
  theme(plot.title = element_text(size = 12))

p13<- groupwiseDistanceHeatmap(
  group_labels = group_labels, 
  dist_mat = othr_dist_output$dist_list$SIMLR_25_dims,
  title = "(e) SIMLR (25 Connected Components)",
  symmetrize = SYMMETRIZE,
  do_hclust_axes = FALSE,
  preset_levels = as.character(1:20),
  global_group_label_df = group_global_group_map) + 
  theme(axis.text.x = element_text(size=X_SIZE),
        axis.text.y = element_text(size=Y_SIZE),
        title = element_text(size=TITLE_SIZE),
        legend.title = element_text(size = LEGEND_SIZE)) + 
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = LEGEND_TEXT_SIZE)) +
  theme(plot.title = element_text(size = 12))




### arrange =========================================
splat_fig1 <-
  ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, 
            common.legend = TRUE, legend="right") 

splat_fig2 <-
  ggarrange(p1, p7, p8, p9, nrow = 2, ncol = 2,
            common.legend = TRUE, legend="right")

splat_fig3 <-
  ggarrange(p1, p10, p11, nrow = 1, ncol = 3,
            common.legend = TRUE, legend="right")

splat_figss4 <-
  ggarrange(p12, p13, nrow = 1, ncol = 2,
            common.legend = TRUE, legend="right")


### save =========================================
ggsave(here("manuscript_files//Figure1.png"),
       plot = splat_fig1 +
         theme(
           panel.background = element_rect(fill = "white"), # bg of the panel
           plot.background = element_rect(fill = "white", color = NA), # bg of the plot
           #panel.grid.major = element_blank(), # get rid of major grid
           #panel.grid.minor = element_blank(), # get rid of minor grid
           legend.background = element_rect(fill = "white"), # get rid of legend bg
           legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
         ),
       width = 13, height = 15,
       device='png', dpi=450)

ggsave(here("manuscript_files/Figure2.png"),
       plot = splat_fig2 +
         theme(
           panel.background = element_rect(fill = "white"), # bg of the panel
           plot.background = element_rect(fill = "white", color = NA), # bg of the plot
           #panel.grid.major = element_blank(), # get rid of major grid
           #panel.grid.minor = element_blank(), # get rid of minor grid
           legend.background = element_rect(fill = "white"), # get rid of legend bg
           legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
         ),
       width = 13, height = 10,
       device='png', dpi=450)


ggsave(here("manuscript_files/Figure3.png"),
       plot = splat_fig3 +
         theme(
           panel.background = element_rect(fill = "white"), # bg of the panel
           plot.background = element_rect(fill = "white", color = NA), # bg of the plot
           #panel.grid.major = element_blank(), # get rid of major grid
           #panel.grid.minor = element_blank(), # get rid of minor grid
           legend.background = element_rect(fill = "white"), # get rid of legend bg
           legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
         ),
       width = 23.5, height = 7,
       device='png', dpi=450)

ggsave(here("manuscript_files/FigureS4.png"),
       plot = splat_figss4 +
         theme(
           panel.background = element_rect(fill = "white"), # bg of the panel
           plot.background = element_rect(fill = "white", color = NA), # bg of the plot
           #panel.grid.major = element_blank(), # get rid of major grid
           #panel.grid.minor = element_blank(), # get rid of minor grid
           legend.background = element_rect(fill = "white"), # get rid of legend bg
           legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
         ),
       width = 16.5, height = 7,
       device='png', dpi=450)



### A few UMAP examples ===============================================


p_umap_sideref <-
  data.frame(main_dist_output$umap_list$side_ref_g300_umap) %>% 
  mutate(global_group = meta_data$global_group_label) %>% 
  ggplot() + 
  geom_point(aes(x=X1, y=X2, colour=global_group),
             size = 0.5) + 
  theme_bw() + 
  labs(x = "UMAP1", y = "UMAP2", colour = "Global Group") + 
  ggtitle("(a) SIDEREF")

p_umap_pca25 <-
  data.frame(main_dist_output$umap_list$pca_wtd25_umap) %>% 
  mutate(global_group = meta_data$global_group_label) %>% 
  ggplot() + 
  geom_point(aes(x=X1, y=X2, colour=global_group),
             size = 0.5) + 
  theme_bw() + 
  labs(x = "UMAP1", y = "UMAP2", colour = "Global Group") +
  ggtitle("(b) PCA (25 Dims.)")

splat_umap_fig <-
  ggarrange(p_umap_sideref, 
            p_umap_pca25, ncol=2, nrow=1, 
            common.legend = TRUE, legend="bottom") 


ggsave(here("manuscript_files/FigureS3.eps"),
       plot = splat_umap_fig,
       width =14, height = 7,
       device='eps', dpi=300)




### SIMLR UMAP ===============================================

# simlr_umap <- 
#   uwot::umap(as.dist(othr_dist_output$SIMLR_10_dims))
# 
# ggplot(data.frame(simlr_umap) %>% 
#          mutate(cell_type = as.character(meta_data$cell_group_high))) + 
#   geom_point(aes(x=X1,y=X2,colour=cell_type))



