library(here)
source(here("R/libraries.R"))
source(here("R/tm_load_cut.R"))
source(here('R/relative_group_dist_comps.R'))

set.seed(1)

### TM Helper Funcs
source(here("R/tm_helper_functions.R"))

### Load results from full tissue distance compute.
loadTMCellTypeSampleDistFile(
  file_tag = "main_dist_comps",
  file_loc = here("output/tab_muris_sc/dist_res/"),
  samp_size = TM_SAMP_SIZE,
  use_droplet_only = USE_DROPLET_ONLY,
  load_meta_data = TRUE)

meta_data_full <- 
  prepTMGroupLabels(meta_data = meta_data_full, include_tissue = TRUE)

###############################################################################
### Immune Tissue-Only Heatmaps
###############################################################################

### Displaying all 80 cell types on the heatmap is unwieldy. In practice, one 
##   has a subset of cell groups for which the global comparison will be useful, 
##   such as broadly defined cell types compared to more narrowly defined cell 
##   types.

## B Cells and Leukocytes
b_cells     <- meta_data_full$cell_id[
  which(meta_data_full$cell_type == "B cell")]
leuko_cells <- meta_data_full$cell_id[
  which(meta_data_full$cell_type == "leukocyte")]

## Global group of endothelial cells.
endo_selected_cells <- meta_data_full %>% 
  dplyr::filter((tissue == "Limb_Muscle" & cell_type == "endothelial cell") |
                  (tissue == "Lung" & cell_type == "lung endothelial cell") | 
                  (tissue == "Heart_and_Aorta" & cell_type == "endocardial cell")) %>%
  pull(cell_id)

## a couple unrelated cell types to provide "far end" of global distance gradient.
unrelated_cells <- meta_data_full %>%
  dplyr::filter((tissue == "Lung" & cell_type == "type II pneumocyte") |
                  (tissue == "Tongue" & cell_type == "basal cell of epidermis")) %>% 
  pull(cell_id)


cells_wout_b <- union(leuko_cells, 
                      union(endo_selected_cells,
                            unrelated_2_cells))

cells_wout_leuko <- union(b_cells, 
                          union(endo_selected_cells,
                                unrelated_2_cells))


cells_w_b <- union(b_cells,
                   cells_wout_b)




############################### ===============================================
### Label Ordering for Heatmaps ===============================================
############################### ===============================================


cells_wout_leuko_levels <- 
  c("Type II Pneumocyte, Lung", 
    "Basal Cell of Epidermis, Tongue",
    "Endocardial Cell, Heart and Aorta",
    "Endothelial Cell, Limb Muscle", 
    "Lung Endothelial Cell, Lung", 
    "B Cell, Spleen", 
    "B Cell, Lung", 
    "B Cell, Limb Muscle", 
    "B Cell, Mammary Gland")

cells_wout_b_levels <- 
  c("Type II Pneumocyte, Lung",
    "Basal Cell of Epidermis, Tongue", 
    "Endocardial Cell, Heart and Aorta",
    "Endothelial Cell, Limb Muscle", 
    "Lung Endothelial Cell, Lung", 
    "Leukocyte, Liver", 
    "Leukocyte, Bladder",
    "Leukocyte, Thymus",
    "Leukocyte, Kidney",
    "Leukocyte, Lung"
  )

cells_w_b_levels <- 
  c("Type II Pneumocyte, Lung",
    "Basal Cell of Epidermis, Tongue", 
    "Endocardial Cell, Heart and Aorta",
    "Endothelial Cell, Limb Muscle", 
    "Lung Endothelial Cell, Lung", 
    "Leukocyte, Liver", 
    "Leukocyte, Bladder",
    "Leukocyte, Lung",
    "Leukocyte, Thymus",
    "Leukocyte, Kidney",
    "B Cell, Spleen", 
    "B Cell, Lung", 
    "B Cell, Limb Muscle", 
    "B Cell, Mammary Gland")



###############################################################################
### Functions to overlay 3 heatmap versions with each other.
############################################################################### 



groupwiseDistanceHeatmapCellType <- 
  function(cell_id_vec, dist_mat, meta_data,
           symmetrize = FALSE,
           title = "",
           title_size = 10,
           do_hclust_axes = FALSE,
           preset_levels = NULL,
           do_numeric_x=TRUE,
           do_numeric_y=FALSE) {
    
    rownames(meta_data) = meta_data$cell_id
    
    dist_mat <- selectFullDistIdx(cell_id_vec,
                                  dist = dist_mat,
                                  meta_data = meta_data)
    
    group_labels <- getGroupLabels(meta_data, cell_id_vec = cell_id_vec)
    
    groupwiseDistanceHeatmap(group_labels, 
                             dist_mat,
                             symmetrize = symmetrize,
                             title = title,
                             do_hclust_axes = do_hclust_axes,
                             preset_levels = preset_levels,
                             do_numeric_x = do_numeric_x,
                             do_numeric_y = do_numeric_y)  +
      theme(plot.title = element_text(size=title_size))
  }



# createHeatmapTM <- 
#   function(cell_vec, full_dist_mat,
#            meta_data, 
#            type_levels,
#            symmetrize = FALSE,
#            title = "(a)",
#            title_size = 8,
#            do_numeric_y = FALSE) {
#     
#     
#     ## Truncating the full distance matrix:
#     fig <-
#       groupwiseDistanceHeatmapCellType(
#         cell_vec = cell_vec, 
#         dist_mat = selectFullDistIdx(cell_vec,
#                                      dist = full_dist_mat,
#                                      meta_data = meta_data),
#         ## TODO: is meta_data a useful arg in the package-level function?
#         meta_data = meta_data,
#         #hclust = TRUE,
#         preset_levels = type_levels,
#         symmetrize = symmetrize,
#         do_numeric_x = TRUE,
#         do_numeric_y = do_numeric_y,
#         title = title,
#         title_size = title_size) + 
#       theme_bw() +
#       theme(axis.title.y = element_blank())
#     
#     fig
#     return(fig)
#     
#   }


## TODO: select figures and adjust comments appropriately.

###############################################################################
### Fig NEW A: BiPartite Graphs
###############################################################################

source(here("R/bipartite_graphs.R"))

## TODO: can make this even more modular for looping through the same set of distances per fig type..

bipartiteTMCellsWithB <- function(dist, title) {
  bipartiteNetworkGraph(
    dist_mat = selectFullDistIdx(
      cell_vec = cells_w_b,
      dist = dist,
      meta_data = meta_data_full),
    group_labels = prepTMGroupLabels(cells_w_b, meta_data_full, 
                                     include_tissue = TRUE),
    preset_levels = cells_w_b_levels,
    source_groups = c("Leukocyte, Bladder", "Leukocyte, Kidney"),
    node_size = 7.5,
    text_size = 3,
    title = title,
    legend_title = "Cell Type, Tissue")
  
}

figbi_side <- 
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$side_ref_g300_dist,
    title = "(a) SIDEREF")

figbi_spect3 <-
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    title = "(b) SIDEREF Spectral (3 Dims.)")

figbi_spect10 <-
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$side_ref_g300_dist_spectral_d10_dist,
    title = "(c) SIDEREF Spectral (10 Dims.)")

figbi_pca25 <-
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$pca_wtd25_dist,
    title = "(d) PCA (25 Dims.)")

figbi_euclid <-
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$euclid_dist,
    title = "(e) Euclidean")

figbi_pear <-
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$pear_dist,
    title = "(f) 1 - |Pearson Correlation|")

figbi_spearman <-
  bipartiteTMCellsWithB(
    dist = full_dist_res$dist_list$spearman_dist,
    title = "(g) 1 - |Spearman Correlation|")

p_tm_dist_fbi <-
  ggarrange(figbi_side, figbi_spect3, figbi_spect10, 
            figbi_pca25, figbi_euclid, 
            ggplot() + theme_void(), 
            figbi_pear, figbi_spearman,
            nrow=3,
            ncol=3,
            #widths = c(1, 1, 1),
            common.legend = TRUE,
            legend = "right")

p_tm_dist_fbi

# ggsave(here("manuscript_files/revision_figures/FigureBiPartite.eps"),
#        plot = p_tm_dist_fbi,
#        width = 16.5, height = 9,
#        device='eps', dpi=250)

ggsave(here("manuscript_files/revision_figures/FigureBiPartite2.png"),
       plot = p_tm_dist_fbi,
       dpi=250,
       width = 16.5, height = 9)

###############################################################################
### Fig NEW B: Symmetric Heatmap.
###############################################################################

symmetricTMCellsWithB <- function(dist, title) {
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = dist,
    symmetrize = TRUE,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = title,
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE) + 
    theme(plot.title = element_text(size = FIG_TITLE_SIZE))  
  
}

figsym_side <- 
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$side_ref_g300_dist,
    title = "(a) SIDEREF")

figsym_spect3 <-
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    title = "(b) SIDEREF Spectral (3 Dims.)")

figsym_spect10 <-
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$side_ref_g300_dist_spectral_d10_dist,
    title = "(c) SIDEREF Spectral (10 Dims.)")

figsym_pca25 <-
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$pca_wtd25_dist,
    title = "(d) PCA (25 Dims.)")

figsym_euclid <-
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$euclid_dist,
    title = "(e) Euclidean")

figsym_pear <-
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$pear_dist,
    title = "(f) 1 - |Pearson Correlation|")

figsym_spearman <-
  symmetricTMCellsWithB(
    dist = full_dist_res$dist_list$spearman_dist,
    title = "(g) 1 - |Spearman Correlation|")

p_tm_dist_sym <-
  ggarrange(figsym_side, figsym_spect3, figsym_spect10, 
            figsym_pca25, figsym_euclid, 
            ggplot() + theme_void(), 
            figsym_pear, figsym_spearman,
            nrow=3,
            ncol=3,
            #widths = c(1, 1, 1),
            common.legend = TRUE,
            legend = "right")

p_tm_dist_sym

ggsave(here("manuscript_files/revision_figures/Figure4.png"),
       plot = p_tm_dist_sym,
       width = 21, height = 15,
       dpi=250)
       #device='eps', dpi=250)

###############################################################################
### Fig 1: Leukocytes, Endo, 2 Others
###############################################################################

rownames(meta_data_full) <- meta_data_full$cell_id

FIG_TITLE_SIZE <- 8

fig1_side <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(a) SIDEREF",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE) + 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig1_spect3 <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(b) SIDEREF Spectral (3 Dims.)",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig1_spect10 <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d10_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(c) SIDEREF Spectral (10 Dims.)",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))


fig1_pca <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$pca_wtd25_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(d) PCA (25 Dims.)",
    title_size = FIG_TITLE_SIZE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))


fig1_euclid <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$euclid_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(e) Euclidean",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig1_pear <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$pear_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(f) 1 - |Pearson Correlation|",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig1_spearman <-
  createHeatmapTM(
    cell_vec = cells_wout_b, 
    full_dist_mat = full_dist_res$dist_list$spearman_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_b_levels,
    title = "(g) 1 - |Spearman Correlation|",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))



###############################################################################
### Fig 2: B Cells, Endo, 2 Others
###############################################################################

fig2_side <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(a) SIDEREF",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig2_spect3 <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(b) SIDEREF Spectral (3 Dims.)",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig2_spect10 <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d10_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(c) SIDEREF Spectral (10 Dims.)",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))


fig2_pca <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$pca_wtd25_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(d) PCA (25 Dims.)",
    title_size = FIG_TITLE_SIZE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig2_euclid <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$euclid_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(e) Euclidean",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig2_pear <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$pear_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(f) 1 - |Pearson Correlation|",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig2_spearman <-
  createHeatmapTM(
    cell_vec = cells_wout_leuko, 
    full_dist_mat = full_dist_res$dist_list$spearman_dist,
    meta_data = meta_data_full, 
    type_levels = cells_wout_leuko_levels,
    title = "(g) 1 - |Spearman Correlation|",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))


###############################################################################
### Fig 3: B + Leuko, Endo, 2 Others
###############################################################################

fig3_side <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(a) SIDEREF",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig3_spect3 <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(b) SIDEREF Spectral (3 Dims.)",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig3_spect10 <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d10_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(c) SIDEREF Spectral (10 Dims.)",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))


fig3_pca <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$pca_wtd25_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(d) PCA (25 Dims.)",
    title_size = FIG_TITLE_SIZE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig3_euclid <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$euclid_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(e) Euclidean",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig3_pear <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$pear_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(f) 1 - |Pearson Correlation|",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = FALSE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

fig3_spearman <-
  createHeatmapTM(
    cell_vec = cells_w_b, 
    full_dist_mat = full_dist_res$dist_list$spearman_dist,
    meta_data = meta_data_full, 
    type_levels = cells_w_b_levels,
    title = "(g) 1 - |Spearman Correlation|",
    title_size = FIG_TITLE_SIZE,
    do_numeric_y = TRUE)+ 
  theme(plot.title = element_text(size = FIG_TITLE_SIZE))

p_tm_dist_f1 <-
  ggarrange(fig1_side, fig1_spect3, fig1_spect10, 
            fig1_pca, fig1_euclid, ggplot() + theme_minimal(), 
            fig1_pear, fig1_spearman,
            nrow=3,
            ncol=3,
            widths = c(1.8, 1, 1),
            common.legend = TRUE,
            legend = "bottom")


p_tm_dist_f2 <-
  ggarrange(fig2_side, fig2_spect3, fig2_spect10, 
            fig2_pca, fig2_euclid, ggplot() + theme_minimal(), 
            fig2_pear, fig2_spearman,
            nrow=3,
            ncol=3,
            widths = c(1.85, 1, 1),
            common.legend = TRUE,
            legend = "bottom")


p_tm_dist_f3 <-
  ggarrange(fig3_side, fig3_spect3, fig3_spect10, 
            fig3_pca, fig3_euclid, ggplot() + theme_minimal(), 
            fig3_pear, fig3_spearman,
            nrow=3,
            ncol=3,
            widths = c(1.65, 1, 1),
            common.legend = TRUE,
            legend = "bottom")


ggsave(here("manuscript_files/Figure3.eps"),
       plot = p_tm_dist_f1,
       width = 9, height = 9,
       device='eps', dpi=250)

ggsave(here("manuscript_files/Figure4.eps"),
       plot = p_tm_dist_f3,
       width = 10.2, height = 9.75,
       device='eps', dpi=250)

ggsave(here("manuscript_files/FigureS1.eps"),
       plot = p_tm_dist_f2,
       width = 9, height = 8.75,
       device='eps', dpi=250)


## clear space
rm(list=ls()[grepl("fig", ls())])
rm(full_dist_res)
gc()



###############################################################################
### Gene Set Enrichment for B-Cells within Kidney, Thymus Leukocytes.
###############################################################################



### Relaod TM Data
KEEP_DROPLET <- TRUE
KEEP_FACS <- FALSE
data_pct <- ""
file_name <- "tab_muris_full" %p% ifelse(data_pct == "", "", "_" %p% data_pct)


tm_data_list <- 
  tmLoadCut(file_name = file_name, data_pct = data_pct, 
            sample_size = 20,
            sample_cell_types = TRUE,
            keep_droplet = KEEP_DROPLET, keep_facs = KEEP_FACS,
            sample_pct = FALSE)

## extract objects
for(n in names(tm_data_list)) {
  assign(n, tm_data_list[[n]])
}

data_full <- tm_data_list$data_store_list[[1]][[2]]
meta_data_full <- tm_data_list$data_store_list[[1]][[3]]
rownames(meta_data_full) <- meta_data_full$cell_id

rm(tm_data_list, data_store_list)
gc()


tm_preproc <-
  TabMurisProcSubset(full_data = data_full, 
                     meta_data = meta_data_full, 
                     var_features= VAR_FEATURES,
                     return_seurat = TRUE)[[1]]



## DE Genes among the B-cell groups
getMarkerGSList <- function(assay, cells_in, cells_vs, 
                            cell_marker_list_sizes = 
                              c(10, 50, 100, 300),
                            type_name = "B_markers") {
  cell_marker_res <- 
    FindMarkers(assay,
                cells.1 = cells_in,
                cells.2 = setdiff(cells_vs, cells_in),
                logfc.threshold = 0)
  
  cell_marker_gs_list <- 
    vector(mode = "list", 
           length = length(cell_marker_list_sizes))
  
  names(cell_marker_gs_list) <- type_name %p% "_top" %p% 
    cell_marker_list_sizes
  
  for(i in seq_len(length(cell_marker_list_sizes))) {
    cell_marker_gs_list[[i]] <- 
      rownames(cell_marker_res)[seq_len(cell_marker_list_sizes[i])]
  }
  
  return(cell_marker_gs_list)
}


## Gene set enrichment analysis among the leukocyte groups ====================


## Functions for GSE analysis
prepRanks <- function(de_marker_res, all_gene_names = c()) {
  ## reorganize so that downregulated genes are at the bottom of the list 
  de_marker_res_reorg <- 
    de_marker_res %>%
    mutate(rank = nrow(de_marker_res):1,
           rank = rank * ifelse(avg_log2FC > 0, 1, -1)) %>% 
    arrange(-rank)
  
  de_marker_gene_order <- rownames(de_marker_res_reorg)
  
  de_marker_res_ranks <- de_marker_res_reorg$rank
  names(de_marker_res_ranks) <- de_marker_gene_order 
  
  missing <- setdiff(all_gene_names, names(de_marker_res_ranks))
  
  missing_ranks <- rep(median(de_marker_res_ranks), length(missing))
  names(missing_ranks) <- missing
  
  de_marker_res_ranks <- c(de_marker_res_ranks, missing_ranks)
  de_marker_res_ranks <- sort(de_marker_res_ranks, decreasing = TRUE)
  
  
  return(de_marker_res_ranks)
}


getGSEAScores <- function(marker_gs_list,
                          cell_type_name_vec, cells_vs,
                          assay, meta_data,
                          logfc_thres = 0.01) {
  de_marker_list <- vector(mode="list", length = length(cell_type_name_vec))
  
  names(de_marker_list) <- cell_type_name_vec
  
  
  for(i in seq_len(length(cell_type_name_vec))) {
    print(i)
    
    cells <- meta_data %>% 
      dplyr::filter(cell_type %p% ", " %p% 
                      tissue == cell_type_name_vec[i]) %>% 
      pull(cell_id)
    
    de_marker_list[[i]] <-
      FindMarkers(assay,
                  cells.1 = cells,
                  cells.2 = setdiff(cells_vs, cells),
                  logfc.threshold = logfc_thres)
    
    de_marker_list[[i]] <- 
      prepRanks(de_marker_list[[i]],
                all_gene_names = 
                  marker_gs_list[[length(marker_gs_list)]])
    
    print(length(de_marker_list[[i]]))
    
    ### Run FGSEA function
    de_marker_list[[i]] <- fgsea(pathways = marker_gs_list, 
                                 stats = de_marker_list[[i]],
                                 eps=0)
    
  }
  return(de_marker_list)
}


## Gene Enrichment Scores with Canonical B-Cell Genes ======================
##   for each Leukocyte Group ==============================================


## What does SIDEREF Distance Matrix Indicate to us?

### 1: Kidney, Thymus group composition more similar to B-Cell than
### 2: Lung, Bladder, Liver.

leuko_names <- unique(meta_data_full[leuko_cells, ]$cell_type %p% ", " %p%
                        meta_data_full[leuko_cells, ]$tissue)


############################
### B-marker GS lists

### vs all except leukocytes
b_marker_gs_list_vs_non_leuko <-
  getMarkerGSList(assay = tm_preproc@assays$RNA,
                  cells_in = b_cells, 
                  cells_vs = setdiff(colnames(tm_preproc@assays$RNA),
                                     leuko_cells))


############################
### Compute GSEA scores 

### 3. B-non leuko DE list vs Leuko-non B DE list
gsea_res_vs_all_other_cells <- 
  getGSEAScores(marker_gs_list = b_marker_gs_list_vs_non_leuko,
                cell_type_name_vec = leuko_names, 
                cells_vs = setdiff(meta_data_full$cell_id,
                                   b_cells),
                assay = tm_preproc@assays$RNA, 
                meta_data = meta_data_full)


### Visualize GSEA Results

BarplotGSEA <- function(gsea_res_list,
                        preset_levels = NULL) {
  
  gsea_res_list_clean <- 
    lapply(seq_len(length(gsea_res_list)), function(x) {
      return(gsea_res_list[[x]] %>% 
               mutate(celltype = names(gsea_res_list)[x]))
    })
  
  df <- do.call("rbind", gsea_res_list_clean) %>% 
    mutate(num_de_markers = gsub("^.*(top)(.*)", "Top \\2", pathway))
  
  if(!is.null(preset_levels)) {
    df <-
      df %>% 
      mutate(num_de_markers = 
               factor(num_de_markers,
                      levels = preset_levels))
  }
  p <-
    df %>% 
    mutate(celltype = str_to_title(celltype)) %>%
    ggplot(aes(x = num_de_markers, 
               y = abs(NES),
               fill = celltype)) + 
    geom_bar(stat="identity", position="dodge", color="black") + 
    theme_bw() + 
    labs(y = "Abs. Norm. Enrichment Score", x = "B Cell DE Gene Set", 
         fill = "Cell Type, Organ")
  
  return(p)
  
}

p_gsea_bar <-
  BarplotGSEA(gsea_res_vs_all_other_cells, 
              preset_levels = c("Top 10", 
                                "Top 50", 
                                "Top 100",
                                "Top 300"))


p_gsea_bar

ggsave(here("manuscript_files/Figure5.eps"),
       plot = p_gsea_bar,
       width = 6, height = 4,
       device='eps', dpi=250)

gsea_table <-
  do.call(rbind, 
          lapply(seq_len(length(gsea_res_vs_all_other_cells)),
                 function(i) return(
                   gsea_res_vs_all_other_cells[[i]] %>% 
                     mutate(celltype = names(gsea_res_vs_all_other_cells)[i])
                 ))) %>% 
  dplyr::select(pathway, pval, celltype) %>% 
  arrange(pathway, celltype)


