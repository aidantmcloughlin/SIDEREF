library(here)
source(here("R/libraries.R"))

source(here("R/tm_load_cut.R"))

source(here('R/dist_and_clust_pipelines.R'))
source(here('R/subgroup_composition_gene_contrib.R'))
source(here('R/hcclust.R'))


############################# CONSTANTS ####################################
## Data subset Constants
KEEP_DROPLET <- TRUE
KEEP_FACS <- FALSE


data_pct <- ""
file_name <- "tab_muris_full" %p% ifelse(data_pct == "", "", "_" %p% data_pct)
print(file_name)
REPS <- 1
FOLDS <- 75
FOLD_SEEDS <- 100:(100+REPS-1)

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

rm(tm_data_list, data_store_list)

## Distance setting constants:
SPECTRAL_DIMS <- 25
G <- c(150, 300)
R <- 100
PCs <- c(3, 10, 25)
N_CLUST <- 20 ## k-means clusters during SIDEREF
VAR_FEATURES <- 3000

## Other:
N_CORES <- 4


## Analysis re-run options
RERUN_CELLTYPE_DIST <- FALSE


###############################################################################
### Immune Tissue-Only Heatmaps
###############################################################################

getCellTypeSampleFileString <- function(file_tag) {
  return(here("output/tab_muris_sc/dist_res/full_dist_res__" %p% 
                file_tag %p% "_fold1" %p%
                "_samp_ctTRUE" %p%
                "_samp_pctFALSE" %p%
                "_samp_size20" %p% ".RData"))
}

### Load results from full tissue distance compute.
load(file = getCellTypeSampleFileString("droplet_only"),
     verbose = TRUE) 


### Function to collect random cell types outside the set of interest
sampleRandomCellType <- function(n, exclude_cells = c(),
                                 meta_data) {
  
  types = meta_data$cell_type
  types_tissues = meta_data$cell_type %p% meta_data$tissue
  
  types_cant_sample <- 
    unique(types[meta_data$cell_id %in% exclude_cells])
  
  types_cant_sample_w_tiss <- 
    expand.grid(types_cant_sample, unique(meta_data$tissue))
  
  types_cant_sample_w_tiss <- 
    types_cant_sample_w_tiss$Var1 %p% types_cant_sample_w_tiss$Var2
  
  types_tissues_can_sample <- 
    setdiff(types_tissues, types_cant_sample_w_tiss)
  
  types_tissues_can_sample <- 
    types_tissues_can_sample[!grepl("NA", types_tissues_can_sample)]
  
  return(sample(types_tissues_can_sample, n, replace=FALSE))
}

getCellIdsRandomCellTypes <- function(n, exclude_cell_ids) {
  return(
    meta_data_full$cell_id[
      which(meta_data_full$cell_type %p% 
              meta_data_full$tissue %in% 
              sampleRandomCellType(n, exclude_cells = exclude_cell_ids, 
                                   meta_data = meta_data_full))]
  )
}

### Collect cell groups of interest:
marrow_cells             <- meta_data_full$cell_id[which(meta_data_full$tissue == "Marrow")]
b_cells                  <- meta_data_full$cell_id[which(meta_data_full$cell_type == "B cell")]
leuko_cells              <- meta_data_full$cell_id[which(meta_data_full$cell_type == "leukocyte")]
immune_w_b_cells         <- union(immune_cells, b_cells)
immune_w_leuko_cells     <- union(immune_cells, leuko_cells) 
immune_w_b_w_leuko_cells <- union(immune_w_b_cells, leuko_cells)

endo_cells <- 
  meta_data_full$cell_id[which(grepl("endo", meta_data_full$cell_type, 
                                     ignore.case=TRUE) & 
                                 !grepl("endocrine", meta_data_full$cell_type,
                                        ignore.case=TRUE))]


immune_cells_regex <- "B cell|T cell|macrophage|leukocyte|monocyte|" %p%
  "erythro|hemato|granulo|basophil|blood|mast|natural killer"
immune_cells <- meta_data_full$cell_id[which(grepl(immune_cells_regex, 
                                                   meta_data_full$cell_type))]


## Random sampling:
set.seed(5)
no_immune_no_endo_3 <- getCellIdsRandomCellTypes(3, union(immune_cells,
                                                          endo_cells))
no_immune_5  <- getCellIdsRandomCellTypes(5, immune_cells)
no_b_5       <- getCellIdsRandomCellTypes(5, b_cells)
no_leuko_5    <- getCellIdsRandomCellTypes(5, leuko_cells)

no_immune_10 <- 
  c(no_immune_5, getCellIdsRandomCellTypes(5, c(immune_cells, no_immune_5)))

no_b_10 <- 
  c(no_b_5, getCellIdsRandomCellTypes(5, c(b_cells, no_b_5)))

no_leuko_10 <- 
  c(no_leuko_5, getCellIdsRandomCellTypes(5, c(leuko_cells, no_leuko_5)))


b_no_b_5_cells  <- union(b_cells, no_b_5)
b_no_b_10_cells <- union(b_cells, no_b_10)

b_no_immune_5_cells <- union(b_cells, no_immune_5)
b_no_immune_10_cells <- union(b_cells, no_immune_10)

leuko_no_immune_5_cells <- union(leuko_cells, no_immune_5)
leuko_no_immune_10_cells <- union(leuko_cells, no_immune_10)

leuko_b_no_immune_5_cells <- union(leuko_cells, b_no_immune_5_cells)

b_endo_no_immune_3_cells <- union(b_cells,
                                  union(endo_cells, no_immune_no_endo_3))
leuko_endo_no_immune_3_cells <- union(leuko_cells,
                                      union(endo_cells, no_immune_no_endo_3))

## TODO: SPECTRAL DIMS option as constant at the top
## TODO: removing unused computations.
## TODO: re-run from the top.

### Distance compute pipeline for One Tissue.
ComputeDistCellVec <- function(cell_vec) {
  return(
    TabMurisDistResPipeline(data_full[,cell_vec], 
                            meta_data_full, 
                            filename = NULL, 
                            spectral_dims = c(3, 5, 10), 
                            n_clust = N_CLUST,
                            n_cores = N_CORES,
                            R = R,
                            G = 300,
                            PCs = 25,
                            var_features = VAR_FEATURES)
  )
}


if(RERUN_CELLTYPE_DIST) {
  
  immune_cells_full_dist_res <- 
    ComputeDistCellVec(immune_cells)
  
  leuko_b_no_immune_5_cells_full_dist_res <- 
    ComputeDistCellVec(leuko_b_no_immune_5_cells)
  
  b_no_immune_5_cells_full_dist_res <- 
    ComputeDistCellVec(b_no_immune_5_cells)
  leuko_no_immune_5_cells_full_dist_res <-
    ComputeDistCellVec(leuko_no_immune_5_cells)
  
  leuko_endo_no_immune_3_cells_full_dist_res <- 
    ComputeDistCellVec(leuko_endo_no_immune_3_cells)
  
  b_endo_no_immune_3_cells_full_dist_res <- 
    ComputeDistCellVec(b_endo_no_immune_3_cells)
  
  ## 10 cells and combined 5:
  leuko_b_no_immune_5_cells_full_dist_res <- 
    ComputeDistCellVec(leuko_b_no_immune_5_cells)
  b_no_immune_10_cells_full_dist_res <- 
    ComputeDistCellVec(b_no_immune_10_cells)
  leuko_no_immune_10_cells_full_dist_res <-
    ComputeDistCellVec(leuko_no_immune_10_cells)
  
  cell_type_dist_res_list <- 
    list(
      immune_cells_full_dist_res = immune_cells_full_dist_res,
      
      b_cells_full_dist_res = b_cells_full_dist_res,
      b_no_immune_5_cells_full_dist_res = b_no_immune_5_cells_full_dist_res,
      leuko_no_immune_5_cells_full_dist_res = leuko_no_immune_5_cells_full_dist_res,
      leuko_b_no_immune_5_cells_full_dist_res = leuko_b_no_immune_5_cells_full_dist_res,
      b_no_immune_10_cells_full_dist_res = b_no_immune_10_cells_full_dist_res,
      leuko_no_immune_10_cells_full_dist_res = leuko_no_immune_10_cells_full_dist_res,
      
      
      leuko_endo_no_immune_3_cells_full_dist_res = leuko_endo_no_immune_3_cells_full_dist_res,
      b_endo_no_immune_3_cells_full_dist_res = b_endo_no_immune_3_cells_full_dist_res
      
    )
  
  save(cell_type_dist_res_list, 
       file = here("output/tab_muris_sc/dist_res/cell_type_dist_res.RData"))
  
}

load(here("output/tab_muris_sc/dist_res/cell_type_dist_res.RData"), 
     verbose = TRUE)

############################### ===============================================
### Label Ordering for Heatmaps ===============================================
############################### ===============================================


### Cell Type Ordering for Heatmaps:
leuko_endo_no_immune_3_levels <- 
  c(## No Immune / Endo 3:
    "Bladder Urothelial Cell, Bladder",
    "Skeletal Muscle Satellite Cell, Limb Muscle",
    "Mesenchymal Cell, Trachea",
    ## Endo Cells:
    "Endocardial Cell, Heart and Aorta",
    "Endothelial Cell of Hepatic Sinusoid, Liver",
    "Kidney Capillary Endothelial Cell, Kidney",
    "Lung Endothelial Cell, Lung", 
    "Endothelial Cell, Bladder",
    "Endothelial Cell, Heart and Aorta",
    "Endothelial Cell, Limb Muscle",
    "Endothelial Cell, Mammary Gland",
    "Endothelial Cell, Trachea", 
    ## Leukocytes:
    "Leukocyte, Liver", 
    "Leukocyte, Bladder",
    "Leukocyte, Thymus",
    "Leukocyte, Kidney",
    "Leukocyte, Lung"
  )

b_endo_no_immune_3_levels <- 
  c(## No Immune / Endo 3:
    "Bladder Urothelial Cell, Bladder",
    "Skeletal Muscel Satellite Cell, Limb Muscle",
    "Mesenchymal Cell, Trachea",
    ## Endo Cells:
    "Endocardial Cell, Heart and Aorta",
    "Endothelial Cell of Hepatic Sinusoid, Liver",
    "Kidney Capillary Endothelial Cell, Kidney",
    "Lung Endothelial Cell, Lung", 
    "Endothelial Cell, Bladder",
    "Endothelial Cell, Heart and Aorta",
    "Endothelial Cell, Limb Muscle",
    "Endothelial Cell, Mammary Gland",
    "Endothelial Cell, Trachea", 
    ## B Cells:
    "B Cell, Spleen", 
    "B Cell, Lung", 
    "B Cell, Limb Muscle", 
    "B Cell, Mammary Gland"
  )


b_no_immune_5_levels <- c("Type II Pneumocyte, Lung", 
                          "Basal Cell of Epidermis, Tongue",
                          "Endocardial Cell, Heart and Aorta",
                          "Endothelial Cell, Limb Muscle", 
                          "Lung Endothelial Cell, Lung", 
                          "B Cell, Spleen", 
                          "B Cell, Lung", 
                          "B Cell, Limb Muscle", 
                          "B Cell, Mammary Gland")

leuko_no_immune_5_levels <- c(
  "Type II Pneumocyte, Lung",
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

leuko_b_no_immune_5_levels <- 
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
  function(cell_vec, dist_mat, meta_data,
           title = "",
           include_tissue = FALSE,
           hclust = FALSE,
           preset_levels = NULL,
           numeric_x=TRUE,
           numeric_y=FALSE) {
    
    rownames(meta_data) = meta_data$cell_id
    group_labels <- 
      meta_data[cell_vec, ]$cell_type
    
    if(include_tissue) {
      group_labels <- group_labels %p% ", " %p% meta_data[cell_vec, ]$tissue
    }
    
    ## cleaning of labels.
    group_labels <- str_to_title(gsub("_", " ", group_labels))
    ## title case corrections
    group_labels <- gsub(" And", " and", group_labels)
    group_labels <- gsub(" Of", " of", group_labels)
    group_labels <- gsub(" Ii", " II", group_labels)
    
    group_labels <- as.factor(group_labels)
    
    groupwiseDistanceHeatmap(group_labels, 
                             dist_mat,
                             title = title,
                             hclust = hclust,
                             preset_levels = preset_levels,
                             numeric_x = numeric_x,
                             numeric_y = numeric_y)  
  }




createHeatmapGrid <- 
  function(cell_vec, dist_mat, full_dist_mat,
           meta_data, type_levels,
           fig1_w_ratio = 1.5,
           titles = c("(a)", "(b)", "(c)"),
           main_title = NULL,
           main_title_size = 8,
           do_legend = TRUE) {
  
  ## 1. distance and heatmap computation restricted to cell types of interest.
  fig_1 <-
    groupwiseDistanceHeatmapCellType(
      cell_vec = cell_vec, 
      dist_mat = dist_mat,
      meta_data = meta_data_full,
      include_tissue = TRUE,
      #hclust = TRUE,
      preset_levels = type_levels,
      numeric_x = TRUE,
      title = titles[1]) + 
    theme_bw()
  
  ## 2. truncating the full distance matrix
  fig_2 <-
    groupwiseDistanceHeatmapCellType(
      cell_vec = cell_vec, 
      dist_mat = selectFullDistIdx(cell_vec,
                                   dist = full_dist_mat,
                                   meta_data = meta_data),
      meta_data = meta_data,
      include_tissue = TRUE,
      #hclust = TRUE,
      preset_levels = type_levels,
      numeric_x = TRUE,
      numeric_y = TRUE,
      title = titles[2]) + 
    theme_bw() +
    theme(axis.title.y = element_blank())
  
  ## 3. truncating the full heatmap
  fig_3 <-
    groupwiseDistanceHeatmapCellType(
      cell_vec = meta_data$cell_id, 
      dist_mat = full_dist_mat,
      meta_data = meta_data,
      include_tissue = TRUE,
      #hclust = TRUE,
      preset_levels = type_levels,
      numeric_x = TRUE,
      numeric_y = TRUE,
      title = titles[3]) + 
    theme_bw() + 
    theme(axis.title.y = element_blank())
  
  full_fig <-
    ggarrange(fig_1, fig_2, fig_3, ncol=3, 
              common.legend = do_legend,
              legend = ifelse(do_legend, "bottom", "none"),
              widths = c(fig1_w_ratio, 1, 1))
  
  if(!is.null(main_title)) {
    full_fig <- 
      annotate_figure(full_fig, top = text_grob(
      main_title, face = "bold", size = main_title_size))
  }
  
  full_fig
  return(full_fig)
  
}


### Load results from full tissue distance compute.
load(file = getCellTypeSampleFileString("droplet_only"), verbose = TRUE) 

rownames(meta_data_full) <- meta_data_full$cell_id


###############################################################################
### Fig 1: Leukocytes, 5 Other Cell Types
###############################################################################

FIG1_RATIO = 2
FIG1_TITLE_SIZE = 10

fig1_side <-
  createHeatmapGrid(
    cell_vec = leuko_no_immune_5_cells, 
    dist_mat = leuko_no_immune_5_cells_full_dist_res$dist_list$
      side_ref_g300_dist, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_no_immune_5_levels,
    fig1_w_ratio = FIG1_RATIO,
    main_title = "SIDEREF",
    main_title_size = FIG1_TITLE_SIZE,
    do_legend = FALSE)

fig1_spect3 <-
  createHeatmapGrid(
    cell_vec = leuko_no_immune_5_cells, 
    dist_mat = leuko_no_immune_5_cells_full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_no_immune_5_levels,
    fig1_w_ratio = FIG1_RATIO,
    title = c("(d)", "(e)", "(f)"),
    main_title = "SIDEREF Spectral (3 Dims.)",
    main_title_size = FIG1_TITLE_SIZE,
    do_legend = FALSE)

fig1_spect5 <-
  createHeatmapGrid(
    cell_vec = leuko_no_immune_5_cells, 
    dist_mat = leuko_no_immune_5_cells_full_dist_res$dist_list$side_ref_g300_dist_spectral_d5_dist, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d5_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_no_immune_5_levels,
    fig1_w_ratio = FIG1_RATIO,
    titles = c("(g)", "(h)", "(i)"),
    main_title = "SIDEREF Spectral (5 Dims.)",
    main_title_size = FIG1_TITLE_SIZE)

ggarrange(fig1_side, fig1_spect3, fig1_spect5, nrow=3,
          heights = c(1, 1, 1.2))

fig1_pca <-
  createHeatmapGrid(
    cell_vec = leuko_no_immune_5_cells, 
    dist_mat = leuko_no_immune_5_cells_full_dist_res$dist_list$
      pca_wtd25_dist, 
    full_dist_mat = full_dist_res$dist_list$pca_wtd25_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_no_immune_5_levels,
    fig1_w_ratio = FIG1_RATIO,
    main_title = "PCA (25 Dims.)",
    main_title_size = FIG1_TITLE_SIZE,
    do_legend = FALSE)

fig1_pear <-
  createHeatmapGrid(
    cell_vec = leuko_no_immune_5_cells, 
    dist_mat = leuko_no_immune_5_cells_full_dist_res$dist_list$pear_dist, 
    full_dist_mat = full_dist_res$dist_list$pear_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_no_immune_5_levels,
    fig1_w_ratio = FIG1_RATIO,
    title = c("(d)", "(e)", "(f)"),
    main_title = "1 - |Pearson Correlation|",
    main_title_size = FIG1_TITLE_SIZE,
    do_legend = FALSE)

fig1_spearman <-
  createHeatmapGrid(
    cell_vec = leuko_no_immune_5_cells, 
    dist_mat = leuko_no_immune_5_cells_full_dist_res$dist_list$spearman_dist, 
    full_dist_mat = full_dist_res$dist_list$spearman_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_no_immune_5_levels,
    fig1_w_ratio = FIG1_RATIO,
    titles = c("(g)", "(h)", "(i)"),
    main_title = "1 - |Spearman Correlation|",
    main_title_size = FIG1_TITLE_SIZE)

ggarrange(fig1_pca, fig1_pear, fig1_spearman, nrow=3,
          heights = c(1, 1, 1.2))

###############################################################################
### Set 3: B+ Leukocytes, 5 Other Cell Types
###############################################################################


FIG2_2_RATIO = 1.75
FIG2_TITLE_SIZE = 10

fig2_sideref <-
  createHeatmapGrid(
    cell_vec = leuko_b_no_immune_5_cells, 
    dist_mat = leuko_b_no_immune_5_cells_full_dist_res$dist_list$
      side_ref_g300_dist, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_b_no_immune_5_levels,
    fig1_w_ratio = FIG2_2_RATIO,
    titles = c("(a)", "(b)", "(c)"),
    main_title = "SIDEREF",
    main_title_size = FIG2_TITLE_SIZE,
    do_legend = FALSE)

fig2_sideref_spect3 <-
  createHeatmapGrid(
    cell_vec = leuko_b_no_immune_5_cells, 
    dist_mat = leuko_b_no_immune_5_cells_full_dist_res$dist_list$
      side_ref_g300_dist_spectral_d3_dist, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d3_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_b_no_immune_5_levels,
    fig1_w_ratio = FIG2_2_RATIO,
    titles = c("(d)", "(e)", "(f)"),
    main_title = "SIDEREF Spectral (3 Dims.)",
    main_title_size = FIG2_TITLE_SIZE)


fig2_sideref_spect5 <-
  createHeatmapGrid(
    cell_vec = leuko_b_no_immune_5_cells, 
    dist_mat = leuko_b_no_immune_5_cells_full_dist_res$dist_list$
      side_ref_g300_dist_spectral_d5_dist, 
    full_dist_mat = full_dist_res$dist_list$side_ref_g300_dist_spectral_d5_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_b_no_immune_5_levels,
    fig1_w_ratio = FIG2_2_RATIO,
    titles = c("(d)", "(e)", "(f)"),
    main_title = "SIDEREF Spectral (5 Dims.)",
    main_title_size = FIG2_TITLE_SIZE)

ggarrange(fig2_sideref, fig2_sideref_spect3, nrow=2,
          heights = c(1, 1.25))



fig2_pca <-
  createHeatmapGrid(
    cell_vec = leuko_b_no_immune_5_cells, 
    dist_mat = leuko_b_no_immune_5_cells_full_dist_res$dist_list$
      pca_wtd25_dist, 
    full_dist_mat = full_dist_res$dist_list$pca_wtd25_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_b_no_immune_5_levels,
    fig1_w_ratio = FIG2_2_RATIO,
    titles = c("(d)", "(e)", "(f)"),
    main_title = "PCA (25 Dims.)",
    main_title_size = FIG2_TITLE_SIZE,
    do_legend = FALSE)

fig2_spearman <-
  createHeatmapGrid(
    cell_vec = leuko_b_no_immune_5_cells, 
    dist_mat = leuko_b_no_immune_5_cells_full_dist_res$dist_list$
      spearman_dist, 
    full_dist_mat = full_dist_res$dist_list$spearman_dist,
    meta_data = meta_data_full, 
    type_levels = leuko_b_no_immune_5_levels,
    fig1_w_ratio = FIG2_2_RATIO,
    titles = c("(g)", "(h)", "(i)"),
    main_title = "1 - |Spearman Correlation|",
    main_title_size = FIG2_TITLE_SIZE)

ggarrange(fig2_pca, fig2_spearman, nrow=2,
          heights = c(1, 1.25))




###############################################################################
### Gene Set Enrichment for B-Cells within Kidney, Thymus Leukocytes.
###############################################################################

tm_preproc <-
  TabMurisProcSubset(full_data = data_full, 
                     meta_data = meta_data_full, 
                     var_features= VAR_FEATURES,
                     return_seurat = TRUE)[[1]]



## DE Genes among the B-cell groups
getMarkerGSList <- function(assay, cells_in, cells_vs, 
                            cell_marker_list_sizes = c(10, 50, 100),
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
                          assay, meta_data) {
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
                  logfc.threshold = 0)
    
    de_marker_list[[i]] <- 
      prepRanks(de_marker_list[[i]],
                all_gene_names = 
                  marker_gs_list[[length(marker_gs_list)]])
    
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
    labs(y = "Abs. Norm. Enrichment Score", x = "B Cell D.E. Gene Set", 
         fill = "Cell Type, Organ")
  
  return(p)
  
}

p_gsea_bar <-
  BarplotGSEA(gsea_res_vs_all_other_cells, 
              preset_levels = c("Top 10", "Top 50", "Top 100"))


p_gsea_bar

df %>% dplyr::select(pathway, pval, celltype) 
