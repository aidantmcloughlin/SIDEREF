library(here)
source(here("R/libraries.R"))
source(here('R/relative_group_dist_comps.R'))
source(here('R/distance_compute_helpers.R'))
source(here("R/tm_helper_functions.R"))
source(here("R/bipartite_graphs.R"))

set.seed(1)


############ PLOTTING CONSTANTS ###############
BIP_FIG_TITLE_SIZE <- 8
HM_FIG_TITLE_SIZE <- 8
TWO_COLOR_GRAD <- FALSE

FILE_TAGS <- 
  c("main_dist_comps", "othr_scrna_seq_dist_comps", "spectral_dist_comps")

## Get PCs:
load(
  getTMCellTypeSampleDistFileName(
    file_tag = FILE_TAGS[1],
    file_loc = here("output/tab_muris_sc/dist_res/"),
    samp_size = TM_SAMP_SIZE,
    use_droplet_only = USE_DROPLET_ONLY),
  verbose = TRUE)

dist_list_names <- names(main_dist_output$dist_list)
pcs <- dist_list_names[which(grepl("pca_wtd", dist_list_names))]
pcs <- as.numeric(str_extract(pcs, '[:digit:]+'))

rm(main_dist_output)
gc()

## Get titles from dist list
getDistTitles <- function(
  pcs, spect_dims = TM_SPECTRAL_DIMS) {
  
  title_list_main <- c("SIDEREF", "PCA (" %p% pcs %p% " Dims.)",
                       "Euclidean", "1 - |Pearson Correlation|",
                       "1 - |Spearman Correlation|")
  
  
  title_list_other <- c("RAFSIL", "SIMLR (C = " %p% spect_dims %p% ")")
  
  title_list_spect <- c("SIDEREF Spectral (" %p% spect_dims %p% " Dims.)")
  for(i in pcs) {
    title_list_spect <- append(
      title_list_spect,
      c("PCA (" %p% i %p% " Dims.)" %p% " Spectral (" %p% 
          spect_dims %p% " Dims.)")
    )
  }
  
  return(list(title_list_main = title_list_main,
              title_list_other = title_list_other,
              title_list_spect = title_list_spect))
}

plot_titles_list <-
  getDistTitles(pcs = pcs, spect_dims = TM_SPECTRAL_DIMS)

############################### ===============================================
### Factor Label Orders for Heatmaps ==========================================
############################### ===============================================


cells_leuko_levels <- 
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

cells_b_levels <- 
  c("Type II Pneumocyte, Lung", 
    "Basal Cell of Epidermis, Tongue",
    "Endocardial Cell, Heart and Aorta",
    "Endothelial Cell, Limb Muscle", 
    "Lung Endothelial Cell, Lung", 
    "B Cell, Spleen", 
    "B Cell, Lung", 
    "B Cell, Limb Muscle", 
    "B Cell, Mammary Gland")

cells_b_leuko_levels <- 
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

cell_type_levels_list <- 
  list(cells_leuko = cells_leuko_levels,
       cells_b = cells_b_levels,
       cells_b_leuko = cells_b_leuko_levels
  )



#################################### ==========================================
### Selecting the IDs of these cells types from given meta data set
#################################### ==========================================

### Displaying all 80 cell types on the heatmap is unwieldy. In practice, one 
##   has a subset of cell groups for which the global comparison will be useful, 
##   such as broadly defined cell types compared to more narrowly defined cell 
##   types.

GetTMHeatmapsCellIDs <- function(meta_data_full) {
  ## B Cells and Leukocytes
  b_cells     <- meta_data_full$cell_id[
    which(meta_data_full$cell_type == "B cell")]
  leuko_cells <- meta_data_full$cell_id[
    which(meta_data_full$cell_type == "leukocyte")]
  
  ## Endothelial Cells
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
  
  
  cells_leuko <- union(leuko_cells, 
                       union(endo_selected_cells,
                             unrelated_cells))
  
  cells_b <- union(b_cells, 
                   union(endo_selected_cells,
                         unrelated_cells))
  
  cells_b_leuko <- union(cells_leuko,
                         cells_b)
  
  return(list(cells_leuko = cells_leuko,
              cells_b = cells_b,
              cells_b_leuko = cells_b_leuko))
}



##################################################### =========================
###### Plotting Functions =====================================================
##################################################### =========================


##################  bipartite network ##################  
bipartiteTMCells <- 
  function(cell_id_vec,
           dist, 
           meta_data,
           preset_levels,
           source_groups = c("Leukocyte, Bladder", "Leukocyte, Kidney"),
           title = "",
           two_color_grad = TWO_COLOR_GRAD) {
    rownames(meta_data) = meta_data$cell_id
    dist_mat <- selectFullDistIdx(cell_id_vec,
                                  dist = dist,
                                  meta_data = meta_data)
    
    
    p <-
      bipartiteNetworkGraph(
        dist_mat,
        group_labels = as.character(meta_data[cell_id_vec, ]$group_label),
        preset_levels = preset_levels,
        source_groups = source_groups,
        node_size = 7.5,
        text_size = 3,
        s_t_size = 3,
        title = title,
        two_color_grad = two_color_grad,
        legend_title = "Cell Type, Tissue") +
      theme(plot.title = element_text(size=BIP_FIG_TITLE_SIZE))
    
    return(p)
}



##################  groupwise distance heatmap ##################  
groupwiseDistanceHeatmapCellType <- 
  function(cell_id_vec, dist_mat, meta_data,
           preset_levels,
           title = "",
           symmetrize = TRUE,
           do_numeric_y = FALSE) {
    
    rownames(meta_data) = meta_data$cell_id
    
    dist_mat <- selectFullDistIdx(cell_id_vec,
                                  dist = dist_mat,
                                  meta_data = meta_data)
    
    group_labels <- getGroupLabels(meta_data, cell_id_vec = cell_id_vec)
    
    p <-
      groupwiseDistanceHeatmap(group_labels, 
                               dist_mat,
                               symmetrize = symmetrize,
                               title = title,
                               do_hclust_axes = FALSE,
                               preset_levels = preset_levels,
                               do_numeric_x=TRUE,
                               do_numeric_y=do_numeric_y)  +
      theme(plot.title = element_text(size=HM_FIG_TITLE_SIZE))
    
    return(p)
  }



###############################################################################
###### MAIN    ================================================================
###############################################################################

### Load results from full tissue distance compute.
processTMHeatmaps <- function(
  i,
  file_tags,
  plot_titles_list,
  cell_type_levels_list) {
  load(
    getTMCellTypeSampleDistFileName(
      file_tag = file_tags[i],
      file_loc = here("output/tab_muris_sc/dist_res/"),
      samp_size = TM_SAMP_SIZE,
      use_droplet_only = USE_DROPLET_ONLY),
    verbose = TRUE)
  
  load(
    getTMCellTypeSampleDistFileName(
      file_tag = "meta_data",
      file_loc = here("output/tab_muris_sc/dist_res/"),
      samp_size = TM_SAMP_SIZE,
      use_droplet_only = USE_DROPLET_ONLY),
    verbose = TRUE)
  
  meta_data_full <- 
    prepTMGroupLabels(meta_data = meta_data_full, include_tissue = TRUE) 
  
  hm_cell_ids <- GetTMHeatmapsCellIDs(meta_data_full)
  stopifnot(all.equal(names(cell_type_levels_list),
                      names(hm_cell_ids)))
  
  
  ### MAIN PLOT LOOP =================
  if(i==1) {
    ### FILTER
    non_unwtd_pca_idx <- 
      which(!grepl("pca_[[:digit:]]", names(main_dist_output$dist_list)))
    main_dist_output[[1]] <- 
      main_dist_output[[1]][non_unwtd_pca_idx]
    main_dist_output[[2]] <- 
      main_dist_output[[2]][non_unwtd_pca_idx]
    
    dist_obj_name <- "main_dist_output"
  } else if(i == 2) {
    
    dist_obj_name <- "othr_dist_output"
    
  } else if(i == 3) {
    ### FILTER:
    spectral_idx <-
      which(grepl("spect", names(spectral_dist_res$dist_list)))
    spectral_dist_res[[1]] <- 
      spectral_dist_res[[1]][spectral_idx]
    spectral_dist_res[[2]] <- 
      spectral_dist_res[[2]][spectral_idx]
    
    dist_obj_name <- "spectral_dist_res" 
  }
  
  n_cell_groups <- length(hm_cell_ids)
  n_dists <- length(get(dist_obj_name)[['dist_list']])
  
  plot_groups <- vector(mode = "list", length = n_cell_groups)
  names(plot_groups) <- names(hm_cell_ids)
  
  for(n in seq_len(n_cell_groups)) {
    ### PLOT:
    plot_list_bip     <- vector(mode = "list", length = n_dists)
    plot_list_hm_asym <- vector(mode = "list", length = n_dists)
    plot_list_hm_sym  <- vector(mode = "list", length = n_dists)
    
    
    names(plot_list_bip) <- names(get(dist_obj_name)[['dist_list']])
    names(plot_list_hm_asym) <- names(get(dist_obj_name)[['dist_list']])
    names(plot_list_hm_sym) <- names(get(dist_obj_name)[['dist_list']])
    
    for(p in seq_len(n_dists)) {
      if(names(hm_cell_ids)[n] == "cells_b_leuko") {
        plot_list_bip[[p]] <- 
          bipartiteTMCells(
            cell_id_vec = hm_cell_ids[[n]],
            dist = get(dist_obj_name)[['dist_list']][[p]], 
            meta_data = meta_data_full,
            preset_levels = cell_type_levels_list[[n]],
            source_groups = c("Leukocyte, Bladder", "Leukocyte, Kidney"),
            title = plot_titles_list[[i]][p])
           
      }
      
      
      plot_list_hm_asym[[p]] <- 
        groupwiseDistanceHeatmapCellType(
          cell_id_vec = hm_cell_ids[[n]], 
          dist_mat = get(dist_obj_name)[['dist_list']][[p]], 
          meta_data = meta_data_full,
          preset_levels = cell_type_levels_list[[n]],
          title = plot_titles_list[[i]][p],
          do_numeric_y = 
            ifelse(plot_titles_list[[i]][p] %in% c("SIDEREF", 
                                                   "PCA (25 Dims.)",
                                                   "1 - |Pearson Correlation|"),
                   FALSE, TRUE),
          symmetrize = FALSE)
      
      plot_list_hm_sym[[p]] <- 
        groupwiseDistanceHeatmapCellType(
          cell_id_vec = hm_cell_ids[[n]], 
          dist_mat = get(dist_obj_name)[['dist_list']][[p]], 
          meta_data = meta_data_full,
          preset_levels = cell_type_levels_list[[n]],
          title = plot_titles_list[[i]][p],
          do_numeric_y = 
            ifelse(plot_titles_list[[i]][p] %in% c("SIDEREF", 
                                                   "PCA (25 Dims.)",
                                                   "1 - |Pearson Correlation|"),
                   FALSE, TRUE),
          symmetrize = TRUE)
    
    }
    #plot_groups <- vector(mode = "list", length = n_cell_groups)
    #names(plot_groups) <- names(hm_cell_ids)
    plot_groups[[n]] <- 
      list(plot_list_bip = plot_list_bip,
           plot_list_hm_asym = plot_list_hm_asym,
           plot_list_hm_sym = plot_list_hm_sym)
  }
  
  return(plot_groups)
}



tm_group_hm_list_main <-
  processTMHeatmaps(i = 1,
                    file_tags = FILE_TAGS,
                    plot_titles_list = plot_titles_list,
                    cell_type_levels_list = cell_type_levels_list)


tm_group_hm_list_other <-
  processTMHeatmaps(i = 2,
                    file_tags = FILE_TAGS,
                    plot_titles_list = plot_titles_list,
                    cell_type_levels_list = cell_type_levels_list)


tm_group_hm_list_spectral <-
  processTMHeatmaps(i = 3,
                    file_tags = FILE_TAGS,
                    plot_titles_list = plot_titles_list,
                    cell_type_levels_list = cell_type_levels_list)




arranged_p_1 <-
  ggarrange(
    tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$side_ref_g300_dist + 
      labs(title = "(a) " %p% tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$side_ref_g300_dist$labels$title),
    tm_group_hm_list_spectral$cells_b_leuko$plot_list_hm_asym$side_ref_g300_dist_spectral_d3_dist +
      labs(title = "(b) " %p% tm_group_hm_list_spectral$cells_b_leuko$plot_list_hm_asym$side_ref_g300_dist_spectral_d3_dist$labels$title),
    tm_group_hm_list_spectral$cells_b_leuko$plot_list_hm_asym$side_ref_g300_dist_spectral_d10_dist +
      labs(title = "(c) " %p% tm_group_hm_list_spectral$cells_b_leuko$plot_list_hm_asym$side_ref_g300_dist_spectral_d10_dist$labels$title),
    tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$pca_wtd25_dist +
      labs(title = "(d) " %p% tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$pca_wtd25_dist$labels$title),
    tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$euclid_dist +
      labs(title = "(e) " %p% tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$euclid_dist$labels$title),
    tm_group_hm_list_other$cells_b_leuko$plot_list_hm_asym$RAFSIL +
      labs(title = "(f) " %p% tm_group_hm_list_other$cells_b_leuko$plot_list_hm_asym$RAFSIL$labels$title),
    #ggplot() + theme_minimal(),
    tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$pear_dist +
      labs(title = "(g) " %p% tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$
             pear_dist$labels$title),
    tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$spearman_dist +
      labs(title = "(h) " %p% tm_group_hm_list_main$cells_b_leuko$plot_list_hm_asym$
             spearman_dist$labels$title),
    tm_group_hm_list_other$cells_b_leuko$plot_list_hm_asym$SIMLR_10_dims +
      labs(title = "(i) " %p% tm_group_hm_list_other$cells_b_leuko$plot_list_hm_asym$
             SIMLR_10_dims$labels$title), 
    nrow=3,
    ncol=3,
    widths = c(1.35, 1, 1),
    common.legend = TRUE,
    legend = "bottom"
  )


### Bipartite and Symmetrized SIDEREF sample
arranged_p_2 <- 
  ggarrange(tm_group_hm_list_main$cells_b_leuko$plot_list_hm_sym$
              side_ref_g300_dist,
            tm_group_hm_list_main$cells_b_leuko$plot_list_bip$
              side_ref_g300_dist +
              labs(title = ""),
            nrow = 2)
            

ggsave(here("manuscript_files/Figure4.eps"),
       plot = arranged_p_1,
       width = 17, height = 15,
       device='eps', dpi=250)


ggsave(here("manuscript_files/FigureS7.png"),
       plot = arranged_p_2,
       width = 8, height = 9.75,
       device='png', dpi=300)


