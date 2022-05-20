library(here)
source(here("R/libraries.R"))
source(here('R/relative_group_dist_comps.R'))
source(here('R/perf_metrics_funcs.R'))
source(here('R/distance_compute_helpers.R'))
source(here("R/tm_helper_functions.R"))

set.seed(1)


##### DATA LOAD, CLEAN, SELECT CLEAR GLOBAL GROUP CELL TYPE VECS.
FILE_TAGS <- 
  c("main_dist_comps", "othr_scrna_seq_dist_comps")


for(f in FILE_TAGS) {
  load(
    getTMCellTypeSampleDistFileName(
      file_tag = f,
      file_loc = here("output/tab_muris_sc/dist_res/"),
      samp_size = TM_SAMP_SIZE,
      use_droplet_only = USE_DROPLET_ONLY),
    verbose = TRUE) 
}

dist_list_to_measure <- 
  list(
    side_ref_g300_dist = main_dist_output$dist_list$side_ref_g300_dist,
    pca_wtd25_dist = main_dist_output$dist_list$pca_wtd25_dist,
    euclid_dist = main_dist_output$dist_list$euclid_dist,
    pear_dist = main_dist_output$dist_list$pear_dist,
    spearman_dist = main_dist_output$dist_list$spearman_dist,
    RAFSIL = othr_dist_output$dist_list$RAFSIL,
    SIMLR_25_dims = othr_dist_output$dist_list$SIMLR_25_dims
  )

load(
  getTMCellTypeSampleDistFileName(
    file_tag = "meta_data",
    file_loc = here("output/tab_muris_sc/dist_res/"),
    samp_size = TM_SAMP_SIZE,
    use_droplet_only = USE_DROPLET_ONLY),
  verbose = TRUE)

meta_data_full <- 
  prepTMGroupLabels(meta_data = meta_data_full, include_tissue = TRUE) 

###############################################################################
### Enumerate Predefined Global Groups.
###############################################################################


### Enumerating the cell type annotations from the data into relevant groupings.
broad_cell_types <- 
  c("leukocyte", "blood cell", "hematopoietic precursor cell", "myeloid cell")


granulocytes <- c("granulocyte", "granulocytopoietic cell", "basophil", "mast cell") 
monocytes <- c("monocyte", "promonocyte", "classical monocyte", "non-classical monocyte",
               "dendritic cell", "alveolar macrophage", "macrophage", 
               "Langerhans cell")
rbcs <- c("erythroblast", "proerythroblast")

nk_cells <- c("natural killer cell")
t_cells <- c("T cell")
b_cells <- c("B cell")
pro_b_cells <- c("early pro-B cell", "Fraction A pre-pro B cell", "late pro-B cell")
stromal <- c("stromal cell", "mesenchymal cell", "mesenchymal stem cell")

endothelial <- c("endothelial cell", "endothelial cell of hepatic sinusoid", 
                 "kidney capillary endothelial cell", "lung endothelial cell")
    

## Other:
basal <- c("basal cell", "basal cell of epidermis")
    
epithelial <- c("epithelial cell",
                "duct epithelial cell", "kidney collecting duct epithelial cell", 
                "luminal epithelial cell of mammary gland", "type II pneumocyte",
                "kidney proximal straight tubule epithelial cell", 
                "kidney loop of Henle ascending limb epithelial cell")


#### SINGLETONS:
others <-
  c("bladder cell", "bladder urothelial cell",
    "cardiac muscle cell",
    "ciliated columnar cell of tracheobronchial tree",
    "endocardial cell",
    "fibroblast", 
    "hepatocyte", 
    "keratinocyte",
    "kidney cell", "mesangial cell",
    "neuroendocrine cell", "skeletal muscle satellite cell")

other_immunes <- c("DN1 thymic pro-T cell")

#### Progenitor, i.e., Cells that will develop into more distinguished types.
progenitor_cells <-
  c("promonocyte",
    "proerythroblast",
    "immature T cell","immature B cell", 
    "hematopoietic precursor cell",
    "stromal cell", "mesenchymal cell", "mesenchymal stem cell")

###############################################################################
### Build global group data set.
###############################################################################

global_group_list <- 
  list(singletons = others,
       broad_cell_types = broad_cell_types,
       granulocytes = granulocytes,
       monocytes = monocytes,
       nk_cells = nk_cells,
       t_cells = t_cells,
       b_cells = b_cells,
       pro_b_cells = pro_b_cells,
       stromal = stromal, ### Note These are progenitors as well
       endothelial = endothelial,
       basal = basal,
       epithelial = epithelial,
       rbcs = rbcs,
       other_immunes = other_immunes)

group_global_group_map <- 
  lapply(seq_len(length(global_group_list)), 
         function(i) return(
           data.frame(cell_group_low = global_group_list[[i]],
                      cell_group_high = names(global_group_list)[i]
                      )
         )) %>% purrr::reduce(rbind) %>% 
  mutate(is_singleton = ifelse(cell_group_high == "singletons", TRUE, FALSE),
         cell_group_high = 
           ifelse(cell_group_high == "singletons", cell_group_low, cell_group_high)) %>% 
  left_join(meta_data_full %>%
              dplyr::distinct(cell_type, tissue) %>% 
              dplyr::select(cell_group_low = cell_type, tissue))



group_global_group_map_all <- 
  group_global_group_map %>% 
  dplyr::filter(
    ## removing other immune cells and broad cells.
    !cell_group_high %in% c("broad_cell_types", 
                            "granulocytes", 
                            "monocytes",
                            "rbcs",
                            "other_immunes"),
    ## removing progenitor cells.
    !cell_group_low %in% c(progenitor_cells)
  ) %>%
  mutate(cell_group_high = 
           ifelse(cell_group_high %in% c(
             "nk_cells",
             "t_cells",
             "b_cells",
             "pro_b_cells",
             "endothelial",
             "basal"), cell_group_high, cell_group_low)) %>% 
  mutate(cell_group_low = cell_group_low %p% ", " %p% tissue) %>%
  dplyr::select(cell_group_low, cell_group_high)


### Only the global cell groups.
group_global_group_map_only_globals <- 
  group_global_group_map %>% 
  dplyr::filter(cell_group_high %in% c(
    "nk_cells",
    "t_cells",
    "b_cells",
    "pro_b_cells",
    "endothelial",
    "basal"
  )) %>%
  mutate(cell_group_low = cell_group_low %p% ", " %p% tissue) %>%
  dplyr::select(cell_group_low, cell_group_high)

## For multi:
group_global_group_map_all_multi <- 
  group_global_group_map_all %>% 
  dplyr::rename(cell_group_high1 = cell_group_high) %>%
  mutate(cell_group_high2 = 
           case_when(
             cell_group_high1 %in% c("nk_cells", "t_cells", "b_cells", "pro_b_cells") ~ "Immune",
             TRUE ~ cell_group_high1
           ))

## For multi:
group_global_group_map_only_globals_multi <- 
  group_global_group_map_only_globals %>% 
  dplyr::rename(cell_group_high1 = cell_group_high) %>%
  mutate(cell_group_high2 = 
           case_when(
             cell_group_high1 %in% c("nk_cells", "t_cells", "b_cells", "pro_b_cells") ~ "Immune",
             TRUE ~ cell_group_high1
           ))



computeGrpwiseStatsTM <- 
  function(meta_data,
           global_group_map_multi,
           dist_list) {
    
    ### prepare group labels, filter:
    meta_data_subset <- 
      meta_data %>% 
      mutate(cell_group_low = cell_type %p% ", " %p% tissue)
    
    meta_idx <- which(meta_data_subset$cell_group_low %in% 
                        global_group_map_multi$cell_group_low)
    
    meta_data_subset <- 
      meta_data_subset %>%
      dplyr::filter(cell_group_low %in% global_group_map_multi$cell_group_low)
    
    
    group_labels <- 
      meta_data_subset$cell_group_low
    
    
    ### Non multi-hierarchy evaluations:
    grpwise_max_vio <- sapply(
      dist_list, 
      function(dist) {
        dist <- dist[meta_idx, meta_idx]
        countGlobalGrpMaxViolationsMulti(
          dist, 
          group_labels = group_labels,
          group_global_group_map = global_group_map_multi %>% 
            dplyr::select(cell_group_low, cell_group_high1))}
    )
    
    grpwise_min_vio <- sapply(
      dist_list,
      function(dist) {
        dist <- dist[meta_idx, meta_idx]
        countGlobalGrpMinViolationsMulti(
          dist, 
          group_labels = group_labels,
          group_global_group_map = global_group_map_multi %>% 
            dplyr::select(cell_group_low, cell_group_high1))}
    )
    
    ### Multi-hierarchy evaluations:
    grpwise_max_multi_vio <- sapply(
      dist_list,
      function(dist) {
        dist <- dist[meta_idx, meta_idx]
        countGlobalGrpMaxViolationsMulti(
          dist, 
          group_labels = group_labels,
          group_global_group_map = global_group_map_multi)}
    )
    
    grpwise_min_multi_vio <- sapply(
      dist_list,
      function(dist) {
        dist <- dist[meta_idx, meta_idx]
        countGlobalGrpMinViolationsMulti(
          dist, 
          group_labels = group_labels,
          group_global_group_map = global_group_map_multi)}
    )
    
    return(list(
      grpwise_max_vio = grpwise_max_vio,
      grpwise_min_vio = grpwise_min_vio,
      grpwise_max_multi_vio = grpwise_max_multi_vio,
      grpwise_min_multi_vio = grpwise_min_multi_vio
    ))
}


res_list_only_globals <- 
  computeGrpwiseStatsTM(meta_data = meta_data_full,
                        global_group_map_multi = group_global_group_map_only_globals_multi,
                        dist_list = dist_list_to_measure)


res_list_all <- 
  computeGrpwiseStatsTM(meta_data = meta_data_full,
                        global_group_map_multi = group_global_group_map_all_multi,
                        dist_list = dist_list_to_measure)

lapply(res_list_only_globals, function(x)x * 100)
lapply(res_list_all, function(x)x * 100)

