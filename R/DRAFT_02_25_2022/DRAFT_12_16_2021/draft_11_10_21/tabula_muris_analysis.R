### Tabula Muris Analysis
#library(GeneFishing)


## Load datasets 

liver_df <- 
  readRDS(here("data/tab_muris_bulk/Liver.rds"))

heart_aorta_df <- 
  readRDS(here("data/tab_muris_bulk/Heart_and_Aorta.rds"))

### clean data by size counts
liver_mat <- 
  as.matrix(t(t(liver_df$counts) / liver_df$size_factor))

heart_aorta_mat <- 
  as.matrix(t(t(heart_aorta_df$counts) / heart_aorta_df$size_factor))



### Cell type metadata

## Heart Aorta
load(here("data/tab_muris_bulk/droplet_Heart_and_Aorta_seurat_tiss.Robj"))

tiss_new <- UpdateSeuratObject(tiss)

colnames(heart_aorta_mat) = gsub("\\-1", "", colnames(heart_aorta_mat)) 

heart_aorta_cell_names <- data.frame(barcode = colnames(heart_aorta_mat))

heart_aorta_cell_names <- 
  heart_aorta_cell_names %>% 
  left_join(
    data.frame(barcode = rownames(tiss_new@meta.data),
               cell_type = tiss_new@meta.data$cell_ontology_class)) %>% 
  mutate(cell_type = as.character(cell_type)) %>% 
  mutate(cell_type = ifelse(is.na(cell_type), "erythrocyte",
                            cell_type))


rm(tiss_new, tiss)


## Liver:
load(here("data/tab_muris_bulk/droplet_Liver_seurat_tiss.Robj"))

tiss_new <- UpdateSeuratObject(tiss)

colnames(liver_mat) = gsub("\\-1", "", colnames(liver_mat)) 

liver_cell_names <- data.frame(barcode = colnames(liver_mat))

liver_cell_names <- 
  liver_cell_names %>% 
  left_join(
    data.frame(barcode = rownames(tiss_new@meta.data),
               cell_type = tiss_new@meta.data$cell_ontology_class)) %>% 
  mutate(cell_type = as.character(cell_type)) %>% 
  mutate(cell_type = ifelse(is.na(cell_type), "unlabeled",
                            cell_type))


rm(tiss_new, tiss)


## GO datasets

## Liver genes
GO_cholest_metab <- 
  fread(here("data/tab_muris_bulk/" %p%
               "GO_term_cholesterol_metabolism.csv"))


## heart, aorta genes

GO_card_musc_contract <- 
  fread(here("data/tab_muris_bulk/" %p% 
               "GO_term_cardiac_muscle_contraction.csv"))

GO_antigen_pres <- 
  fread(here("data/tab_muris_bulk/" %p% 
               "GO_term_antigen_presentation.csv"))

GO_fibro_prolif <- 
  fread(here("data/tab_muris_bulk/" %p% 
               "GO_term_fibroblast_proliferation.csv"))





### GeneFishing Replication


## Get rownames of considered bait genes

## Liver
bait_genes_chol_metab <-
  GO_cholest_metab$Symbol[GO_cholest_metab$Symbol %in% rownames(liver_mat)]

## Heart Aorta
bait_genes_card_musc_contract <-
  GO_card_musc_contract$Symbol[GO_card_musc_contract$Symbol %in% rownames(heart_aorta_mat)]

bait_genes_antigen_pres <-
  GO_antigen_pres$Symbol[GO_antigen_pres$Symbol %in% rownames(heart_aorta_mat)]

bait_genes_fibro_prolif <-
  GO_fibro_prolif$Symbol[GO_fibro_prolif$Symbol %in% rownames(heart_aorta_mat)]


### Cell type metadata

## Heart Aorta
load(here("data/tab_muris_bulk/droplet_Heart_and_Aorta_seurat_tiss.Robj"))

tiss_new <- UpdateSeuratObject(tiss)

colnames(heart_aorta_mat) = gsub("\\-1", "", colnames(heart_aorta_mat)) 

heart_aorta_cell_names <- data.frame(barcode = colnames(heart_aorta_mat))

heart_aorta_cell_names <- 
  heart_aorta_cell_names %>% 
  left_join(
    data.frame(barcode = rownames(tiss_new@meta.data),
               cell_type = tiss_new@meta.data$cell_ontology_class)) %>% 
  mutate(cell_type = as.character(cell_type)) %>% 
  mutate(cell_type = ifelse(is.na(cell_type), "erythrocyte",
                            cell_type))


save(heart_aorta_cell_names, file=here("output/tab_muris/heart_aorta_cell_names.RData"))

rm(tiss_new, tiss)


## Liver:
load(here("data/tab_muris_bulk/droplet_Liver_seurat_tiss.Robj"))

tiss_new <- UpdateSeuratObject(tiss)

colnames(liver_mat) = gsub("\\-1", "", colnames(liver_mat)) 

liver_cell_names <- data.frame(barcode = colnames(liver_mat))

liver_cell_names <- 
  liver_cell_names %>% 
  left_join(
    data.frame(barcode = rownames(tiss_new@meta.data),
               cell_type = tiss_new@meta.data$cell_ontology_class)) %>% 
  mutate(cell_type = as.character(cell_type)) %>% 
  mutate(cell_type = ifelse(is.na(cell_type), "unlabeled",
                            cell_type))


save(liver_cell_names, file=here("output/tab_muris/liver_cell_names.RData"))


rm(tiss_new, tiss)






## Load SIDEREF Results and Plot Heatmap
source(here("R/subgroup_composition_gene_contrib.R"))
load(here("output/tab_muris/heart_aorta_sideref.RData"))
load(here("output/tab_muris/liver_sideref.RData"))


### Loading GeneFishing Modules ----------------------------
# 
# ## Source GF Functions and Report Tightness of Bait Genes
# gf_sourcers <- list.files(here("../GeneFishing/R/"))
# 
# for(g in gf_sourcers) {
#   source(here("../GeneFishing/R/" %p% g))
# }
# 
# 

library(GeneFishing)



heart_aorta_dist_df <- 
  getDistDf(group_labels = heart_aorta_cell_names$cell_type, 
            dist_mat = heart_aorta_sideref_res$dissim_final)



fibroblast_cell_dist_df <-
  computeCellToGroupDists(ref_group = "fibroblast", 
                          group_labels = heart_aorta_cell_names$cell_type,
                          dist_mat = heart_aorta_sideref_res$dissim_final)




## Heatmaps of cell subgroups

liver_tab_muris_sideref_group_heatmap <-
  groupwiseDistanceHeatmap(
    group_labels = liver_cell_names$cell_type, 
    dist_mat = liver_sideref_res$dissim_final, 
    title = "SIDEREF Tabula Muris Liver Cells",
    hclust = TRUE)

liver_tab_muris_sideref_group_heatmap

heart_aorta_tab_muris_sideref_group_heatmap <-
  groupwiseDistanceHeatmap(
    group_labels = heart_aorta_cell_names$cell_type, 
    dist_mat = heart_aorta_sideref_res$dissim_final, 
    title = "SIDEREF Tabula Muris Heart Aorta Cells",
    hclust = TRUE)

heart_aorta_tab_muris_sideref_group_heatmap

## Liver, cholesterol metabolism search.
# probe_res <-
#   GeneFishing::probeFishability(liver_mat,
#                                 potential_bait = bait_genes_chol_metab,
#                                 n_rounds = 50,
#                                 ncores = 5)
# 
# 
#
#
#
# fibro_probe_res <-
#   GeneFishing::probeFishability(heart_aorta_mat,
#                                 potential_bait = bait_genes_fibro_prolif,
#                                 n_rounds = 50,
#                                 ncores = 2, 
#                                 min_tightness = 0.5)
# 
# doParallel::stopImplicitCluster()
# 
# save(fibro_probe_res, file = here('output/tab_muris/fibro_probe_res.RData'))
#
# 
# set.seed(25)
# min_tightness = 1
# fibro_prolif_probe_res <-
#   GeneFishing::probeFishability(heart_aorta_mat,
#                                 potential_bait = bait_genes_fibro_prolif,
#                                 n_rounds = 50,
#                                 ncores = 4, 
#                                 min_tightness = min_tightness)
# 
# 
# doParallel::stopImplicitCluster()
# 
# set.seed(25)
# min_tightness = 1
# fibro_prolif_probe_res2 <-
#   probeFishability(heart_aorta_mat,
#                                 potential_bait = bait_genes_fibro_prolif,
#                                 n_rounds = 50,
#                                 ncores = 4,
#                                 min_tightness = min_tightness)
# 
# doParallel::stopImplicitCluster()

## create null object
# set.seed(25)
# min_tightness = 0.01
# fibro_prolif_probe_res_null <-
#   GeneFishing::probeFishability(heart_aorta_mat,
#                                 potential_bait = bait_genes_fibro_prolif,
#                                 n_rounds = 2,
#                                 ncores = 4, 
#                                 min_tightness = min_tightness)
# 
# 
# doParallel::stopImplicitCluster()


# gf_res <-
#   GeneFishing::geneFishing(liver_mat,
#                            bait_genes = res$best_bait,
#                            fishing_rounds = 50,
#                            n_probing_rounds = 50)
#
#
# plot(gf_res)


## internal functions needed for the following:
source(here("../GeneFishing/R/computeAvgDBIndexUMAP.R"))
source(here("../GeneFishing/R/getUMAPCoordinates.R"))

computeTightnessResults <- 
  function(bait_genes, 
           good_bait,
           main_cell_group, 
           exp_mat,
           cell_types,
           group_dist_df,
           probe_each_time = TRUE,
           limit_samp_size = FALSE,
           n_rounds = 25,
           ncores = 5,
           min_tightness = 0.5) {
    
    stopifnot(main_cell_group %in% cell_types)
    
    dist_sub = group_dist_df %>% 
      dplyr::filter(group2 == main_cell_group & group != main_cell_group)
    
    n_main_cell_group = sum(cell_types == main_cell_group)
    
    cell_group_list = c(main_cell_group, 
                        dist_sub %>% 
                          arrange(dist) %>% pull(group),
                        "all")
    
    tightness_res_vec = vector(mode = "list", length= length(unique(cell_types))+1)
    names(tightness_res_vec) = cell_group_list
    
    i=1
    for(c in cell_group_list) {
      print(c)
      if(c == "all"){
        cell_indices = seq_len(ncol(exp_mat))
      } else{
        cell_indices = which(cell_types %in% c(main_cell_group, c)) 
      }
      
      ## limit each subset to be same number of cells, if desired
      if(limit_samp_size) {
        cell_indices = base::sample(cell_indices, n_main_cell_group, replace = FALSE)
        print(paste0(length(cell_indices), " cells"))
      }
      
      exp_mat_subset = exp_mat[,cell_indices]
      
      
      
      ## need to remove genes with zero expression among the subgroups
      ## TODO: could be considered a better preprocessing step to take away all genes that don't have some shared expression within each subgroup
      exp_mat_subset = exp_mat_subset[which(rowVars(exp_mat_subset) > 0), ]
      
      if(probe_each_time==TRUE) {
        subgroup_bait_genes <- bait_genes[which(bait_genes %in% rownames(exp_mat_subset))]
        
        subgroup_bait_genes <- 
          GeneFishing::probeFishability(exp_mat_subset,
                                        potential_bait = subgroup_bait_genes,
                                        n_rounds = n_rounds,
                                        ncores = ncores,
                                        min_tightness = min_tightness)
        
        doParallel::stopImplicitCluster()
        tightness_res_vec[[i]]$bait_status = "found bait"
        
        subgroup_bait_genes <- subgroup_bait_genes$best_bait
        
        
        ## offered good bait selection.
        if(is.null(subgroup_bait_genes)) {
          subgroup_bait_genes = good_bait
          tightness_res_vec[[i]]$bait_status = "default bait"
        }
        
        
        
        
        
      } else{
        subgroup_bait_genes = bait_genes
        tightness_res_vec[[i]]$bait_status = "all bait"
      }
      
      tightness_res_vec[[i]]$bait_genes = subgroup_bait_genes
      
      ## compute tightness.
      tightness_res <- 
        computeAvgDBIndexUMAP(bait_genes = subgroup_bait_genes, 
                              exp_mat = exp_mat_subset, 
                              method = "spearman", 
                              k = 2, alpha = 5, 
                              n_rounds = n_rounds, 
                              n_neighbors = 15)
      
      tightness_res_vec[[i]]$mean_tightness =
        sapply(tightness_res, function(x) x) %>% mean()
      
      
      
      i = i+1
    }
    return(tightness_res_vec)
  }


## Custom cell groups for a single run:
computeTightnessResultsByDist <- 
  function(result_length,
           bait_genes,
           exp_mat,
           main_cell_group,
           cell_types,
           cell_dists,
           n_rounds = 25,
           ncores = 4,
           min_tightness = 0.5,
           decreasing_dist = TRUE) {
    
    stopifnot(main_cell_group %in% cell_types)
    
    n_main_cell_group = sum(cell_types == main_cell_group)
    
    main_cells = which(cell_types %in% main_cell_group)
    
    ordered_cells <-
      as.integer(names(base::sort(cell_dists, decreasing = decreasing_dist)))
    
    tightness_res_vec = base::vector(mode = "numeric", length = result_length)
    names(tightness_res_vec) = ordered_cells[1:result_length]
    
    i=1
    extra_cells = c()
    for(o in ordered_cells[1:result_length]) {
      print(i)
      extra_cells = c(extra_cells, o)
      cell_indices = c(main_cells, extra_cells)
      
      exp_mat_subset = exp_mat[,cell_indices]
      
      
      ## need to remove genes with zero expression among the subgroups
      ## TODO: could be considered a better preprocessing step to take away all genes that don't have some shared expression within each subgroup
      exp_mat_subset = exp_mat_subset[which(rowVars(exp_mat_subset) > 0), ]
 
      ## compute tightness.
      tightness_res <- 
        computeAvgDBIndexUMAP(bait_genes = bait_genes, 
                              exp_mat = exp_mat_subset, 
                              method = "spearman", 
                              k = 2, alpha = 5, 
                              n_rounds = n_rounds, 
                              n_neighbors = 15)
      
      tightness_res_vec[[i]] =
        sapply(tightness_res, function(x) x) %>% mean()
      
      i = i+1
    }
    return(tightness_res_vec)
  }


# bycell_tight_res <-
#   computeTightnessResultsByDist(result_length = 300,
#                                 bait_genes = fibro_prolif_probe_res$bait_sets[[2]],
#                                 exp_mat = heart_aorta_mat,
#                                 main_cell_group = "fibroblast",
#                                 cell_types = heart_aorta_cell_names$cell_type,
#                                 cell_dists = fibroblast_cell_dist_df,
#                                 n_rounds = 50,
#                                 ncores = 4,
#                                 min_tightness = 0.5,
#                                 decreasing_dist = TRUE)
# 
# 
# 
# 
# tight_res = 
#   computeTightnessResults(bait_genes = bait_genes_fibro_prolif, 
#                           good_bait = fibro_prolif_probe_res$best_bait,
#                           main_cell_group = "fibroblast", 
#                           exp_mat = heart_aorta_mat,
#                           cell_types = heart_aorta_cell_names$cell_type,
#                           group_dist_df = heart_aorta_dist_df,
#                           probe_each_time = TRUE,
#                           n_rounds = 50,
#                           ncores = 4)
#     
# 
# 
# 
# tight_res_same_n_cells = 
#   computeTightnessResults(bait_genes = fibro_prolif_probe_res$bait_sets[[2]], 
#                           good_bait = fibro_prolif_probe_res$bait_sets[[2]],
#                           main_cell_group = "fibroblast", 
#                           exp_mat = heart_aorta_mat,
#                           cell_types = heart_aorta_cell_names$cell_type,
#                           group_dist_df = heart_aorta_dist_df,
#                           probe_each_time = FALSE,
#                           limit_samp_size = TRUE,
#                           n_rounds = 50,
#                           ncores = 4)



## TODO: separate modules with actual data analysis.


### Load Results and plot

LoadBuildCellTightnessDf <- function(out_path, file_pattern, n_cells) {
  
  bycell_tight_res_df <-
    data.frame(n_cells = seq_len(n_cells))
  
  for(i in list.files(path=here(out_path), 
                      pattern = 'cell_tight_res_' %p% file_pattern)) {
    load(here(out_path %p% i))
    
    join_data <- 
      data.frame(n_cells = seq_len(length(bycell_tight_res)))
    
    join_data[gsub("cell_tight_res_" %p% file_pattern %p% "_", "", 
                   gsub("\\.RData", "", i))] <- bycell_tight_res
    
    bycell_tight_res_df <- 
      bycell_tight_res_df %>% 
      left_join(join_data)
  }
  return(bycell_tight_res_df)
}

cardio_bycell_tight_res_df <- 
  LoadBuildCellTightnessDf(out_path = 'output/tab_muris/',
                           file_pattern = "cardio",
                           n_cells = ncol(heart_aorta_mat))

fibro_bycell_tight_res_df <- 
  LoadBuildCellTightnessDf(out_path = 'output/tab_muris/',
                           file_pattern = "fibro",
                           n_cells = ncol(heart_aorta_mat))


## plots:
cardio_bycell_tight_res_df %>% 
  gather(method, tightness, setdiff(names(cardio_bycell_tight_res_df), "n_cells")) %>% 
  dplyr::filter(!is.na(tightness)) %>%
  dplyr::filter(method != "decr_dist_sideref_full") %>% 
  mutate(method = case_when(
    method == "decr_dist" ~ "Decreasing Distance\n(SIDEREF on Relevant Genes)",
    method == "incr_dist" ~ "Increasing Distance\n(SIDEREF on Relevant Genes)",
    method == "random"    ~ "Random Cells",
    method == "decr_dist_sideref_full" ~ "Decreasing Distance\n(SIDEREF on All Variable Genes)"
  )) %>%
  ggplot(aes(x=n_cells, y =tightness, col=method)) + 
  geom_point(size=0.2) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  labs(x = "Number cells added to Main Subgroup", 
       y = "Tightness", 
       col = "Distance Criteria") + 
  guides(color=guide_legend(override.aes=list(size=2))) + 
  ggtitle("Cardiac Muscle Contraction Probed Genes\n" %p% 
            "Cells Added to Cardiac Muscle Cell Group")


ggsave(here("output/figures/tab_muris/bycell_tightness_cardio.png"), 
       height = 6, width = 10)



fibro_bycell_tight_res_df %>% 
  gather(method, tightness, setdiff(names(fibro_bycell_tight_res_df), "n_cells")) %>%
  dplyr::filter(!is.na(tightness)) %>% 
  #dplyr::filter(method != "decr_dist_sideref_full") %>% 
  mutate(method = case_when(
    method == "decr_dist" ~ "Decreasing Distance\n(SIDEREF on Relevant Genes)",
    method == "incr_dist" ~ "Increasing Distance\n(SIDEREF on Relevant Genes)",
    method == "random"    ~ "Random Cells",
    method == "decr_dist_sideref_full" ~ "Decreasing Distance\n(SIDEREF on All Variable Genes)"
  )) %>%
  ggplot(aes(x=n_cells, y =tightness, col=method)) + 
  geom_point(size=0.2) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  labs(x = "Number cells added to Main Subgroup", 
       y = "Tightness", 
       col = "Distance Criteria") + 
  guides(color=guide_legend(override.aes=list(size=2))) + 
  ggtitle("Fibroblast Proliferation Probed Genes\n" %p% 
            "Cells Added to Fibroblast Cell Group")


ggsave(here("output/figures/tab_muris/bycell_tightness_fibro.png"), 
       height = 6, width = 10)






