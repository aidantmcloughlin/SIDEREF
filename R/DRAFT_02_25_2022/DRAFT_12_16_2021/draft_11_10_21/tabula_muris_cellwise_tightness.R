### Tabula Muris Analysis (cell-wise tightness impact)

raw_data_loc <- "data/tab_muris_bulk/"
output_data_loc <- "output/tab_muris/"

## Load datasets 

liver_df <- 
  readRDS(here(raw_data_loc %p% "Liver.rds"))

heart_aorta_df <- 
  readRDS(here(raw_data_loc %p% "Heart_and_Aorta.rds"))

### clean data by size counts
liver_mat <- 
  as.matrix(t(t(liver_df$counts) / liver_df$size_factor))

heart_aorta_mat <- 
  as.matrix(t(t(heart_aorta_df$counts) / heart_aorta_df$size_factor))



### Cell type metadata

## Heart Aorta
load(here(raw_data_loc %p% "droplet_Heart_and_Aorta_seurat_tiss.Robj"))

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
load(here(raw_data_loc %p% "droplet_Liver_seurat_tiss.Robj"))

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
  fread(here(raw_data_loc %p%
               "GO_term_cholesterol_metabolism.csv"))


## heart, aorta genes

GO_card_musc_contract <- 
  fread(here(raw_data_loc %p% 
               "GO_term_cardiac_muscle_contraction.csv"))

GO_antigen_pres <- 
  fread(here(raw_data_loc %p% 
               "GO_term_antigen_presentation.csv"))

GO_fibro_prolif <- 
  fread(here(raw_data_loc %p% 
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



## Load SIDEREF Results and Plot Heatmap
source(here("R/subgroup_composition_gene_contrib.R"))


## internal functions needed for the following tightness computations:
source(here("../GeneFishing/R/computeAvgDBIndexUMAP.R"))
source(here("../GeneFishing/R/getUMAPCoordinates.R"))




## Custom cell groups for a single run:
computeTightnessResultsByDist <- 
  function(result_length = NULL,
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
    
    if(is.null(result_length)) {
      result_length = length(ordered_cells)
    }
    
    for(o in ordered_cells[1:result_length]) {
      print(i %p% " of " %p% result_length)
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








## load SIDEREF dist RES
print("loading sideref dist mats")
load(here(output_data_loc %p% "cardio_heart_aorta_sideref.RData"))
load(here(output_data_loc %p% "fibro_heart_aorta_sideref.RData"))


### Get pathway specific cell distances
fibroblast_cell_dist_df <-
  computeCellToGroupDists(ref_group = "fibroblast", 
                          group_labels = heart_aorta_cell_names$cell_type,
                          dist_mat = fibro_heart_aorta_sideref_res$dissim_final)


cardio_cell_dist_df <-
  computeCellToGroupDists(ref_group = "cardiac muscle cell", 
                          group_labels = heart_aorta_cell_names$cell_type,
                          dist_mat = cardio_heart_aorta_sideref_res$dissim_final)


## Load the Bait Sets
load(here(output_data_loc %p% "cardiac_musc_probe_res.RData"))
load(here(output_data_loc %p% "fibro_probe_res.RData"))


print("running tightness computations")


## FIBRO DECREASE
bycell_tight_res <-
  computeTightnessResultsByDist(result_length = 300,
                                bait_genes = fibro_probe_res$best_bait,
                                exp_mat = heart_aorta_mat,
                                main_cell_group = "fibroblast",
                                cell_types = heart_aorta_cell_names$cell_type,
                                cell_dists = fibroblast_cell_dist_df,
                                n_rounds = 50,
                                ncores = 1,
                                min_tightness = 0.5,
                                decreasing_dist = TRUE)


save(bycell_tight_res, file = "cell_tight_res_fibro_decr_dist.RData")


## CARDIAC INCREASE
bycell_tight_res <-
  computeTightnessResultsByDist(result_length = 300,
                                bait_genes = cardiac_musc_probe_res$best_bait,
                                exp_mat = heart_aorta_mat,
                                main_cell_group = "cardiac muscle cell",
                                cell_types = heart_aorta_cell_names$cell_type,
                                cell_dists = cardio_cell_dist_df,
                                n_rounds = 50,
                                ncores = 1,
                                min_tightness = 0.5,
                                decreasing_dist = FALSE)


save(bycell_tight_res, file = "cell_tight_res_cardio_incr_dist.RData")



## CARDIAC DECREASE
bycell_tight_res <-
  computeTightnessResultsByDist(result_length = 300,
                                bait_genes = cardiac_musc_probe_res$best_bait,
                                exp_mat = heart_aorta_mat,
                                main_cell_group = "cardiac muscle cell",
                                cell_types = heart_aorta_cell_names$cell_type,
                                cell_dists = cardio_cell_dist_df,
                                n_rounds = 50,
                                ncores = 1,
                                min_tightness = 0.5,
                                decreasing_dist = TRUE)


save(bycell_tight_res, file = "cell_tight_res_cardio_decr_dist.RData")


## FIBRO INCREASE
bycell_tight_res <-
  computeTightnessResultsByDist(result_length = 300,
                                bait_genes = fibro_probe_res$best_bait,
                                exp_mat = heart_aorta_mat,
                                main_cell_group = "fibroblast",
                                cell_types = heart_aorta_cell_names$cell_type,
                                cell_dists = fibroblast_cell_dist_df,
                                n_rounds = 50,
                                ncores = 1,
                                min_tightness = 0.5,
                                decreasing_dist = FALSE)


save(bycell_tight_res, file = "cell_tight_res_fibro_incr_dist.RData")



#### RANDOM RESULTS:
fibroblast_cell_dist_df_rand = fibroblast_cell_dist_df
fibroblast_cell_dist_df_rand[1:length(fibroblast_cell_dist_df_rand)] = 
  runif(n = length(fibroblast_cell_dist_df_rand))

bycell_tight_res_random <-
  computeTightnessResultsByDist(result_length = 75,
                                bait_genes = fibro_prolif_probe_res$bait_sets[[2]],
                                exp_mat = heart_aorta_mat,
                                main_cell_group = "fibroblast",
                                cell_types = heart_aorta_cell_names$cell_type,
                                cell_dists = fibroblast_cell_dist_df_rand,
                                n_rounds = 50,
                                ncores = 4,
                                min_tightness = 0.5,
                                decreasing_dist = TRUE)



########################################## Plot Tightness Results


bycell_tight_res_df <- 
  data.frame(n_cells = seq_len(length(bycell_tight_res)),
             mean_tightness = bycell_tight_res,
             distance = "SIDEREF Distance to Fibroblasts")


bycell_tight_res_df_rand <- 
  data.frame(n_cells = seq_len(length(bycell_tight_res_random)),
             mean_tightness = bycell_tight_res_random,
             distance = "Random")

bind_rows(load(output_data_loc %p% "cell_tight_res_fibro_decr_dist.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "SIDEREF Decr Dist to Fibroblasts"),
          
          load(output_data_loc %p% "cell_tight_res_fibro_incr_dist.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "SIDEREF Incr Dist to Fibroblasts"),
          
          load(output_data_loc %p% "cell_tight_res_fibro_decr_dist_sideref_full.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "SIDEREF Decr Dist to Fibroblasts\nAll Genes"),
          load(output_data_loc %p% "cell_tight_res_fibro_random.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "Random Cells")
          ) %>% 
  ggplot(aes(x = n_cells, y = mean_tightness, color = distance)) + 
  geom_point(size = 0.2) + 
  theme_bw() + 
  labs(x = "Number of Added Cells",
       y = "Mean Tightness",
       color = "Cell Inclusion Method") + 
  ggtitle("Fibroblast Proliferation Genes\nCells added by Distance to Fibroblast Cell Type")





bind_rows(load(output_data_loc %p% "cell_tight_res_cardio_decr_dist.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "SIDEREF Decr Dist to Cardiac Muscle"),
          
          load(output_data_loc %p% "cell_tight_res_cardio_incr_dist.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "SIDEREF Incr Dist to Cardiac Muscle"),
          
          # load(output_data_loc %p% "cell_tight_res_cardio_decr_dist_sideref_full.RData") %>% 
          #   data.frame(n_cells = seq_len(length(bycell_tight_res)),
          #              mean_tightness = bycell_tight_res,
          #              distance = "SIDEREF Decr Dist to Cardiac Muscle\nAll Genes"),
          
          load(output_data_loc %p% "cell_tight_res_cardio_random.RData") %>% 
            data.frame(n_cells = seq_len(length(bycell_tight_res)),
                       mean_tightness = bycell_tight_res,
                       distance = "Random Cells")
) %>% 
  ggplot(aes(x = n_cells, y = mean_tightness, color = distance)) + 
  geom_point(size = 0.02) + 
  theme_bw() + 
  labs(x = "Number of Added Cells",
       y = "Mean Tightness",
       color = "Cell Inclusion Method") + 
  ggtitle("Cardiac Muscle Contraction Genes\nCells added by Distance to Cardiac Muscle Cell Type")



### TODO: RUN subgroup-wise inclusion on cardiac muscles to see if it conforms with cell-wise result we see.





### let's look at underlying UMAP plot
umap_res <- uwot::umap(as.dist(fibro_heart_aorta_sideref_res$dissim_final))

far_cells = as.integer(names(base::sort(fibroblast_cell_dist_df, decreasing = TRUE)))[1:100]

near_cells = as.integer(names(base::sort(fibroblast_cell_dist_df, decreasing = FALSE)))[1:100]


umap_by_side_dist <-
  data.frame(umap_res) %>% 
  mutate(cell_num = row_number(),
         cell_cat = case_when(
           cell_num %in% far_cells ~ "far cell",
           cell_num %in% near_cells ~ "near cell",
           cell_num %in% which(heart_aorta_cell_names$cell_type == "fibroblast") ~ "fibroblast",
           TRUE ~ "other cell")) %>% 
  ggplot(aes(x=X1, y=X2, color=cell_cat)) + 
  geom_point() + 
  labs(x = "UMAP1", y = "UMAP2", color = "Cell Category") + 
  theme_bw()

umap_by_side_dist



umap_by_cell_type <-
  data.frame(umap_res) %>% 
  mutate(cell_num = row_number(),
         cell_cat = heart_aorta_cell_names$cell_type) %>% 
  ggplot(aes(x=X1, y=X2, color=cell_cat)) + 
  geom_point() + 
  labs(x = "UMAP1", y = "UMAP2", color = "Cell Category") + 
  theme_bw()

umap_by_cell_type


## TODO UMAP PLOTS OF BAIT GENES WITH RANDOM GENES.



### TODO: more complicated splatter simulation run and ranking results

## TODO: some brainstorming regarding what to do next about cell selection.



