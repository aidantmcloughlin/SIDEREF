###############################################################################
### Gene Set Enrichment for B-Cells within Kidney, Thymus Leukocytes.
###############################################################################

library(here)
source(here("R/libraries.R"))
source(here("R/tm_helper_functions.R"))
source(here("R/distance_compute_helpers.R"))

set.seed(1)


## Load variable genes full data set to subset via metadata indices.
load(here("output/tab_muris_sc/preproc_data/tab_muris_full.RData"),
     verbose=TRUE)

## Get metadata via helper function
load(
  getTMCellTypeSampleDistFileName(
    file_tag = "meta_data",
    file_loc = here("output/tab_muris_sc/dist_res/"),
    samp_size = TM_SAMP_SIZE,
    use_droplet_only = USE_DROPLET_ONLY),
  verbose = TRUE)

rownames(meta_data_full) <- meta_data_full$cell_id


data_full <- data_full[, meta_data_full$cell_id]

## Convert to Seurat
tm_preproc <-
  TabMurisProcSubset(full_data = data_full, 
                     meta_data = meta_data_full, 
                     var_features= VAR_FEATURES,
                     return_seurat = TRUE)[[1]]

## DE Genes among the B-cell groups
getMarkerGSList <- function(assay, cells_in, cells_vs, 
                            cell_marker_list_sizes = 
                              c(10, 50, 100, 300),
                            type_name = "B_markers",
                            logfc_thres = 0.25) {
  
  
  cell_marker_res <- 
    FindMarkers(assay,
                cells.1 = cells_in,
                cells.2 = setdiff(cells_vs, cells_in),
                ## speeds up computation while we are looking for small
                ##  sets of DE gene cells in our main reference list.
                logfc.threshold = logfc_thres)
  
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
prepRanks <- function(de_marker_res, all_gene_names = c(),
                      keep_stat = FALSE) {
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
                  test.use = "t",
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
b_cells     <- meta_data_full$cell_id[
  which(meta_data_full$cell_type == "B cell")]

leuko_cells <- meta_data_full$cell_id[
  which(meta_data_full$cell_type == "leukocyte")]

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
       width = 8, height = 4,
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

gsea_table%>% 
  write.csv(here("output/tab_muris_sc/gse_pvals.csv"))

