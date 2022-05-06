###############################################################################
### Load, Cut TM Data.
###############################################################################

library(here)
source(here('R/libraries.R'))

tmLoadCut <- function(file_name, 
                      reps, sample_size = 20,
                      sample_cell_types = FALSE,
                      keep_droplet = TRUE,
                      keep_facs = TRUE,
                      seed = 1) {
  
  set.seed(seed)
  
  ### load Tabula Muris data
  load(here('output/tab_muris_sc/preproc_data/' %p% file_name %p% '.RData'),
       verbose = TRUE)
  
  ### prep input files 
  
  droplet_cells <- which(grepl("droplet", meta_data_full$file_source))
  facs_cells <- which(grepl("facs", meta_data_full$file_source))
  
  if(!keep_facs & keep_droplet) {
    data_full <- data_full[,droplet_cells]
    meta_data_full <- meta_data_full[droplet_cells,]
    
    file_tag = "droplet_only"
    
  } else if(keep_facs & !keep_droplet) {
    data_full <- data_full[,facs_cells]
    meta_data_full <- meta_data_full[facs_cells,]
    
    file_tag = "facs_only"
  } else{  
    file_tag = "droplet_and_facs"
  }
  
  
  tm_filenames <- unique(meta_data_full$file_source)
  
  
  ## Evenly split each file source (or celltype).
  if(sample_cell_types == FALSE) {
    
    group_based_idx <- lapply(unique(meta_data_full$file_source),
                             function(x) which(meta_data_full$file_source == x))
  } else{
    group_based_idx <- lapply(unique(meta_data_full$cell_type %p% meta_data_full$tissue),
                             function(x) which(meta_data_full$cell_type %p% 
                                                 meta_data_full$tissue == x))
    
  }
  
  cell_idx_splits <- vector(mode="list", length=reps)
  
  for(r in reps) {
    group_idx_shuf <-
      lapply(group_based_idx, function(x){
        sample(x, size = min(sample_size, length(x)), replace=FALSE)}) 
    
    for(g in seq_len(length(group_idx_shuf))) {
      cell_idx_splits[[r]] <- c(cell_idx_splits[[1]], 
                                group_idx_shuf[[g]])
    } 
  }
  
  
  data_store_list <- vector(mode="list", length=reps)
  data_inner_list <- vector(mode="list", length=3)
  rep_seeds <- 100:(100+reps-1)
  
  
  for(r in seq_len(reps)) {
    data_store_list[[r]] <- data_inner_list
    cells <- cell_idx_splits[[r]]
    
    data_store_list[[r]][[1]] <- r
    data_subset <- data_full[, cells]
    meta_data_subset <- meta_data_full[cells, ]
    
    data_store_list[[r]][[2]] <- data_subset
    
    data_store_list[[r]][[3]] <- meta_data_subset
    
    data_store_list[[r]][[4]] <- rep_seeds[r]
    
  }
  
  
  ### Store as output list.
  return(list(
    data_store_list = data_store_list,
    cell_idx_splits = cell_idx_splits,
    tm_filenames = tm_filenames,
    file_tag = file_tag,
    data_full = data_full,
    meta_data_full = meta_data_full,
    droplet_cells = droplet_cells,
    facs_cells = facs_cells))
  
}



