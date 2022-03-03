###############################################################################
### Load, Cut TM Data.
###############################################################################

library(here)
source(here('R/libraries.R'))

tmLoadCut <- function(file_name, data_pct, 
                      reps, folds, sample_size = 20,
                      keep_droplet = TRUE,
                      keep_facs = TRUE,
                      sample_cell_types = FALSE,
                      sample_pct = TRUE,
                      seed = 1) {
  
  set.seed(seed)
  ## set seeds for each fold
  if(sample_pct == FALSE) {
    reps = 1
    folds = 1
  }
  
  fold_seeds <- 100:(100+reps-1)
  
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
  
  
  ## Evenly split each file source (or celltpye).
  if(sample_cell_types == FALSE) {
    
    group_based_idx <- lapply(unique(meta_data_full$file_source),
                             function(x) which(meta_data_full$file_source == x))
  } else{
    group_based_idx <- lapply(unique(meta_data_full$cell_type %p% meta_data_full$tissue),
                             function(x) which(meta_data_full$cell_type %p% meta_data_full$tissue == x))
    
  }
  
  cell_idx_splits <- vector(mode="list", length=folds)
  
  if(sample_pct == FALSE) {
    group_idx_shuf <-
      lapply(group_based_idx, function(x){
        sample(x, size = min(sample_size, length(x)), replace=FALSE)}) 
    
    
    for(g in seq_len(length(group_idx_shuf))) {
      cell_idx_splits[[1]] <- c(cell_idx_splits[[1]], 
                                group_idx_shuf[[g]])
    }
    
    
  } else{
    group_idx_shuf <-
      lapply(group_based_idx, function(x){
        sample(x, replace=FALSE)})
    
    cell_idx_splits_nested <- lapply(group_idx_shuf, function(x){
      split(x,
            cut(seq_along(x), folds,
                labels = FALSE))
      
    })
    
    for(f in seq_len(folds)) {
      for(g in seq_len(length(cell_idx_splits_nested))) {
        cell_idx_splits[[f]] = c(cell_idx_splits[[f]], 
                                 cell_idx_splits_nested[[g]][[f]])
      }
    }
  }
  
  
  data_store_list <- vector(mode="list", length=reps)
  data_inner_list <- vector(mode="list", length=3)
  
  
  for(f in seq_len(reps)) {
    data_store_list[[f]] <- data_inner_list
    cells <- cell_idx_splits[[f]]
    
    data_store_list[[f]][[1]] <- f
    data_subset <- data_full[, cells]
    meta_data_subset <- meta_data_full[cells, ]
    
    data_store_list[[f]][[2]] <- data_subset
    
    data_store_list[[f]][[3]] <- meta_data_subset
    
    data_store_list[[f]][[4]] <- fold_seeds[f]
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



