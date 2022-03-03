source(here("R/distance_compute_helpers.R"))

getSubDistDf <- 
  function(dist,
           meta_data,
           cell_id_vec = NULL) {
    
    dist <-
      selectFullDistIdx(
        cell_id_vec = cell_id_vec,
        dist = dist,
        meta_data)
    group_labels <- getGroupLabels(meta_data = meta_data, 
                                   cell_id_vec = cell_id_vec)
    
    dist_df <- getDistDf(group_labels, 
                         dist, 
                         gather = TRUE)
    
    return(dist_df)
  }


sampleNPerClust <- 
  function(
    assay_data,
    meta_data, 
    n
  ) {
    group_based_idx <- 
      lapply(unique(meta_data$group_label),
             function(x) which(meta_data$group_label == x))
    
    
    group_idx_shuf <-
      lapply(group_based_idx, function(x){
        sample(x, size = min(n, length(x)), replace=FALSE)})
    
    cell_idx <- c()
    
    for(g in seq_len(length(group_idx_shuf))) {
      cell_idx <- c(cell_idx, 
                    group_idx_shuf[[g]])
    }
    
    assay_data_sample <- assay_data[, cell_idx]
    meta_data_sample <- meta_data[cell_idx, ]
    
    return(list(assay_data_sample = assay_data_sample,
                meta_data_sample  = meta_data_sample))
    
  }

