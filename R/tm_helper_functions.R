loadTMCellTypeSampleDistFile <- function(file_tag,
                                         file_loc = here("output/tab_muris_sc/dist_res/"),
                                         samp_size = TM_SAMP_SIZE,
                                         use_droplet_only = USE_DROPLET_ONLY,
                                         load_meta_data = TRUE) {
  load(file_loc %p% "samp_ctTRUE_samp_pctFALSE_samp_size" %p% samp_size %p%
         "_" %p% file_tag %p% ".RData",
       verbose=TRUE)
  
  if(load_meta_data) {
    load(file_loc %p% "meta_data_samp_ctTRUE_samp_pctFALSE_samp_size" %p% 
           samp_size %p% ".RData",
         verbose=TRUE)
  }
  
  
}

prepTMGroupLabels <- 
  function(meta_data,
           include_tissue = TRUE) {
    
    
    if(include_tissue) {
      meta_data <- meta_data %>%
        mutate(group_label = cell_type %p% ", " %p% tissue)
    }
    
    ## cleaning of labels.
    meta_data$group_label <- str_to_title(gsub("_", " ", meta_data$group_label))
    ## title case corrections
    meta_data$group_label <- gsub(" And", " and", meta_data$group_label)
    meta_data$group_label <- gsub(" Of", " of", meta_data$group_label)
    meta_data$group_label <- gsub(" Ii", " II", meta_data$group_label)
    
    meta_data$group_label <- as.factor(meta_data$group_label)
    
    return(meta_data)
  }







