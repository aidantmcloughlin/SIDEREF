library(here)
source(here("R/libraries.R"))
source(here('R/relative_group_dist_comps.R'))

############ PERFORMANCE METRICS FUNCTIONS ###############

#########################################################################
### Global group violations over multiple cell type grouping levels.
#########################################################################



countGlobalGrpMaxViolationsMulti <- function(
  dist, 
  group_labels, 
  group_global_group_map) {
  ### group_global_group_map contains columns: 
  ##   (1) cell_group_low
  ##   columns with the following naming convention:
  ##   (2) cell_group_high1, cell_group_high2, ..., cell_group_highM
  ##      later columns are coarser groupings of prior columns.
  
  ## TODO: error raising condition on cell_group_high columns.
  
  
  n_g = nrow(group_global_group_map)
  n_multi = sum(grepl('cell_group_high', names(group_global_group_map)))
  ### Generate Relative Distances DF (e.g. SIDEREF):
  grp_dist_mat <- 
    getDistDf(
      group_labels = group_labels, 
      dist_mat = dist,
      symmetrize = FALSE,
      gather = FALSE) %>% 
    dplyr::rename(cell_group_low = target_group) %>% 
    left_join(
      group_global_group_map %>%
        group_by(cell_group_high1) %>% 
        mutate(global_group_size = n()) %>%
        ungroup() %>%
        mutate(cell_group_low = as.character(cell_group_low)) %>%
        dplyr::select(cell_group_low, contains("cell_group_high"), global_group_size))
  
  ## loop thru rows and count violations
  glob_grp_vio_count <- 0
  glob_grp_nonvio_count <- 0
  for(i in seq_len(n_g)) {
    
    if(grp_dist_mat$global_group_size[i] > 1) {
      row_vio <- rep(FALSE, n_g)
      col_vio <- rep(FALSE, n_g)
      for(m in seq_len(n_multi)) {
        glob_grp_idx = which(grp_dist_mat[,"cell_group_high" %p% m] == 
                               grp_dist_mat[i, "cell_group_high" %p% m])
        
        ### one "cross" at a time:
        max_within_dist = max(
          max(grp_dist_mat[glob_grp_idx, i]),
          max(grp_dist_mat[i, glob_grp_idx]))
          
        
        row_vio[setdiff(1:n_g, glob_grp_idx)] <- 
          pmax(row_vio[setdiff(1:n_g, glob_grp_idx)],
               grp_dist_mat[setdiff(1:n_g, glob_grp_idx), i] <= max_within_dist)
        col_vio[setdiff(1:n_g, glob_grp_idx)] <- 
          pmax(col_vio[setdiff(1:n_g, glob_grp_idx)],
               grp_dist_mat[i, setdiff(1:n_g, glob_grp_idx)] <= max_within_dist)
        
      }
      ### subset to what has actually been evaluated
      granular_idx <- which(grp_dist_mat[,"cell_group_high1"] == 
                               grp_dist_mat[i, "cell_group_high1"])
      row_vio <- row_vio[setdiff(1:n_g, granular_idx)]
      col_vio <- col_vio[setdiff(1:n_g, granular_idx)]
      row_nonvio = ifelse(row_vio == 1, 0, 1)
      col_nonvio = ifelse(col_vio == 1, 0, 1)
          
      ## tabulate the violations 
      glob_grp_vio_count <- glob_grp_vio_count + 
        sum(row_vio) + sum(col_vio) 
      
      glob_grp_nonvio_count <- glob_grp_nonvio_count +
        sum(row_nonvio) + sum(col_nonvio) 
            
      
    }
  }
  ## return success rate:
  return(glob_grp_vio_count / (glob_grp_nonvio_count + glob_grp_vio_count))
  
}




countGlobalGrpMinViolationsMulti <- function(
  dist, 
  group_labels, 
  group_global_group_map) {
  ### group_global_group_map contains columns: 
  ##   (1) cell_group_low
  ##   columns with the following naming convention:
  ##   (2) cell_group_high1, cell_group_high2, ..., cell_group_highM
  ##      later columns are coarser groupings of prior columns.
  
  ## TODO: error raising condition on cell_group_high columns.
  
  
  n_g = nrow(group_global_group_map)
  n_multi = sum(grepl('cell_group_high', names(group_global_group_map)))
  ### Generate Relative Distances DF (e.g. SIDEREF):
  grp_dist_mat <- 
    getDistDf(
      group_labels = group_labels, 
      dist_mat = dist,
      symmetrize = FALSE,
      gather = FALSE) %>% 
    dplyr::rename(cell_group_low = target_group) %>% 
    left_join(
      group_global_group_map %>%
        group_by(cell_group_high1) %>% 
        mutate(global_group_size = n()) %>%
        ungroup() %>%
        mutate(cell_group_low = as.character(cell_group_low)) %>%
        dplyr::select(cell_group_low, contains("cell_group_high"), global_group_size))
  
  ## loop thru rows and count violations
  glob_grp_vio_count <- 0
  glob_grp_nonvio_count <- 0
  for(i in seq_len(n_g)) {
    
    if(grp_dist_mat$global_group_size[i] > 1) {
      row_vio <- rep(FALSE, n_g)
      col_vio <- rep(FALSE, n_g)
      for(m in seq_len(n_multi)) {
        glob_grp_idx = which(grp_dist_mat[,"cell_group_high" %p% m] == 
                               grp_dist_mat[i, "cell_group_high" %p% m])
        
        ### one "cross" at a time:
        min_without_dist = min(
          min(grp_dist_mat[setdiff(1:n_g, glob_grp_idx), i]),
          min(grp_dist_mat[i, setdiff(1:n_g, glob_grp_idx)]))
        
        
        row_vio[glob_grp_idx] <- 
          pmax(row_vio[glob_grp_idx],
               grp_dist_mat[glob_grp_idx, i] >= min_without_dist)
        col_vio[glob_grp_idx] <- 
          pmax(col_vio[glob_grp_idx],
               grp_dist_mat[i, glob_grp_idx] >= min_without_dist)
        
      }
      ### subset to what has actually been evaluated
      coarse_idx <- which(grp_dist_mat[,"cell_group_high" %p% n_multi] == 
                            grp_dist_mat[i, "cell_group_high" %p% n_multi])
      row_vio <- row_vio[coarse_idx]
      col_vio <- col_vio[coarse_idx]
      row_nonvio = ifelse(row_vio == 1, 0, 1)
      col_nonvio = ifelse(col_vio == 1, 0, 1)
      
      ## tabulate the violations 
      glob_grp_vio_count <- glob_grp_vio_count + 
        sum(row_vio) + sum(col_vio) 
      
      glob_grp_nonvio_count <- glob_grp_nonvio_count +
        sum(row_nonvio) + sum(col_nonvio) 
      
      
    }
  }
  ## return success rate:
  return(glob_grp_vio_count / (glob_grp_nonvio_count + glob_grp_vio_count))
  
}



