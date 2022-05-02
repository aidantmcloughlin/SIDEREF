


## average distances of distance matrix between cluster groups.
computeGroupDists <- function(ref_group, group_labels, dist_mat,
                              normalize = TRUE) {
  
  G <- length(unique(group_labels))

  ref_cells <- which(group_labels == ref_group)

  group_dists <- rep(NA, G)
  
  for(g in seq_len(G)) {
    target_cells <- which(group_labels == g)
    
    cell_dists <- rep(NA, length(ref_cells))
    i <- 1
    
    for(c in ref_cells) {
      clust_dists <- dist_mat[setdiff(target_cells, c), c]
      if(g == ref_group) {
        cell_dists[i] <- 1/(length(target_cells) - 1) * sum(clust_dists)
      } else{
        cell_dists[i] <- 1/(length(target_cells)) * sum(clust_dists)}
      
      i <- i + 1
    }
    
    group_dists[g] <- mean(cell_dists)
    
  }
  
  
  ## normalize
  if(normalize == TRUE) {
    group_dists <- (group_dists - min(group_dists)) / (max(group_dists) - min(group_dists)) 
  }
  
  return(group_dists)
}


computeCellToGroupDists <- function(ref_group, group_labels, dist_mat,
                                    normalize = TRUE) {
  
  ref_cells <- which(group_labels == ref_group)
  
  other_cells <- setdiff(seq_len(ncol(dist_mat)), ref_cells)
  
  cell_dists <- rep(NA, length(other_cells))
  names(cell_dists) = other_cells
  i <- 1
  for(o in other_cells) {
    cell_dists[i] = mean(dist_mat[ref_cells, o])
    
    i <- i+1
  }
  
  ## normalize
  if(normalize == TRUE) {
    cell_dists <- (cell_dists - min(cell_dists)) / (max(cell_dists) - min(cell_dists)) 
  }
  
  return(cell_dists)
}


### TODO: axis labeling functions could be in their own R module

hclustBasedAxisLabels <- function(group_labels, dist_mat,
                                  symmetrize = FALSE) {
  distance_df <- 
    getDistDf(group_labels, 
              dist_mat,
              symmetrize = symmetrize,
              gather = FALSE) %>% 
    dplyr::select(-target_group)
  
  distance_df_sym <- 
    distance_df
  
  distance_df_sym <- 0.5*(distance_df_sym + t(distance_df_sym))
  
  cluster     <- hclust(dist(distance_df_sym), method="ward.D2")
  dendrogram  <- as.dendrogram(cluster)
  Rowv        <- rowMeans(distance_df, na.rm = T)
  dendrogram  <- reorder(dendrogram, -Rowv)
  
  
  rowInd <- rev(order.dendrogram(dendrogram))
  colInd <- rowInd
  distance_df <- distance_df[rowInd, colInd]
  
  ## change levels
  distance_df <- 
    distance_df %>% 
    mutate(target_group = factor(names(distance_df)))
  
  distance_df$target_group <- 
    forcats::fct_inorder(distance_df$target_group)
  
  distance_df <- distance_df %>%
    gather(source_group, dist, 
           setdiff(names(distance_df), c("target_group")))
  
  ## change levels (source_group)
  distance_df$source_group <-
    factor(distance_df$source_group, 
           levels = levels(distance_df$target_group))
  
  return(distance_df)
}



presetAxisLabels <- function(distance_df, preset_levels) {
  
  ## TODO: error checking on preset_levels entry?!
  distance_df$target_group = factor(distance_df$target_group, 
                                    levels = preset_levels)
  distance_df$source_group = factor(distance_df$source_group, 
                                    levels = preset_levels)
  
  ## remove groups not in the preset levels.
  distance_df <- distance_df %>% 
    dplyr::filter(!is.na(target_group) & !is.na(source_group))
  
  return(distance_df)
}


### TODO: I think this function needs to be vetted.
### TODO: Debug setting of numeric axes labels.
numericAxisLabels <- function(
  distance_df, 
  do_numeric_x, 
  do_numeric_y) {
  
  if(do_numeric_x) {
    
    distance_df$target_group = 
      factor(as.character(as.numeric(distance_df$target_group)), 
             levels = as.character(seq_len(max(
               as.numeric(distance_df$target_group)))))
    
    
    distance_df <- distance_df %>% 
      dplyr::mutate(new_group = factor(
        source_group %p% ", " %p% as.numeric(source_group),
        levels = levels(source_group) %p% ", " %p% 
          seq_len(length(levels(source_group)))))
    
    if(do_numeric_y) {
      distance_df <- distance_df %>% 
        mutate(new_group = 
                 factor(as.character(as.numeric(new_group)),
                        levels = as.character(seq_len(
                          max(as.numeric(new_group))))
                 )
        )
    }
    distance_df <- distance_df %>% 
      mutate(source_group = new_group) %>% 
      dplyr::select(-new_group)
  } 
    
    return(distance_df)
}

prepareGroupDistDF <- function(
  group_labels, 
  dist_mat,
  symmetrize = FALSE,
  do_hclust_axes = FALSE,
  preset_levels = NULL,
  do_numeric_x = FALSE,
  do_numeric_y = FALSE) {
  ## HClust axis labels?
  if(do_hclust_axes) {
    distance_df <- hclustBasedAxisLabels(group_labels, dist_mat)
  } else{
    distance_df <- getDistDf(group_labels, dist_mat, 
                             symmetrize = symmetrize,
                             gather = TRUE)
  }
  
  ## Preset levels for Axis?
  if(!is.null(preset_levels)) {
    distance_df <- presetAxisLabels(distance_df, preset_levels)
  } 
  
  ## Append numbers to Axis Labels?
  distance_df <- 
    numericAxisLabels(distance_df, 
                      do_numeric_x, do_numeric_y)
  
  return(distance_df)
}

groupwiseDistanceHeatmap <- function(group_labels, dist_mat, 
                                     title = "",
                                     symmetrize = FALSE,
                                     do_hclust_axes = FALSE,
                                     preset_levels = NULL,
                                     do_numeric_x = FALSE,
                                     do_numeric_y = FALSE,
                                     global_group_label_df = NULL,
                                     text_size = 3.88) {
  
  n_groups <- length(unique(group_labels))
  
  ## Prepare Axis Labels:
  distance_df <- prepareGroupDistDF(
    group_labels, dist_mat,
    symmetrize = symmetrize,
    do_hclust_axes = do_hclust_axes, 
    preset_levels = preset_levels, 
    do_numeric_x = do_numeric_x,
    do_numeric_y = do_numeric_y
  )
  
  ## Create GGPlot
  p <-
    ggplot() + 
    geom_tile(data = distance_df,
              aes(x = target_group, y = source_group, fill = dist)) + 
    scale_fill_gradient2(low = "royalblue4", 
                         high = "orangered2", 
                         mid = "white", 
                         midpoint = 0.5, limit = c(0,1)) + 
    theme_minimal()
  
  if(symmetrize) {
    p <- p + 
      labs(x = "", y = "", fill = "Group-Wise\nDistance")
  } else{
    p <- p + 
      labs(x = "Target", y = "Source", fill = "Group-Wise\nDistance")
  }
    
    
  
  ## Add Top Row Labels, If specified
  if(!is.null(global_group_label_df)) {
    global_group_label_df <- global_group_label_df %>% 
      arrange(cell_group_low) %>% 
      mutate(x = 1:n_groups, y = n_groups+1)
    
    
    
    p <- p + 
      geom_text(data = global_group_label_df,
                aes(x = x, 
                    y = y,
                    label = hm_global_group_label,
                    colour = global_group_label),
                size = text_size) + 
      expand_limits(y = n_groups+2) + 
      theme(panel.grid.major = element_blank()) + 
      guides(
        colour = guide_legend(
          title = "Global Group Label:",
          #title.theme = element_text(size=12),
          override.aes = list(size = 0)))
  }
  
  p <- p + ggtitle(title)
  
  return(p)  
}



getDistDf <- function(
  group_labels, 
  dist_mat,
  symmetrize = FALSE,
  gather = TRUE) {
  n_groups <- length(unique(group_labels))
  
  if(class(group_labels) %in% c("factor", "char", "character")){
    group_labels = as.factor(group_labels)
    group_indices = as.integer(group_labels)
  } else{group_indices = group_labels}
  
  distance_df <-
    data.frame(
      sapply(seq_len(n_groups), 
             function(g) {
               computeGroupDists(g, group_indices, dist_mat)}
             ))
  
  if(class(group_labels) %in% c("factor", "char", "character")){
    names(distance_df) <- levels(group_labels)
  } else{names(distance_df) <- "group" %p% seq_len(n_groups)} 
  
  if(symmetrize == TRUE) {
    symmetric_vals <- 0.5 * (as.matrix(distance_df) + t(as.matrix(distance_df))) %>% 
      data.frame()
    
    distance_df[1:n_groups, 1:n_groups] <-
      symmetric_vals[1:n_groups, 1:n_groups]
  }
  
  distance_df <- distance_df %>%
    mutate(target_group = names(distance_df))
  
  if(gather == TRUE){
    distance_df <- distance_df %>% 
      gather(source_group, dist, setdiff(names(distance_df), c("target_group")))
    
  }
  
  
  return(distance_df)
}





