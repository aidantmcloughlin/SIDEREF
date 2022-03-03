library(here)
source(here('R/libraries.R'))
source(here('R/subgroup_composition_gene_contrib.R'))

library(network)
library(GGally)
library(ggnewscale)
library(gridExtra)
library(gtable)
library(grid)


## TODO: some plot element location args just be adjustable to the user.


getBipartiteSegments <- function(distance_df,
                                 source_groups) {
  
  n_groups <- length(unique(distance_df$target_group))
  
  ## build segment x-axis
  x2 <- seq_len(n_groups)
  
  x1_inds <- which(unique(distance_df$target_group) %in% source_groups)
  x1 <- x2[x1_inds]
  
  ## build segment y-axis
  y1 <- rep(1, length(x1))
  y2 <- rep(0, length(x2))
  
  ## construct DF
  edges <- expand.grid(x1 = x1, x2 = x2) %>% 
    cbind(expand.grid(y1 = y1, y2 = y2))
  
  ## add distances 
  
  edges <- edges %>% 
    arrange(x1) %>% 
    mutate(
      dist = 
        distance_df[distance_df$source_group %in% source_groups, ]$dist)
    
  
  return(edges)
  
} 

bipartiteNetworkGraph <- 
  function(dist_mat,
           group_labels,
           source_groups,
           symmetrize = FALSE,
           do_hclust_axes = FALSE,
           preset_levels = NULL,
           title = NULL,
           legend_title = "Cluster Name",
           node_size = 5,
           text_size = 8,
           node_fill = "grey75",
           node_color = "black") {
    
    
    n_groups <- length(unique(group_labels))
    
    ### Relative Distances DF
    distance_df <- prepareGroupDistDF(
      group_labels, dist_mat,
      symmetrize = symmetrize,
      do_hclust_axes = do_hclust_axes, 
      preset_levels = preset_levels, 
      do_numeric_x = F,
      do_numeric_y = F
    )
    
    distance_df <- distance_df %>% 
      arrange(source_group, target_group)
    
    ### DF of line segments
    edges <- getBipartiteSegments(
      distance_df,
      source_groups = source_groups)
    
    p <- 
      ggplot() + 
      geom_segment(
        data = edges, 
        aes(x = x1, y = y1, 
            xend = x2, yend = y2,
            colour = dist)) +
      # scale_colour_gradient2(low = "royalblue4",
      #                       high = "orangered2",
      #                       mid = "white",
      #                       midpoint = 0.5, limit = c(0,1),
      #                       guide = guide_colourbar(
      #                         title = 'Group-Wise\nDistance',
      #                         title.theme = element_text(size=8.5),
      #                         direction = "horizontal"))
      scale_colour_gradient(low = "royalblue4",
                             high = "white",
                             #midpoint = 0.5, 
                            limit = c(0,1),
                             guide = guide_colourbar(
                               title = 'Group-Wise\nDistance',
                               title.theme = element_text(size=8.5),
                               direction = "horizontal"))
    
    p <- p +
      ggnewscale::new_scale_colour() +
      geom_point(data = edges %>% distinct(x1, y1), 
                 aes(x=x1, y=y1),
                 size = node_size,
                 shape = 21,
                 colour = node_color,
                 fill = node_fill) +
      geom_point(data = edges %>% distinct(x2, y2), 
                 aes(x=x2, y=y2),
                 size = node_size,
                 shape = 21,
                 colour = node_color,
                 fill = node_fill) +
      theme_void() +
      expand_limits(y = -0.25) + expand_limits(y = 1.25) +
      expand_limits(x = -1) + expand_limits(x = n_groups + 0.25)
    
    ## Labels
    keys <- 
      unique(edges$x2) %p% ", " %p% 
      unique(distance_df$target_group)
    
    p <- p + 
      geom_text(data = data.frame(x=unique(edges$x1),
                                  y=unique(edges$y1),
                                  label=unique(edges$x1)),
                aes(x=x,y=y,label=label),
                size = text_size) + 
      geom_text(
        data = data.frame(
          x=unique(edges$x2),
          y=unique(edges$y2),
          label=unique(edges$x2)) %>% 
          mutate(key = factor(keys, levels = keys)),
        aes(x=x,y=y,label=label, 
            colour=key),
        size = text_size) + 
      scale_colour_manual(values=rep("black", n_groups))
    
    ## Add "Source" and "Target" monikers
    
    ## Old location of labels:
    # med_x <- median(seq_len(n_groups))
    # p <- p + 
    #   geom_text(data = data.frame(x=med_x, y=1.2,label="Source"),
    #              aes(x=x,y=y,label=label)) +
    #   geom_text(data = data.frame(x=med_x, y=-0.2,label="Target"),
    #             aes(x=x,y=y,label=label))
    
    ## Vertically aligned with nodes:
    p <- p + 
      geom_text(data = data.frame(x=-0.5, y=1,label="Source:"),
                aes(x=x,y=y,label=label)) +
      geom_text(data = data.frame(x=-0.5, y=0,label="Target:"),
                aes(x=x,y=y,label=label))
    
    ## null points with legend
    p <- 
      p + 
      guides(
        colour = guide_legend(
          title = legend_title %p% ":",
          title.theme = element_text(size=8.5),
          override.aes = list(size = 0)))
    
    
    if(!is.null(title)) {
      p <- p + ggtitle(title)
    }
    
    p
    
}






