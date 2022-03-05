###############################################################################
### Computing subgroup similarities and gene contribution scores
###############################################################################

runDimReds <- function(counts, group_labels, 
                       r = 50, g=c(50,150,300), 
                       pcs=c(3, 25, 50),
                       n_clust = 5, n_cores = detectCores()-2) {
  ### 1. SIDEREF --------------------------------------------------------------
  side_ref = list()
  side_ref_dist = list()
  side_ref_umap = list()
  side_ref_plots = list()
  for(i in seq_len(length(g))) {
    cat("sideref rand" %p% i %p% "...")
    
    side_ref[[i]] <- 
      SIDEseqRefSet(expr_matrix = counts, 
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_top_genes = g[i],
                    ## cell reference set selection parameters
                    selection_method = "cell_embed_sample",
                    size_ref_set = r,
                    n_clust = n_clust,
                    B = 0,
                    D = 1,
                    ## parallelization parameters
                    parallelize = TRUE,
                    n_cores = n_cores,
                    ## other
                    verbose = TRUE)
    
    side_ref_dist[[i]] <- side_ref[[i]]$dissim_final
    
    side_ref_umap[[i]] <-
      uwot::umap(as.dist(side_ref_dist[[i]]),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist)
    
    side_ref_plots[[i]] <- 
      splatUMAPPlot(data  = as.dist(side_ref[[i]]$dissim_final),
                    title = "SIDEREF (PCA Embed Samp) " %p% r %p% " Cells " %p% g[i] %p% " Top Genes",
                    celltypes = group_labels,
                    plot_sil_score = FALSE)
  }
  
  ### 2. PCA ------------------------------------------------------------------
  pca_dist  = list()
  pca_umap  = list()
  pca_plots = list()
  
  pca_wtd_dist  = list()
  pca_wtd_umap  = list()
  pca_wtd_plots = list()
  
  pca_res_full <- prcomp(t(counts))
  for(i in seq_len(length(pcs))) {
    cat(pcs[i] %p% " PCs...")
    
    ## principal component stdevs
    pca_sdevs <- pca_res_full$sdev[1:(min(dim(counts)[2], pcs[i]))]
    
    ## Embedding
    pca_res <- pca_res_full$x[, 1:(min(dim(counts)[2], 
                                  pcs[i]))]
    
    pca_wtd_res <- pca_res
    for(c in seq_len(pcs[i])) {
      pca_wtd_res[,c] <- pca_res[,c] * pca_sdevs[c]
    }
    
    pca_res     <- dist(pca_res)
    pca_wtd_res <- dist(pca_wtd_res)
    
    
    mat <- matrix(0, nrow = ncol(counts),
                  ncol = ncol(counts))
    
    ## PCA (regular)
    mat[lower.tri(mat)]  <- c(pca_res)
    
    mat[upper.tri(mat)] <- 
      t(mat)[upper.tri(mat)]
    
    pca_dist[[i]] <- mat
    
    ## PCA (weighted)
    mat[lower.tri(mat)]  <- c(pca_wtd_res)
    
    mat[upper.tri(mat)] <- 
      t(mat)[upper.tri(mat)]
    
    pca_wtd_dist[[i]] <- mat
    
    
    
    pca_umap[[i]] <-
      uwot::umap(as.dist(pca_dist[[i]]),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist) 
    
    pca_wtd_umap[[i]] <-
      uwot::umap(as.dist(pca_wtd_dist[[i]]),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist) 
    
    pca_plots[[i]] <-
      splatUMAPPlot(as.dist(pca_dist[[i]]), 
                    title = "Euclidean Applied to " %p% pcs[i] %p% " PCs", 
                    celltypes = group_labels,
                    plot_sil_score = FALSE)
    
    pca_wtd_plots[[i]] <-
      splatUMAPPlot(as.dist(pca_wtd_dist[[i]]), 
                    title = "Euclidean Applied to " %p% pcs[i] %p% "Wtd. PCs", 
                    celltypes = group_labels,
                    plot_sil_score = FALSE)
    
    
  }
  

  ### 3. Euclid UMAP ------------------------------------------------------------
  
  cat("euclidean and correlation plots...")
  
  euclid_dist <- dist(t(counts))
  
  mat <- matrix(0, nrow = ncol(counts),
                ncol = ncol(counts))
  
  mat[lower.tri(mat)]  <- c(euclid_dist)
  
  mat[upper.tri(mat)] <- 
    t(mat)[upper.tri(mat)]
  
  euclid_dist <- mat
  
  
  euclid_umap <-
    uwot::umap(as.dist(euclid_dist),
               n_neighbors = n_neighbors,
               min_dist = min_dist) 
  
  p_splat_euclid <-
    splatUMAPPlot(as.dist(euclid_dist), 
                  title = "Euclidean Distance", 
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  
  ### 4. Pearson UMAP -----------------------------------------------------------
  
  pear_dist <- 1 - abs(cor(counts))
  
  pear_umap <-
    uwot::umap(as.dist(pear_dist),
               n_neighbors = n_neighbors,
               min_dist = min_dist)
  
  p_splat_pear <-
    splatUMAPPlot(as.dist(pear_dist), title = "Pearson Distance", 
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  
  ### 5. Spearman UMAP ----------------------------------------------------------
  
  spearman_dist <- 1 - abs(cor(counts, method = "spearman"))
  
  spearman_umap <-
    uwot::umap(as.dist(spearman_dist),
               n_neighbors = n_neighbors,
               min_dist = min_dist)
  
  p_splat_spearman <-
    splatUMAPPlot(as.dist(spearman_dist), title = "Spearman Distance", 
                  celltypes = group_labels,
                  plot_sil_score = FALSE)
  
  

  names_core <- 
    c("side_ref_g" %p% g,
      "pca_" %p% pcs,
      "pca_wtd" %p% pcs,
      "euclid", "pear", "spearman")
  
  names_dist <- names_core %p% "_dist"
  names_umap <- names_core %p% "_umap"
  names_plot <- names_core %p% "_plot_list"
  
  output <-
    list(
      dist_list = append(side_ref_dist, 
                         append(pca_dist,
                                append(pca_wtd_dist,
                                       list(euclid_dist,
                                            pear_dist,
                                            spearman_dist)))),
      umap_list = append(side_ref_umap,
                         append(pca_umap,
                                append(pca_wtd_umap,
                                       list(euclid_umap,
                                            pear_umap,
                                            spearman_umap)))),
      
      plot_list = append(side_ref_plots,
                         append(pca_plots,
                                append(pca_wtd_plots,
                                       list(p_splat_euclid, 
                                            p_splat_pear, 
                                            p_splat_spearman))))
    )
  names(output$dist_list) = names_dist 
  names(output$umap_list) = names_umap 
  names(output$plot_list) = names_plot
  
  return(output)
}




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



groupwiseDistanceHeatmap <- function(group_labels, dist_mat, 
                                     title = "",
                                     hclust = FALSE,
                                     preset_levels = NULL,
                                     numeric_x=FALSE,
                                     numeric_y=FALSE) {
  n_groups <- length(unique(group_labels))
  
  
  if(hclust == TRUE){
    distance_df <- 
      getDistDf(group_labels, 
                dist_mat,
                gather = FALSE) %>% 
      dplyr::select(-group)
    
    distance_df_sym <- 
      distance_df
    
    distance_df_sym <- 0.5*(distance_df_sym + t(distance_df_sym))
    
    cluster     = hclust(dist(distance_df_sym), method="ward.D2")
    dendrogram  = as.dendrogram(cluster)
    Rowv        = rowMeans(distance_df, na.rm = T)
    dendrogram  = reorder(dendrogram, -Rowv)
    
    
    rowInd <- rev(order.dendrogram(dendrogram))
    colInd <- rowInd
    distance_df <- distance_df[rowInd, colInd]
    
    ## change levels
    distance_df <- distance_df %>% mutate(group = factor(names(distance_df)))
    distance_df$group <- forcats::fct_inorder(distance_df$group)
    
    distance_df <- distance_df %>%
      gather(group2, dist, setdiff(names(distance_df), c("group")))
    
    ## change levels (group2)
    distance_df$group2 = factor(distance_df$group2, levels = levels(distance_df$group))
    
  } else{
    distance_df <-
      getDistDf(group_labels, dist_mat, gather = TRUE)}
  
  
  if(!is.null(preset_levels)){
    ## TODO: error checking on preset_levels entry?!
    distance_df$group = factor(distance_df$group, levels = preset_levels)
    distance_df$group2 = factor(distance_df$group2, levels = preset_levels)
    
    ## remove groups not in the preset levels.
    distance_df <- distance_df %>% 
      dplyr::filter(!is.na(group) & !is.na(group2))
  }
  
  if(numeric_x) {
    
    distance_df$group = 
      factor(as.character(as.numeric(distance_df$group)), 
             levels = as.character(seq_len(max(
               as.numeric(distance_df$group)))))
    
    
    distance_df <- distance_df %>% 
      dplyr::mutate(new_group = factor(
        group2 %p% ", " %p% as.numeric(group2),
        levels = levels(group2) %p% ", " %p% 
          seq_len(length(levels(group2)))))
    if(numeric_y) {
      distance_df <- distance_df %>% 
        mutate(new_group = 
                 factor(as.character(as.numeric(new_group)),
                        levels = as.character(seq_len(
                          max(as.numeric(new_group))))
                 )
        )
    }  
    distance_df <- distance_df %>% 
      mutate(group2=new_group)
  } 
  
  p <-
    distance_df %>% 
    ggplot(aes(x=group, y=group2, fill=dist)) + 
    geom_tile() + 
    scale_fill_gradient2(low = "royalblue4", 
                         high = "orangered2", 
                         mid = "white", 
                         midpoint = 0.5, limit = c(0,1)) + 
    labs(x = "Target", y = "Source", fill = "Group-Wise\nDistance") + 
    ggtitle(title)
  
  return(p)  
}



getDistDf <- function(group_labels, dist_mat, gather = TRUE) {
  n_groups <- length(unique(group_labels))
  
  if(class(group_labels) %in% c("factor", "char", "character")){
    group_labels = as.factor(group_labels)
    group_indices = as.integer(group_labels)
  } else{group_indices = group_labels}
  
  distance_df <-
    data.frame(sapply(seq_len(n_groups), 
                      function(g){computeGroupDists(g, group_indices, dist_mat)}))
  
  if(class(group_labels) %in% c("factor", "char", "character")){
    names(distance_df) <- levels(group_labels)
  } else{names(distance_df) <- "group" %p% seq_len(n_groups)} 
  
  distance_df <- distance_df %>%
    mutate(group = names(distance_df))
  
  if(gather == TRUE){
    distance_df <- distance_df %>% 
      gather(group2, dist, setdiff(names(distance_df), c("group")))
    
  }
  
  
  return(distance_df)
}



### TODO: should update functions to enforce proper ordering of labels in DF
### TODO: the distance DF should have from_group and to_group as the names, not group1, group2.

computeRankAccuracy <- function(group_labels, 
                                dist_mat,
                                ground_truth,
                                dist_name = NULL) {
  
  n_groups <- length(unique(group_labels))
  
  distance_df <- getDistDf(group_labels, dist_mat) %>%
    mutate(ground_truth_dist = ground_truth) %>% 
    dplyr::group_by(group2) %>% 
    mutate(rank = rank(dist),
           ground_truth_rank = rank(ground_truth_dist))
  
  
  overall_rank_acc <- mean(distance_df$rank==distance_df$ground_truth_rank)
  
  overall_mean_rank_error <- mean(abs(distance_df$rank - distance_df$ground_truth_rank))
  
  
  groupwise_rank_stats <- 
    distance_df %>% 
    group_by(group2) %>% 
    dplyr::summarise(spearman = 1 - 6*sum((rank-ground_truth_rank)^2)/(n_groups^3-n_groups),
                     rank_acc = mean(rank == ground_truth_rank),
                     mean_rank_error = mean(abs(rank - ground_truth_rank)))
  
  
  
  return(list(overall_rank_acc          = overall_rank_acc,
              overall_mean_rank_error   = overall_mean_rank_error,
              overall_mean_spearman     = groupwise_rank_stats %>% dplyr::pull(spearman) %>% mean(),
              groupwise_rank_acc        = groupwise_rank_stats %>% dplyr::pull(rank_acc),
              groupwise_mean_rank_error = groupwise_rank_stats %>% dplyr::pull(mean_rank_error),
              groupwise_spearman        = groupwise_rank_stats %>% dplyr::pull(spearman),
              group_name_vec = groupwise_rank_stats$group2,
              dist_name = dist_name))
  
}



### Compute MSEs
computeMSEGroupDist <- function(group_labels, 
                                dist_mat,
                                ground_truth) {
  
  n_groups <- length(unique(group_labels))
  
  distance_df <- getDistDf(group_labels, dist_mat) 
  
  full_mse = mean((ground_truth - distance_df$dist)^2) 
  
  group_wise_mse = rep(0, n_groups)
  for(i in seq_len(n_groups)) {
    group_wise_mse[i] = mean((ground_truth[((i-1)*n_groups+1):(n_groups*i)] - 
                                distance_df %>% 
                                dplyr::filter(group2 == "group" %p% i) %>% 
                                pull(dist))^2) 
  }
  
  return(list(full_mse = full_mse, group_wise_mse = group_wise_mse))
  
}


subgroupSilWidth <- function(group_indices, subgroup_indices, dist_mat,
                             dist_name = NULL, ind_stat = "mean", group_stat = "mean") {
  n_cells = length(group_indices)
  
  sil_widths = rep(NA, n_cells)
  
  for(c in seq_len(n_cells)) {
    
    subgroup = subgroup_indices[c]
    
    ## First get full subgroup
    subgroup_cells = which(subgroup_indices == subgroup)
    
    ## Other group of cells
    other_cells = setdiff(seq_len(n_cells), subgroup_cells)
    
    ## Remove cell's own group from the subgroup
    subgroup_cells = setdiff(subgroup_cells, which(group_indices==group_indices[c]))
    
    if(length(subgroup_cells)==0){next}
    
    other_groups = unique(group_indices[other_cells])
    
    a = ifelse(ind_stat == "mean", mean(dist_mat[c, subgroup_cells]),
               median(dist_mat[c, subgroup_cells]))
    
    b = sapply(other_groups, function(x) {
      b_cells = which(group_indices == x)
      return(ifelse(ind_stat == "mean", 
                    mean(dist_mat[c, b_cells]),
                    median(dist_mat[c, b_cells]))
      )}) %>% 
      min()
    
    sil_widths[c] = (b-a) / max(a,b)
    
  }
  
  n_groups = length(unique(group_indices))
  
  group_sil_widths = rep(0, n_groups)
  
  group_sil_widths = sapply(seq_len(n_groups), 
                            function(x) ifelse(group_stat == "mean",
                                               mean(sil_widths[which(group_indices==x)],
                                                    na.rm = TRUE),
                                               median(sil_widths[which(group_indices==x)],
                                                      na.rm = TRUE)))
  avg_sil_width = ifelse(group_stat == "mean", 
                         mean(sil_widths, na.rm = TRUE),
                         median(sil_widths, na.rm = TRUE))
  
  
  sil_widths_df <-
    data.frame(group = c("all", "group" %p% seq_len(n_groups)),
               sil_width = c(avg_sil_width, group_sil_widths))
  
  if(!is.null(dist_name)) {
    sil_widths_df <- sil_widths_df %>% mutate(method = dist_name)
  }
  
  return(sil_widths_df)
  
}


