


predictKNN <- function(knn_list, data_set) {
  ## Args:
  ##    knn_list: list of KNN for each observation
  ##    x: training data for each observation.
  knn_preds <-
      do.call(rbind,
              lapply(seq_len(N),
                     function(i) {
                       return(colMeans(data_set[knn_list[[i]], ]))
                     }))
  
  return(knn_preds)
}


weightedNearestNeighbors <- function(dist_mats, data_sets,
                                     k = 20) {
  
  N <- nrow(x)
  M <- length(dist_mats)
  
  theta <- list()
  exp_s <- list()
  weights <- list()
  
  
  for(m1 in seq_len(M)) {
    theta[[m1]] <- list()
    
    ## First get cell kernel bandwidths
    jaccard_vec <-
      apply(combn(nrow(data_sets[[m1]]), 2), 2, function(x) {
        a <- knn_list[[x[1]]]
        b <- knn_list[[x[2]]]
        length(intersect(a, b)) / 
          length(union(a, b))
      })
    
    jaccard_mat <- matrix(0, nrow = N, ncol = N)
    
    jaccard_mat[lower.tri(jaccard_mat)] <- jaccard_vec
    
    
    jaccard_mat[upper.tri(jaccard_mat)] <- 
      t(jaccard_mat)[upper.tri(jaccard_mat)]
    
    
    low_jaccard_sim <- 
      lapply(seq_len(N),
             function(i) {
               low_sim_cells <- 
                 setdiff(order(jaccard_mat[,i])[-1],
                         which(jaccard_mat[,i] < 1e-8))[1:k]
               return(low_sim_cells)
             })
    
    cell_bwidth <- 
      sapply(seq_len(N),
             function(i) {
               mean(sqrt(rowSums( (data_sets[[m1]][i, ] - 
                                     data_sets[[m1]][low_jaccard_sim[[i]], ])^2 )))
             })
    
    dist_nn <- 
      sapply(seq_len(N),
             function(i) {
               nn <- knn_list[[i]][1]
               return(sqrt( sum( (data_sets[[m1]][i, ] - data_sets[[m1]][nn, ])^2 )))
             })
    
    ## Compute KNN prediction distance from actual
    for(m2 in seq_len(M)) {
      knn_list <- 
        lapply(seq_len(dim(dist_mats[[m2]])[1]),
               function(i) {
                 return(order(dist_mats[[m2]][,i])[-1][1:k])
               })
      
      
      
      knn_preds <- predictKNN(knn_list, data_sets[[m1]])
      
      
      
      pred_act_dist <- 
        sapply(seq_len(N),
               function(i) {
                 return(sqrt( sum( (data_sets[[m1]][i, ] - knn_preds[i, ])^2 ) ))
               }) 
      
      
      
      
      theta[[m1]][[m2]] <- exp(-pmax(pred_act_dist - dist_nn, 0) / 
                                 (cell_bwidth - dist_nn)
      )
    }
    
    
  }
  
  
  ## compute ratio of within modality affinity to across modality affinity
  for(m in seq_len(M)) {
    others_mods   <- setdiff(seq_len(M), m)   
    exp_s[[m]] <- exp(theta[[m]][[m]]  / (mean(theta[[m]][[others_mods]]) + 1e-4))
  }
  
  ## normalize ratios under softmax
  for(m in seq_len(M)) {
    weights[[m]] <- exp_s[[m]] / Reduce("+", exp_s)
  }
  
  return(weights)
  
}




