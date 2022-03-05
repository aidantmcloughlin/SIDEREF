spectralEmbed <- function(dissim,
                          is_dissim = TRUE,
                          ndim = 25) {
  if(is_dissim == TRUE) {
    ## convert to similarity (adjacency)
    dissim_norm = dissim / max(dissim)
    sim = 1 - dissim_norm
    diag(sim) = 0
  } else{sim = dissim}
  
  ## Spectral Embedding Procedure:
  n = ncol(sim)
  S = rowSums(sim)
  #D = diag(S)
  D_sqrt_inv = diag(1/S^0.5)
  #L = D - sim
  ## normalized graph laplacian
  L = diag(1, n, n) - D_sqrt_inv %*% sim %*% D_sqrt_inv
  
  
  cols_extract <- (n-1):(n-ndim)
  ## collect top eigenvalues
  eig_weights <- 1/eigen(L)$values[cols_extract]
  
  U = (eigen(L)$vectors)[, cols_extract]
  row_norms = apply(U, 1, function(x) sqrt(sum(x^2)))
  for(i in 1:n) {U[i,] = U[i,] / row_norms[i]}
  
  return(list(U=U, eig_weights=eig_weights))
  
}


addSpectralDistRes <- function(dist_res, 
                               spectral_dims = c(10),
                               do_spectral_pca=TRUE) {
  
  print(spectral_dims %p% " spectral dims")
  
  if(do_spectral_pca) {
    dist_names <- names(dist_res$dist_list)
  } else{
    dist_names <- names(dist_res$dist_list[
      !grepl("pca", names(dist_res$dist_list))])
  }
  
  for(dist_name in dist_names) {
    for(s in seq_len(length(spectral_dims))) {
      n <- ncol(dist_res$dist_list[[1]])
      
      l <- length(dist_res$dist_list)
      
      dissim <- dist_res$dist_list[[dist_name]]
      
      spec_embed_res <- spectralEmbed(dissim, ndim = spectral_dims[s]) 
      
      ## standard:
      dissim <- dist(spec_embed_res$U)
      
      mat <- matrix(0, nrow = n, ncol = n)
      
      mat[lower.tri(mat)]  <- c(dissim)
      
      mat[upper.tri(mat)] <-
        t(mat)[upper.tri(mat)]
      
      dist_res$dist_list[[l+1]] <- mat
      names(dist_res$dist_list)[l+1] <- 
        dist_name %p% "_spectral_d" %p% spectral_dims[s] %p% "_dist"
      
      dist_res$umap_list[[l+1]] <- 
        uwot::umap(as.dist(dist_res$dist_list[[l+1]]),
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
      
      names(dist_res$umap_list)[l+1] <- names(dist_res$dist_list)[l+1]
      ## weighted:
      U_wtd = spec_embed_res$U
      for(c in seq_len(ncol(spec_embed_res$U))) {
        U_wtd[,c] <- U_wtd[,c] * spec_embed_res$eig_weights[c]
      }
      
      dissim <- dist(U_wtd)
      
      mat <- matrix(0, nrow = n, ncol = n)
      
      mat[lower.tri(mat)]  <- c(dissim)
      
      mat[upper.tri(mat)] <-
        t(mat)[upper.tri(mat)]
      
      dist_res$dist_list[[l+2]] <- mat
      names(dist_res$dist_list)[l+2] <- 
        dist_name %p% "_spectral_wtd_d" %p% spectral_dims[s] %p% "_dist"
      
      dist_res$umap_list[[l+2]] <- 
        uwot::umap(as.dist(dist_res$dist_list[[l+2]]),
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
      
      names(dist_res$umap_list)[l+2] <- names(dist_res$dist_list)[l+2] 
    }
  }
  return(dist_res)
}
