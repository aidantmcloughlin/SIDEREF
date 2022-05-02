SIDEREF <- function(expr_matrix,
                    n_top_genes = 300,
                    R = 100,
                    ## cell reference set selection parameters
                    selection_method = "cell_embed_sample",
                    diff_expr_method = "diff_expr_norm",
                    similarity_method = "n_intersect",
                    n_pcs = NULL, 
                    n_clust = NULL,
                    D = 1,
                    return_sim = FALSE, 
                    n_neighbors = 15, 
                    min_dist = 0.01,
                    pca_elbow_chg_thres = 0.025,
                    kmeans_elbow_chg_thres = 0.025,
                    ## parallelization parameters
                    parallelize = TRUE,
                    n_cores = detectCores()-1,
                    ## other
                    verbose = TRUE) {
  #' Computes SIDEREF distance matrix between cells in scRNA seq array.  
  #' SIDEREF uses concordance of ranked differential expression gene lists 
  #' to measure dissimilarity between cells. A globally representative reference 
  #' set of cells forms the basis for the ranked DE gene lists.
  #' @param expr_matrix Gene x Cell matrix of scRNA seq data.
  #' @param n_top_genes Number of genes comprising each DE gene list.
  #' @param R Size of cell reference set.
  #' @param selection_method Whether to use random sampling or stratified sampling of k-means generated clusters to form the reference set.
  #' @param n_pcs Number of PCs to use in PCA embedding for cluster stratified sampling of the reference set.
  #' @param n_clust Number of clusters to use in K-means for cluster stratified sampling of the reference set.
  #' @param D Number of runs of SIDEREF to average results over.
  #' @param return_sim If true, returns the similarity matrix.
  #' @param n_neighbors Number of neighbors parameter to use in UMAP embedding for cluster stratified sampling of the reference set.
  #' @param min_dist Min distance parameter to use in UMAP embedding for cluster stratified sampling of the reference set.
  #' @param pca_elbow_chg_thres Elbow criteria variance change threshold to use for adaptively selecting PCs
  #' @param kmeans_elbow_chg_thres Elbow criteria variance change threshold to use for adaptively selecting number of clusters.
  #' @param parallelize Whether to parallelize computations.
  #' @param n_cores How many cores to use in the parallel computation.
  #' @param verbose Whether to print status messages during SIDEREF computation.
  #' @return dissim_final: a Cell x Cell SIDEREF dissimilarity (similarity, if specified) matrix.
  
  
  ## initializations:
  dissim_final <- NULL
  
  ## set constants, parallel backend
  J <- dim(expr_matrix)[1]
  N <- dim(expr_matrix)[2]
  
  similarity_idx_method <- "SIDEREF"
  cell_pairs <- combn(N, 2)
  
  if(parallelize == TRUE) {
    cl <- makeCluster(n_cores) 
    if(verbose) {cat("Cores in usage: " %p% n_cores)}
  }
  
  
  ## Set Differential Expression measure
  if(diff_expr_method == "diff_expr_norm") {
    diffExprMeasure <- function(i, j, expr_matrix) {
      return(abs(expr_matrix[, i] - expr_matrix[, j]) / 
               sqrt(expr_matrix[, i] + expr_matrix[, j]))
    }
  }

  ## Set Similarity measure
  similarityMeasure <- 
    setSimilarityMeasure(similarity_method = similarity_method, 
                         index_method = similarity_idx_method)
  
  ### Dissimilarity Matrix Computation Workflow:
  computeDissim <- function(return_sim = FALSE) {
    
    if(class(ref_cells) == "logical") {ref_cells <- which(ref_cells)}
    
    ## create ref set combn list
    ref_set_pairs <-
      matrix(c(rep(seq_len(N), length(ref_cells)),
               sapply(ref_cells, rep, N)), 
             byrow = TRUE,
             nrow  = 2)
    
    ## Compute Differential Expression Vector for each pair of cells 
    diff_expr_matrix <-
      computeDiffExprMat(cell_pairs = ref_set_pairs,
                         expr_matrix = expr_matrix, 
                         diff_expr_measure = diffExprMeasure,
                         n_top_genes = n_top_genes,
                         parallelize = parallelize,
                         cl = cl) 
    
    
    ## Compute similarity scores
    if(parallelize == TRUE) {
      clusterExport(cl, varlist=c("diff_expr_matrix", "N",
                                  "ref_cells",
                                  "similarityMeasure"),
                    envir=environment())
      similarity_vec <-
        parApply(cl = cl, X = cell_pairs, 2, 
                 FUN = function(x) {
                   return(
                     similarityMeasure(i = x[1], j = x[2], N = N,
                                       ref_set = ref_cells,
                                       d_e_mat = diff_expr_matrix)
                   )})
    } else{
      similarity_vec <-
        apply(cell_pairs, 2, 
              function(x) {
                similarityMeasure(i = x[1], j = x[2], N = N, 
                                  ref_set = ref_cells,
                                  d_e_mat = diff_expr_matrix)})
      
    }
    
    ## dissimilarity measure is max score - similarity_score
    if(return_sim == FALSE) {
      mat_fill_vec <- max(similarity_vec) - similarity_vec  
    } else{mat_fill_vec = similarity_vec}
    
    
    ## convert vector to cell-by-cell matrix
    fill_matrix <- matrix(0, nrow = N, ncol = N)
    
    fill_matrix[lower.tri(fill_matrix)] <- mat_fill_vec
    
    
    fill_matrix[upper.tri(fill_matrix)] <- 
      t(fill_matrix)[upper.tri(fill_matrix)]
    
    return(fill_matrix)
  }

  
  ## Get Reference Cells for each averaging iteration
  ref_cell_list <- 
    lapply(seq_len(D),
           function(d) {
             return(selectRefSet(expr_matrix,
                                 selection_method = selection_method,
                                 R = R,
                                 dissim_matrix = NULL,
                                 n_pcs = n_pcs,
                                 n_clust = n_clust,
                                 n_neighbors = n_neighbors,
                                 min_dist = min_dist,
                                 pca_elbow_chg_thres = pca_elbow_chg_thres,
                                 kmeans_elbow_chg_thres = kmeans_elbow_chg_thres
                                 ))
           })
  
  
  ## Run averaging loop
  for(d in seq_len(D)) {
    if(verbose == TRUE) {cat("Running averaging iteration", d, "of", D, "...")}
    
    ref_cells <- ref_cell_list[[d]]
    dissim_curr <- computeDissim(return_sim = return_sim)
    
    if(is.null(dissim_final)){dissim_final <- dissim_curr}
    ## incorporate new matrix into running average.
    dissim_final <-
      1/(d+1) * (d * dissim_final + dissim_curr)
    
  }
  
  
  ## close workers if needed
  if(parallelize == TRUE) {
    stopCluster(cl)
  }
  
  return(dissim_final)
  
}



sampleClusters <- function(clust_res, samp_size, N) {
  ## error checking
  if(samp_size >= N) {stop("Sample size is larger than number of cells")}
  
  ## sample an even proportion of clusters
  fracs <- table(clust_res) / N
  samps <- round(fracs * samp_size)
  
  samples <-
    unlist(
      sapply(seq_len(max(clust_res)), 
             function(c) {
               ## random stratified sampling from each cluster.
               in_clust <- which(clust_res == c)
               return(sample(in_clust, samps[c], replace = FALSE))
             }))
  
  is_sampled <- sapply(seq_len(N), function(x) x %in% samples)
  
  return(is_sampled)
}



selectRefSet <- function(expr_matrix,
                         selection_method = "random",
                         R = 100,
                         dissim_matrix = NULL,
                         n_pcs = NULL,
                         n_clust = NULL,
                         ## UMAP params:
                         n_neighbors = 15, 
                         min_dist = 0.01,
                         ## Other:
                         max_clust_try = 25,
                         max_pcs_store = 100,
                         pca_elbow_chg_thres = 0.025,
                         kmeans_elbow_chg_thres = 0.025
) {
  
  N <- dim(expr_matrix)[2]
  
  cells <- NULL
  
  if(R >= N) {
    cells <- seq_len(N)
  }
  else if(selection_method == "random") {
    cells <- sample(N, R, replace = FALSE)
  }
  else if(selection_method == "cell_embed_sample") {
    ## UMAP + KMeans.
    ## first get num PCs for UMAP
    ## adaptively select number of PCs, if needed:
    if(is.null(n_pcs)) {
      pca_res <- irlba::prcomp_irlba(expr_matrix, n=max_pcs_store, 
                                     retx = TRUE, 
                                     center = TRUE,
                                     scale = FALSE)
      
      n_pcs <- quick_elbow(pca_res$sdev**2, 
                           low = pca_elbow_chg_thres, 1)
      
    }
    if(is.null(dissim_matrix)) {
      
      umap_embed <- 
        uwot::umap(t(expr_matrix),
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   pca = n_pcs)
    } else{
      umap_embed <-
        uwot::umap(stats::as.dist(dissim_matrix),
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
    }
    
    ## adaptively select number of clusters, if needed:
    if(is.null(n_clust)) {
      wss <- rep(0, max_clust_try-1)
      for(k in seq(2, max_clust_try)) {
        wss[k-1] <- kmeanspp(umap_embed, k = k)$tot.withinss
      }
      
      ## quick_elbow from WSS to decide number of clusters
      n_clust <- quick_elbow(
        wss, 
        low = kmeans_elbow_chg_thres, 
        1) + 1
      
    }
    
    ## now with clusters decided, run k-means again
    clust_res <- kmeanspp(umap_embed, k = n_clust)$cluster
    
    
    ## sampling from cluster results
    cells <- sampleClusters(clust_res, 
                            samp_size = R,
                            N = N)
  }
  
  
  ## choose cell ref set from ranks if needed
  if(is.null(cells)) {
    cells <- which(ranks <= R)
  }
  
  return(cells)
  
}
