
selectRefSet <- function(expr_matrix,
                         selection_method = "random",
                         size_ref_set = 500,
                         gamma = NULL,
                         sim_matrix = NULL,
                         n_clust = 10) {
  
  ## sample_weight_decay: if provided, cells are sampled with weights: weight_i = gamma^((rank_i - 1)/N)
  N <- dim(expr_matrix)[2]
  
  cells <- NULL
  
  if(size_ref_set >= N) {
    cells <- seq_len(N)
  }
  else if(selection_method == "random") {
    cells <- sample(N, size_ref_set, replace = FALSE)
  }
  else if(selection_method == "total_reads") {
    ranks <- frank(-colSums(expr_matrix), 
                   ties.method = "random")
  }
  else if(selection_method == "n_genes_detected") {
    ranks <- frank(-colSums(expr_matrix > 0),
                   ties.method = "random")
  }
  else if(selection_method == "n_genes_low_expr") {
    q <- quantile(expr_matrix[expr_matrix > 0], 0.1) + 1e-7 #comparison stability
    ranks <- frank(-colSums(expr_matrix <= q),
                   ties.method = "random")
  }
  else if(selection_method == "n_highly_variable_genes") {
    high_var_genes <- getHighlyVariableGenes(expr_matrix)
    ranks <- frank(-colSums(expr_matrix[high_var_genes,] > 0),
                   ties.method = "random")
  }
  else if(selection_method == "cell_embed_sample") {
    ## compute similarity scores if needed
    if(is.null(sim_matrix)) {
      ## Spearman correlation
      sim_matrix <- rcorr(as.matrix(expr_matrix),
                          type = "spearman")$r
      diag(sim_matrix) <- 0
    }
    
    ## spectral clustering
    spect_clust_res <- spectralCluster(sim_matrix, 
                                       n_clust = n_clust)
    
    ## sampling from cluster results
    cells <- sampleClusters(spect_clust_res, 
                            samp_size = size_ref_set)
  }
  
  
  ## choose cell ref set from ranks if needed
  if(is.null(cells)) {
    ## top s ranks if not sampling weight provided
    if(is.null(gamma)){
      cells <- which(ranks <= size_ref_set)
    } else{
      ## check gamma in (0,1)
      if(gamma < 0 | gamma > 1) stop("Please provide gamma in range (0,1).")
      
      ## compute rank-based sampling weights and sample.
      weights <- gamma ^ ((ranks - 1)/(N-1))
      cells <- sample(N, size_ref_set, prob = weights,
                      replace = FALSE)
    }
  }
  
  return(cells)
  
}

SIDEseqRefSet <- function(expr_matrix, 
                          diff_expr_method = "diff_expr_norm",
                          similarity_method = "top_n_gene_intersect",
                          n_top_genes = 300,
                          ## cell reference set selection parameters
                          selection_method = "random",
                          gamma = NULL,
                          size_ref_set,
                          n_clust = 10,
                          B = 5,
                          prior_mat_method = "new", ## c("new", "average")
                          ## parallelization parameters
                          parallelize = TRUE,
                          n_cores = detectCores()-1,
                          ## other
                          verbose = TRUE) {
  
  ## set constants, parallel backend
  J <- dim(expr_matrix)[1]
  N <- dim(expr_matrix)[2]
  
  SIDESeq_method <- "ref_set"
  
  cell_pairs <- combn(N, 2)
  
  if(parallelize == TRUE) {
    cl <- makeCluster(n_cores) 
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
                         side_seq_method = SIDESeq_method)
  
  
  
  ## Iterate procedure B times or until convergence criteria is reached
  b <- 1
  prior_mat <- NULL
  sim_matrix <- NULL
  
  ## TODO: populate this:
  convergence_res <- data.frame(b = seq_len(B),
                                frob_norm_dist         = NA,
                                spectral_embed_corr    = NA,
                                avg_spectral_sil_score = NA)
  ## TODO: convergence criterion
  while(b <= B) {
    if(verbose == TRUE) {cat("Running iteration", b, "of", B, "...")}
    ## extract reference set
    ref_cells <-
      selectRefSet(expr_matrix,
                   selection_method = selection_method,
                   size_ref_set = size_ref_set,
                   gamma = gamma,
                   sim_matrix = sim_matrix,
                   n_clust)
    
    ## TODO: clean this up..
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
    
    ## convert vector to cell-by-cell matrix
    sim_matrix <- matrix(0, nrow = N, ncol = N)
    
    sim_matrix[lower.tri(sim_matrix)] <- similarity_vec
    
    
    sim_matrix[upper.tri(sim_matrix)] <- 
      t(sim_matrix)[upper.tri(sim_matrix)]
    
    
    # TODO: rename this default update method?
    if(prior_mat_method == "new") {
      ## sim_matrix stays as is.
    } else if(prior_mat_method == "average") {
      if(!is.null(prior_mat)){
        ## incorporate new matrix into running average.
        sim_matrix <-
          1/(b+1) * (b * prior_mat + sim_matrix)
      }
    } else{stop("Provide 'average' or 'update' for arg prior_mat_method")}
    
    
    ## compare to prior matrix
    if(!is.null(prior_mat)) {
      ## compute normalized Frob norm of P - S
      norm_prior <- norm(prior_mat, type = "F")
      norm_current <- norm(sim_matrix, type = "F")
      
      ## storing 1st metric:
      convergence_res$frob_norm_dist[b] <- 
        norm(prior_mat / norm_prior - sim_matrix / norm_current, type = "F")
      
      ## Spectral embedding similarity + convergence metrics:
      spect_clust_res_prior   <- spectralCluster(prior_mat,
                                                 n_clust = n_clust)
      spect_clust_res_current <- spectralCluster(sim_matrix, 
                                                 n_clust = n_clust)
      
      paired_dists_prior   <- stats::dist(spect_clust_res_prior$G)
      paired_dists_current <- stats::dist(spect_clust_res_current$G)
      
      ## storing 2nd metric:
      convergence_res$spectral_embed_corr[b] <-
        cor(paired_dists_prior, paired_dists_current)
      
      ## silhouette scores
      sil_scores_current <- cluster::silhouette(spect_clust_res_current$cluster, 
                                                paired_dists_current)[,3]
      
      ## storing 3rd metric:
      convergence_res$avg_spectral_sil_score[b] <- mean(sil_scores_current)
        
      
    }
    
    ## update prior
    prior_mat <- sim_matrix
    
    ## end of while loop
    b <- b + 1
  }
  
  
  
  ## close workers if needed
  if(parallelize == TRUE) {
    stopCluster(cl)
  }
  
  return(list(sim_matrix = sim_matrix, 
              convergence_res = convergence_res))
  
}



getHighlyVariableGenes <- function(expr_matrix, p_val_thres = 1e-3) {
  ## Source: http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
  gene_means <- rowMeans(expr_matrix)
  gene_vars  <- apply(expr_matrix, 1, var)
  gene_cv2s  <- gene_vars / gene_means^2
  
  ## fit regression line
  minMeanForFit <- unname( quantile( gene_means[ which( gene_cv2s > .3 ) ], .95 ) )
  useForFit <- gene_means >= minMeanForFit # & spikeins
  fit <- statmod::glmgam.fit(cbind(a0 = 1, 
                                   a1tilde = 1/gene_means[useForFit]),
                             gene_cv2s[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  
  ## add fit line
  # xg <- exp(seq( min(log(gene_means[gene_means>0])), max(log(gene_means)), length.out=1000 ))
  # vfit <- a1/xg + a0
  
  
  ## get predictions and distant genes from fit line 
  afit <- a1/gene_means+a0
  varFitRatio <- gene_vars/(afit*gene_means^2)
  varorder <- order(varFitRatio,decreasing=T)
  
  ## compute FDR controlled p-values to select subset of genes
  df <- ncol(expr_matrix) - 1
  pval <- pchisq(varFitRatio*df,
                 df=df,
                 lower.tail=F)
  adj.pval <- p.adjust(pval, "fdr")
  sigVariedGenes <- adj.pval < p_val_thres
  
  #table(sigVariedGenes)
  # par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(gene_means),log(gene_cv2s))
  # lines( log(xg), log(vfit), col="black", lwd=3 )
  # points(log(gene_means[sigVariedGenes]),
  #        log(gene_cv2s[sigVariedGenes]),col=2)
  
  return(which(sigVariedGenes))
}


spectralCluster <- function(sim_matrix, 
                            n_clust = NULL,
                            G_col = NULL,
                            nonzero_thres = 1e-3) {
  ## Args:
  ## G_col is the dimension of the spectral embedding
  if(is.null(n_clust)) {
    n_clust = which.max(Spectrum::estimate_k(sim_matrix, maxk = 10, showplots=FALSE)$z)
  }
  if(is.null(G_col)) {
    G_col = n_clust
  }
  # TODO: decide if using this
  # clust_test <- 
  #   Spectrum::cluster_similarity(sim_matrix, 
  #                                k = which.max(Spectrum::estimate_k(sim_matrix, maxk = 10, showplots=FALSE)$z), 
  #                                clusteralg = "GMM", 
  #                                specalg = "Ng")
  ## NG spectral clustering code from Spectrum package
  # dv <- 1/sqrt(rowSums(sim_matrix))
  # l <- dv * sim_matrix %*% diag(dv)
  # decomp <- eigen(l)
  # xi <- decomp$vectors[, 1:2]
  # yi <- xi/sqrt(rowSums(xi^2))
  
  m <- dim(sim_matrix)[1]
  
  ## compute normalized graph laplacian
  degree_vec <- rowSums(sim_matrix)
  
  normal_graph_laplace <- diag(1, m) - 
    diag(degree_vec^(-0.5)) %*% sim_matrix %*% 
    diag(degree_vec^(-0.5)) 
  
  eig <- eigen(normal_graph_laplace)
  
  nonzero_ind <- max(which(eig$values > nonzero_thres))
  
  #plot(rev(eig$values))
  
  ## get eigenvectors corresponding to smallest positive eigenvalues.
  G <- eig$vectors[,rev((nonzero_ind+1-G_col):nonzero_ind)]
  ## Normalize rows
  G <- G/sqrt(rowSums(G^2))
  clust_res <- kmeanspp(G, k = n_clust)
  
  clust_res$centers
  
  
  ## compute distances from cluster centroids
  centroid_dist <-
    sapply(seq_len(nrow(G)), function(i){
      clust = clust_res$cluster[i]
      sqrt(sum((G[i,] - clust_res$centers[clust,])^2)) 
    })
  
  ## compute euclidean distance matrix
  dists <- stats::dist(G, method = "euclidean")
  
  sil_scores <- cluster::silhouette(clust_res$cluster, dists)[,3]
  
  mean_sil_score <- mean(sil_scores)
  
  # G <- G %>% 
  #   as.data.frame() %>%
  #   mutate(cluster       = clust_res$cluster,
  #          centroid_dist = centroid_dist,
  #          sil_score     = sil_scores)
  
  return(list(G = G, 
              clusters       = clust_res$cluster,
              centroid_dist  = centroid_dist,
              sil_scores     = sil_scores,
              mean_sil_score = mean_sil_score))
}


  
sampleClusters <- function(clust_res, samp_size) {
  ## Args:
  ## clust_res: output of spectralCluster()
  
  N <- dim(clust_res$G)[1]
  
  ## error checking
  if(samp_size >= N) {stop("Sample size is larger than number of cells")}
  
  ## sample an even proportion of clusters
  fracs <- table(clust_res$clusters) / N
  samps <- round(fracs * samp_size)
  
  samples <-
    unlist(
      sapply(seq_len(max(clust_res$clusters)), 
             function(c) {
               ## Choose cluster points with largest sil. width
               in_clust <- which(clust_res$clusters == c)
               clust_ranks <- frank(-clust_res$sil_scores[in_clust])
               return(in_clust[which(clust_ranks <= samps[c])])
             }))
  
  is_sampled <- sapply(seq_len(N), function(x) x %in% samples)
  
  return(is_sampled)
}
  




## TODO Test of cluster sampling pipeline -------------------------------------
## Spearman correlation
# sim_matrix <- rcorr(as.matrix(expr_matrix),
#                     type = "spearman")$r
# diag(sim_matrix) <- 0
# 
# spect_clust_res <- spectralCluster(sim_matrix, n_clust = 10)
# 
# is_sampled <- sampleClusters(spect_clust_res, samp_size = 50)
# 
# G <- spect_clust_res$G %>% 
#   data.frame() %>% 
#   mutate(cluster = spect_clust_res$clusters,
#          is_sampled = is_sampled)
# 
# 
# ggplot() + 
#   geom_point(data=G, mapping = aes(x = X1, y = X2, 
#                                    color=as.factor(cluster), 
#                                    shape=is_sampled)) + 
#   labs(x = "First Spectral Eigenvector", y = "Second Spectral Eigenvector",
#        color = "Cluster", shape = "Is Sampled")
# 



