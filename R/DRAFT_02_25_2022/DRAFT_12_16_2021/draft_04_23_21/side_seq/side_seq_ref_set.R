
selectRefSet <- function(expr_matrix,
                         selection_method = "random",
                         size_ref_set = 500,
                         gamma = NULL,
                         dissim_matrix = NULL,
                         n_clust = 10,
                         n_neighbors = 15, 
                         min_dist = 0.01) {
  
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
    if(is.null(dissim_matrix)) {
      ## Spearman correlation
      sim_matrix <- rcorr(as.matrix(expr_matrix),
                          type = "pearson")$r
      diag(sim_matrix) <- 0
      dissim_matrix <- max(sim_matrix) - sim_matrix
      diag(dissim_matrix) <- 0
    }
    
    ## UMAP + kMeans
    umap_embed <-
      uwot::umap(stats::as.dist(dissim_matrix),
                 n_neighbors = n_neighbors,
                 min_dist = min_dist)
    
    clust_res <- kmeanspp(umap_embed, k = n_clust)$cluster
    
    ## sampling from cluster results
    cells <- sampleClusters(clust_res, 
                            samp_size = size_ref_set,
                            N = N)
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
                          similarity_method = "n_intersect",
                          rbo_p = NULL,
                          n_top_genes = 300,
                          ## cell reference set selection parameters
                          selection_method = "random",
                          gamma = NULL,
                          size_ref_set,
                          n_clust = 10,
                          B = 3,
                          D = 5,
                          purity_tol = 0.85, ## minimum purity score for algorithm to complete.
                          n_neighbors = 15, 
                          min_dist = 0.01,
                          clust_repeats = 5,
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
                         side_seq_method = SIDESeq_method,
                         rbo_p = rbo_p)
  
  
  
  ## Iterate procedure B times or until convergence criteria is reached
  b <- 1
  prior_mat <- NULL
  dissim_matrix <- NULL
  
  ## storing the convergence criterion
  purity_scores <- rep(0, B)
  
  ### Dissimilarity Matrix Computation Workflow:
  computeDissim <- function() {
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
                                  "similarityMeasure", "rbo"),
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
    dissimilarity_vec <- max(similarity_vec) - similarity_vec
    
    ## convert vector to cell-by-cell matrix
    dissim_matrix <- matrix(0, nrow = N, ncol = N)
    
    dissim_matrix[lower.tri(dissim_matrix)] <- dissimilarity_vec
    
    
    dissim_matrix[upper.tri(dissim_matrix)] <- 
      t(dissim_matrix)[upper.tri(dissim_matrix)]
    
    return(dissim_matrix)
  }
  
  ### Run cluster convergence loop
  while(b <= B & purity_scores[max(1, b-1)] < purity_tol) {
    if(verbose == TRUE) {cat("Running convergence iteration", b, "of", B, "...")}
    ## extract reference set
    ref_cells <-
      selectRefSet(expr_matrix,
                   selection_method = selection_method,
                   size_ref_set = size_ref_set,
                   gamma = gamma,
                   dissim_matrix = dissim_matrix,
                   n_clust,
                   n_neighbors = n_neighbors,
                   min_dist = min_dist)
    
    ## compute dissim matrix
    dissim_matrix <- computeDissim()
    
    ## compare to prior matrix
    if(!is.null(prior_mat)) {
      
      internal_purity_vec <- c()
      for(r in seq_len(clust_repeats)) {
        ## UMAP embedding + Kmeans
        umap_embed_current <-
          uwot::umap(stats::as.dist(dissim_matrix),
                     n_neighbors = n_neighbors,
                     min_dist = min_dist)
        
        umap_embed_prior <- 
          uwot::umap(stats::as.dist(prior_mat),
                     n_neighbors = n_neighbors,
                     min_dist = min_dist)
        
        clust_res_current <- kmeanspp(umap_embed_current, k = n_clust)$cluster
        clust_res_prior   <- kmeanspp(umap_embed_prior,   k = n_clust)$cluster
        
        ## storing purity of current cluster with previous clustering:
        internal_purity_vec[r] <-
          clusterPurity(clust_res_current, 
                        clust_res_prior) 
      }
      purity_scores[b] <- mean(internal_purity_vec)
      if(verbose == TRUE){cat("Current clustering purity score", purity_scores[b])}
    }
    
    ## update prior
    prior_mat <- dissim_matrix
    
    ## end of while loop
    b <- b + 1
  }

  
  dissim_final <- NULL
  
  ## Run final averaging loop
  ref_cell_list <- 
    lapply(seq_len(D),
           function(d) {
             return(selectRefSet(expr_matrix,
                                 selection_method = selection_method,
                                 size_ref_set = size_ref_set,
                                 gamma = gamma,
                                 dissim_matrix = dissim_matrix,
                                 n_clust,
                                 n_neighbors = n_neighbors,
                                 min_dist = min_dist))
           })
  
  
  
  ## clear space
  rm(prior_mat,dissim_matrix)
  
  for(d in seq_len(D)) {
    if(verbose == TRUE) {cat("Running averaging iteration", d, "of", D, "...")}
    
    ref_cells <- ref_cell_list[[d]]
    dissim_curr <- computeDissim()
    
    if(is.null(dissim_final)){dissim_final <- dissim_curr}
    ## incorporate new matrix into running average.
    dissim_final <-
      1/(d+1) * (d * dissim_final + dissim_curr)
    
  }
  
  
  ## close workers if needed
  if(parallelize == TRUE) {
    stopCluster(cl)
  }
  
  return(list(dissim_final = dissim_final, 
              purity_scores = purity_scores))
  
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





sampleClusters <- function(clust_res, samp_size, N) {
  ## Args:
  ## clust_res: output of spectralCluster()
  
  
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




clusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}
