### UWOT Umap package changes:
custom_umap <-
function (X, n_neighbors = 15, n_components = 2, metric = "euclidean", 
          n_epochs = NULL, learning_rate = 1, scale = FALSE, init = "spectral", 
          init_sdev = NULL, spread = 1, min_dist = 0.01, set_op_mix_ratio = 1, 
          local_connectivity = 1, bandwidth = 1, repulsion_strength = 1, 
          negative_sample_rate = 5, a = NULL, b = NULL, nn_method = NULL, 
          n_trees = 50, search_k = 2 * n_neighbors * n_trees, approx_pow = FALSE, 
          y = NULL, target_n_neighbors = n_neighbors, target_metric = "euclidean", 
          target_weight = 0.5, pca = NULL, pca_center = TRUE, pcg_rand = TRUE, 
          fast_sgd = FALSE, ret_model = FALSE, ret_nn = FALSE, ret_extra = c(), 
          n_threads = NULL, n_sgd_threads = 0, grain_size = 1, tmpdir = tempdir(), 
          verbose = getOption("verbose", TRUE)) {
  uwot(X = X, n_neighbors = n_neighbors, n_components = n_components, 
       metric = metric, n_epochs = n_epochs, alpha = learning_rate, 
       scale = scale, init = init, init_sdev = init_sdev, spread = spread, 
       min_dist = min_dist, set_op_mix_ratio = set_op_mix_ratio, 
       local_connectivity = local_connectivity, bandwidth = bandwidth, 
       gamma = repulsion_strength, negative_sample_rate = negative_sample_rate, 
       a = a, b = b, nn_method = nn_method, n_trees = n_trees, 
       search_k = search_k, method = "umap", approx_pow = approx_pow, 
       n_threads = n_threads, n_sgd_threads = n_sgd_threads, 
       grain_size = grain_size, y = y, target_n_neighbors = target_n_neighbors, 
       target_weight = target_weight, target_metric = target_metric, 
       pca = pca, pca_center = pca_center, pcg_rand = pcg_rand, 
       fast_sgd = fast_sgd, ret_model = ret_model || "model" %in% 
         ret_extra, ret_nn = ret_nn || "nn" %in% ret_extra, 
       ret_fgraph = "fgraph" %in% ret_extra, tmpdir = tempdir(), 
       verbose = verbose)
}









# Function that does all the real work
uwot <- function(X, n_neighbors = 15, n_components = 2, metric = "euclidean",
                 n_epochs = NULL,
                 alpha = 1, scale = FALSE,
                 init = "spectral", init_sdev = NULL,
                 spread = 1, min_dist = 0.01,
                 set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                 bandwidth = 1.0, gamma = 1.0,
                 negative_sample_rate = 5.0, a = NULL, b = NULL,
                 nn_method = NULL, n_trees = 50,
                 search_k = 2 * n_neighbors * n_trees,
                 method = "umap",
                 perplexity = 50, approx_pow = FALSE,
                 y = NULL, target_n_neighbors = n_neighbors,
                 target_metric = "euclidean",
                 target_weight = 0.5,
                 n_threads = NULL,
                 n_sgd_threads = 0,
                 grain_size = 1,
                 kernel = "gauss",
                 ret_model = FALSE, ret_nn = FALSE, ret_fgraph = FALSE,
                 pca = NULL, pca_center = TRUE,
                 pcg_rand = TRUE,
                 fast_sgd = FALSE,
                 tmpdir = tempdir(),
                 verbose = getOption("verbose", TRUE)) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (method == "umap" && (is.null(a) || is.null(b))) {
    ab_res <- find_ab_params(spread = spread, min_dist = min_dist)
    a <- ab_res[1]
    b <- ab_res[2]
    tsmessage("UMAP embedding parameters a = ", formatC(a), " b = ", formatC(b))
  }
  
  if (n_neighbors < 2) {
    stop("n_neighbors must be >= 2")
  }
  
  if (set_op_mix_ratio < 0.0 || set_op_mix_ratio > 1.0) {
    stop("set_op_mix_ratio must be between 0.0 and 1.0")
  }
  
  if (local_connectivity < 1.0) {
    stop("local_connectivity cannot be < 1.0")
  }
  
  if (!is.null(y) && is.numeric(y) && any(is.na(y))) {
    stop("numeric y cannot contain NA")
  }
  
  if (!is.numeric(n_components) || n_components < 1) {
    stop("'n_components' must be a positive integer")
  }
  
  if (!is.null(pca)) {
    if (!is.numeric(pca) || pca < 1) {
      stop("'pca' must be a positive integer")
    }
    if (pca < n_components) {
      stop("'pca' must be >= n_components")
    }
    if (pca > min(nrow(X), ncol(X))) {
      stop("'pca' must be <= min(nrow(X), ncol(X))")
    }
  }
  
  if (fast_sgd) {
    n_sgd_threads <- "auto"
    pcg_rand <- FALSE
    approx_pow <- TRUE
  }
  
  if (n_threads < 0) {
    stop("n_threads cannot be < 0")
  }
  if (n_threads %% 1 != 0) {
    n_threads <- round(n_threads)
    tsmessage("Non-integer 'n_threads' provided. Setting to ", n_threads)
  }
  if (n_sgd_threads == "auto") {
    n_sgd_threads <- n_threads
  }
  if (n_sgd_threads < 0) {
    stop("n_sgd_threads cannot be < 0")
  }
  if (n_sgd_threads %% 1 != 0) {
    n_sgd_threads <- round(n_sgd_threads)
    tsmessage("Non-integer 'n_sgd_threads' provided. Setting to ", n_sgd_threads)
  }
  
  # Store categorical columns to be used to generate the graph
  Xcat <- NULL
  # number of original columns in data frame (or matrix)
  # will be used only if using df or matrix and ret_model = TRUE
  norig_col <- NULL
  if (is.null(X)) {
    if (!is.list(nn_method)) {
      stop("If X is NULL, must provide NN data in nn_method")
    }
    if (is.character(init) && tolower(init) %in% c("spca", "pca")) {
      stop("init = 'pca' and 'spca' can't be used with X = NULL")
    }
    n_vertices <- x2nv(nn_method)
  }
  else if (methods::is(X, "dist")) {
    if (ret_model) {
      stop("Can only create models with dense matrix or data frame input")
    }
    checkna(X)
    n_vertices <- attr(X, "Size")
    tsmessage("Read ", n_vertices, " rows")
  }
  else if (methods::is(X, "sparseMatrix")) {
    if (ret_model) {
      stop("Can only create models with dense matrix or data frame input")
    }
    checkna(X)
    n_vertices <- nrow(X)
    if (ncol(X) != n_vertices) {
      stop("Sparse matrices are only supported as distance matrices")
    }
    tsmessage("Read ", n_vertices, " rows of sparse distance matrix")
  }
  else {
    cat_ids <- NULL
    norig_col <- ncol(X)
    if (methods::is(X, "data.frame") || methods::is(X, "matrix")) {
      if (methods::is(X, "matrix")) {
        X <- data.frame(X)
      }
      cat_res <- find_categoricals(metric)
      metric <- cat_res$metrics
      cat_ids <- cat_res$categoricals
      # Convert categorical columns to factors if they aren't already
      if (!is.null(cat_ids)) {
        X[, cat_ids] <- lapply(X[, cat_ids, drop = FALSE], factor)
        Xcat <- X[, cat_ids, drop = FALSE]
      }
      
      indexes <- which(vapply(X, is.numeric, logical(1)))
      if (length(indexes) == 0) {
        stop("No numeric columns found")
      }
      X <- as.matrix(X[, indexes])
    }
    checkna(X)
    n_vertices <- nrow(X)
    tsmessage(
      "Read ", n_vertices, " rows and found ", ncol(X),
      " numeric columns",
      appendLF = is.null(cat_ids)
    )
    if (length(cat_ids) > 0) {
      tsmessage(" and ", pluralize("categorical column", length(cat_ids)),
                time_stamp = FALSE
      )
    }
    X <- scale_input(X,
                     scale_type = scale, ret_model = ret_model,
                     verbose = verbose
    )
  }
  
  if (method == "largevis" && kernel == "knn") {
    n_neighbors <- perplexity
  }
  
  if (n_neighbors > n_vertices) {
    # If nn_method is a list, we will determine n_neighbors later
    if (!is.list(nn_method)) {
      # Otherwise,for LargeVis, n_neighbors normally determined from perplexity
      # not an error to be too large
      if (method == "largevis") {
        tsmessage("Setting n_neighbors to ", n_vertices)
        n_neighbors <- n_vertices
      }
      else {
        stop("n_neighbors must be smaller than the dataset size")
      }
    }
  }
  
  if (!is.list(metric)) {
    metrics <- list(c())
    names(metrics) <- metric
  }
  else {
    metrics <- metric
  }
  
  # For typical case of numeric matrix X and not using hamming distance, save
  # PCA results here in case initialization uses PCA too
  pca_models <- NULL
  pca_shortcut <- FALSE
  if (!is.null(pca) && length(metric) == 1 && metric != "hamming" &&
      is.matrix(X) && ncol(X) > pca) {
    tsmessage("Reducing X column dimension to ", pca, " via PCA")
    pca_res <- pca_scores(X,
                          ncol = pca, center = pca_center,
                          ret_extra = ret_model, verbose = verbose
    )
    if (ret_model) {
      X <- pca_res$scores
      pca_models[["1"]] <- pca_res[c("center", "rotation")]
      pca_res <- NULL
    }
    else {
      X <- pca_res
    }
    pca_shortcut <- TRUE
  }
  
  d2sr <- data2set(X, Xcat, n_neighbors, metrics, nn_method,
                   n_trees, search_k,
                   method,
                   set_op_mix_ratio, local_connectivity, bandwidth,
                   perplexity, kernel,
                   n_threads, grain_size,
                   ret_model,
                   pca = pca, pca_center = pca_center,
                   n_vertices = n_vertices,
                   tmpdir = tmpdir,
                   verbose = verbose
  )
  
  V <- d2sr$V
  nns <- d2sr$nns
  if (is.null(pca_models)) {
    pca_models <- d2sr$pca_models
  }
  
  if (!is.null(y)) {
    tsmessage("Processing y data")
    
    if (!is.list(target_metric)) {
      target_metrics <- list(c())
      names(target_metrics) <- target_metric
    }
    else {
      target_metrics <- target_metric
    }
    
    ycat <- NULL
    ycat_ids <- NULL
    if (methods::is(y, "data.frame")) {
      ycat_res <- find_categoricals(target_metric)
      target_metric <- ycat_res$metrics
      ycat_ids <- ycat_res$categoricals
      if (!is.null(ycat_ids)) {
        ycat <- y[, ycat_ids, drop = FALSE]
      }
      else {
        ycindexes <- which(vapply(y, is.factor, logical(1)))
        if (length(ycindexes) > 0) {
          ycat <- (y[, ycindexes, drop = FALSE])
        }
      }
      
      yindexes <- which(vapply(y, is.numeric, logical(1)))
      
      if (length(yindexes) > 0) {
        y <- as.matrix(y[, yindexes])
      }
      else {
        y <- NULL
      }
    }
    else if (is.list(y)) {
      nn_method <- y
    }
    else if (is.numeric(y)) {
      y <- as.matrix(y)
    }
    else if (is.factor(y)) {
      ycat <- data.frame(y)
      y <- NULL
    }
    
    
    if (!is.null(y)) {
      yd2sr <- data2set(y, ycat, target_n_neighbors, target_metrics, nn_method,
                        n_trees, search_k,
                        method,
                        set_op_mix_ratio = 1.0,
                        local_connectivity = 1.0,
                        bandwidth = 1.0,
                        perplexity = perplexity, kernel = kernel,
                        n_threads = n_threads, grain_size = grain_size,
                        ret_model = FALSE,
                        pca = pca,
                        n_vertices = n_vertices,
                        tmpdir = tmpdir,
                        verbose = verbose
      )
      
      tsmessage(
        "Intersecting X and Y sets with target weight = ",
        formatC(target_weight)
      )
      V <- set_intersect(V, yd2sr$V, target_weight, reset = TRUE)
      yd2sr$V <- NULL
      yd2sr$nns <- NULL
    }
    else if (!is.null(ycat)) {
      V <- categorical_intersection_df(ycat, V,
                                       weight = target_weight,
                                       verbose = verbose
      )
    }
  }
  
  if (!(ret_model || ret_nn)) {
    nns <- NULL
    gc()
  }
  
  if (methods::is(init, "matrix")) {
    if (nrow(init) != n_vertices || ncol(init) != n_components) {
      stop("init matrix does not match necessary configuration for X: ", "should
           have dimensions (", n_vertices, ", ", n_components, ")")
    }
    tsmessage("Initializing from user-supplied matrix")
    embedding <- init
  }
  else if (!(methods::is(init, "character") && length(init) == 1)) {
    stop("init should be either a matrix or string describing the ",
         "initialization method")
  }
  else {
    init <- match.arg(tolower(init), c(
      "spectral", "random", "lvrandom", "normlaplacian",
      "laplacian", "spca", "pca", "inormlaplacian", "ispectral",
      "agspectral"
    ))
    
    if (init_is_spectral(init)) {
      connected <- connected_components(V)
      if (connected$n_components > 1) {
        tsmessage("Found ", connected$n_components, " connected components, ", appendLF = FALSE)
        if (is.null(X)) {
          tsmessage("falling back to random initialization", time_stamp = FALSE)
          init <- "random"
        }
        else {
          tsmessage("falling back to 'spca' initialization with init_sdev = 1",
                    time_stamp = FALSE)
          init <- "spca"
          init_sdev <- 1
        }
      }
    }
    
    # Don't repeat PCA initialization if we've already done it once
    if (pca_shortcut && init %in% c("spca", "pca") && pca >= n_components) {
      embedding <- X[, 1:n_components]
      if (init == "spca") {
        tsmessage("Initializing from scaled PCA")
      }
      else {
        tsmessage("Initializing from PCA")
      }
    }
    else {
      embedding <- switch(init,
                          spectral = spectral_init(V, ndim = n_components, verbose = verbose),
                          random = rand_init(n_vertices, n_components, verbose = verbose),
                          lvrandom = rand_init_lv(n_vertices, n_components, verbose = verbose),
                          normlaplacian = normalized_laplacian_init(V,
                                                                    ndim = n_components,
                                                                    verbose = verbose
                          ),
                          laplacian = laplacian_eigenmap(V, ndim = n_components, verbose = verbose),
                          # we handle scaling pca below
                          spca = pca_init(X, ndim = n_components, verbose = verbose),
                          pca = pca_init(X, ndim = n_components, verbose = verbose),
                          ispectral = irlba_spectral_init(V, ndim = n_components, verbose = verbose),
                          inormlaplacian = irlba_normalized_laplacian_init(V,
                                                                           ndim = n_components,
                                                                           verbose = verbose
                          ),
                          agspectral = agspectral_init(V,
                                                       n_neg_nbrs = negative_sample_rate,
                                                       ndim = n_components, verbose = verbose
                          ),
                          stop("Unknown initialization method: '", init, "'")
      )
    }
    if (!is.null(init_sdev) || init == "spca") {
      if (is.null(init_sdev)) {
        init_sdev <- 1e-4
      }
      embedding <- shrink_coords(embedding, init_sdev)
    }
  }
  
  if (is.null(n_epochs) || n_epochs < 0) {
    if (method == "largevis") {
      n_epochs <- lvish_epochs(n_vertices, V)
    }
    else {
      if (n_vertices <= 10000) {
        n_epochs <- 500
      }
      else {
        n_epochs <- 200
      }
    }
  }
  
  if (n_epochs > 0) {
    V@x[V@x < max(V@x) / n_epochs] <- 0
    V <- Matrix::drop0(V)
    
    epochs_per_sample <- make_epochs_per_sample(V@x, n_epochs)
    
    positive_head <- V@i
    positive_tail <- Matrix::which(V != 0, arr.ind = TRUE)[, 2] - 1
    
    tsmessage(
      "Commencing optimization for ", n_epochs, " epochs, with ",
      length(positive_head), " positive edges",
      pluralize("thread", n_sgd_threads, " using")
    )
    
    embedding <- t(embedding)
    if (tolower(method) == "umap") {
      embedding <- optimize_layout_umap(
        head_embedding = embedding,
        tail_embedding = NULL,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        a = a, b = b, gamma = gamma,
        initial_alpha = alpha, negative_sample_rate,
        approx_pow = approx_pow,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        move_other = TRUE,
        verbose = verbose
      )
    }
    else if (method == "tumap") {
      embedding <- optimize_layout_tumap(
        head_embedding = embedding,
        tail_embedding = NULL,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        initial_alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        move_other = TRUE,
        verbose = verbose
      )
    }
    else {
      embedding <- optimize_layout_largevis(
        head_embedding = embedding,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        gamma = gamma,
        initial_alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        verbose = verbose
      )
    }
    embedding <- t(embedding)
    gc()
    # Center the points before returning
    embedding <- scale(embedding, center = TRUE, scale = FALSE)
    tsmessage("Optimization finished")
  }
  
  if (ret_model || ret_nn || ret_fgraph) {
    nblocks <- length(nns)
    res <- list(embedding = embedding)
    if (ret_model) {
      res <- append(res, list(
        scale_info = if (!is.null(X)) { attr_to_scale_info(X) } else { NULL },
        n_neighbors = n_neighbors,
        search_k = search_k,
        local_connectivity = local_connectivity,
        n_epochs = n_epochs,
        alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        method = method,
        a = a,
        b = b,
        gamma = gamma,
        approx_pow = approx_pow,
        metric = metrics,
        norig_col = norig_col,
        pcg_rand = pcg_rand
      ))
      if (nblocks > 1) {
        res$nn_index <- list()
        for (i in 1:nblocks) {
          res$nn_index[[i]] <- nns[[i]]$index
        }
      }
      else {
        res$nn_index <- nns[[1]]$index
        if (is.null(res$metric[[1]])) {
          # 31: Metric usually lists column indices or names, NULL means use all
          # of them, but for loading the NN index we need the number of 
          # columns explicitly (we don't have access to the column dimension of
          # the input data at load time)
          # To be sure of the dimensionality, fetch the first item from the 
          # index and see how many elements are in the returned vector.
          if(!is.null(X)){
            rcppannoy <- get_rcppannoy(res$nn_index)
            res$metric[[1]] <- list(ndim = length(rcppannoy$getItemsVector(0)))
          } else {
            res$metric[[1]] <- list()
          }
        }
      }
      if (!is.null(pca_models)) {
        res$pca_models <- pca_models
      }
    }
    if (ret_nn) {
      res$nn <- list()
      for (i in 1:nblocks) {
        res$nn[[i]] <- list(idx = nns[[i]]$idx, dist = nns[[i]]$dist)
      }
      names(res$nn) <- names(nns)
    }
    if (ret_fgraph) {
      if (method == "largevis") {
        res$P <- V
      }
      else {
        res$fgraph <- V
      }
    }
  }
  else {
    res <- embedding
  }
  
  res
}








data2set <- function(X, Xcat, n_neighbors, metrics, nn_method,
                     n_trees, search_k,
                     method,
                     set_op_mix_ratio, local_connectivity, bandwidth,
                     perplexity, kernel,
                     n_threads, grain_size,
                     ret_model,
                     n_vertices = x2nv(X),
                     tmpdir = tempdir(),
                     pca = NULL, pca_center = TRUE,
                     verbose = FALSE) {
  V <- NULL
  nns <- list()
  nblocks <- length(metrics)
  
  # Check for precalculated NN data in nn_method
  if (is.list(nn_method)) {
    if (is.null(nn_method$idx)) {
      nblocks <- length(nn_method)
      if (nblocks == 0) {
        stop("Incorrect format for precalculated neighbor data")
      }
    }
    else {
      nblocks <- 1
      # wrap nn data in a list so data is always a list of lists
      nn_method <- list(nn_method)
    }
    metrics <- replicate(nblocks, NULL, simplify = FALSE)
    names(metrics) <- rep("precomputed", nblocks)
  }
  
  
  if (nblocks > 1) {
    tsmessage("Found ", nblocks, " blocks of data")
  }
  mnames <- tolower(names(metrics))
  if (is.null(nn_method)) {
    if (n_vertices < 4096 && !ret_model && all(mnames == "euclidean")) {
      tsmessage("Using FNN for neighbor search, n_neighbors = ", n_neighbors)
      nn_method <- "fnn"
    }
    else {
      tsmessage("Using Annoy for neighbor search, n_neighbors = ", n_neighbors)
      nn_method <- "annoy"
    }
  }
  
  pca_models <- list()
  for (i in 1:nblocks) {
    metric <- mnames[[i]]
    metric <- match.arg(metric, c(
      "euclidean", "cosine", "manhattan",
      "hamming", "correlation", "precomputed"
    ))
    # Defaults for this block which can be overridden
    pca_i <- pca
    pca_center_i <- pca_center
    
    subset <- metrics[[i]]
    if (is.null(subset)) {
      Xsub <- X
    }
    else if (is.list(subset)) {
      # e.g. "euclidean" = list(1:10, pca_center = FALSE),
      lsres <- lsplit_unnamed(subset)
      if (is.null(lsres$unnamed)) {
        stop("Error: no subset provided for block ", i)
      }
      if (length(lsres$unnamed) != 1) {
        stop("Error: only one unnamed item should be provided for block ", i)
      }
      subset <- lsres$unnamed[[1]]
      
      # possible overrides
      if (!is.null(lsres$named)) {
        lsnamed <- lsres$named
        lsnames <- names(lsnamed)
        if (!is.null(lsnamed$pca_center)) {
          pca_center_i <- lsnamed$pca_center
        }
        # PCA argument can be NULL, so need to check if it was explicitly provided
        if ("pca" %in% lsnames) {
          pca_i <- lsnamed$pca
        }
      }
      Xsub <- X[, subset, drop = FALSE]
    }
    else {
      Xsub <- X[, subset, drop = FALSE]
    }
    
    if (!is.null(X) && is.matrix(X)) {
      block_size <- ncol(Xsub)
      if (block_size == 0) {
        stop("Block ", i, " has zero size")
      }
      if (nblocks > 1) {
        tsmessage(
          "Processing block ", i, " of ", nblocks,
          " with size ", block_size,
          " using metric '", metric, "'"
        )
      }
    }
    else {
      # X is NULL or dist or something like that
      if (nblocks > 1) {
        tsmessage(
          "Processing block ", i, " of ", nblocks,
          " using metric '", metric, "'"
        )
      }
    }
    
    if (!is.null(pca_i) && is.matrix(X) && metric != "hamming" &&
        ncol(X) > pca_i && nrow(X) > pca_i) {
      tsmessage("Reducing column dimension to ", pca_i, " via PCA")
      pca_res <- pca_scores(Xsub, pca_i,
                            ret_extra = ret_model,
                            center = pca_center_i,
                            verbose = verbose
      )
      if (ret_model) {
        Xsub <- pca_res$scores
        pca_models[[as.character(i)]] <- pca_res[c("center", "rotation")]
        pca_res <- NULL
      }
      else {
        Xsub <- pca_res
      }
    }
    
    nn_sub <- nn_method
    # Extract this block of nn data from list of lists
    if (metric == "precomputed") {
      nn_sub <- nn_method[[i]]
      if (i == 1) {
        n_neighbors <- NULL
      }
      else {
        n_neighbors <- ncol(nn_method[[1]]$idx)
      }
    }
    
    x2set_res <- x2set(Xsub, n_neighbors, metric,
                       nn_method = nn_sub,
                       n_trees, search_k,
                       method,
                       set_op_mix_ratio, local_connectivity, bandwidth,
                       perplexity, kernel,
                       n_threads, grain_size,
                       ret_model,
                       n_vertices = n_vertices,
                       tmpdir = tmpdir,
                       verbose = verbose
    )
    Vblock <- x2set_res$V
    nn <- x2set_res$nn
    nns[[i]] <- nn
    names(nns)[[i]] <- metric
    n_neighbors <- ncol(nn$idx)
    if (is.null(V)) {
      V <- Vblock
    }
    else {
      V <- set_intersect(V, Vblock, weight = 0.5, reset = TRUE)
    }
  }
  
  if (!is.null(Xcat)) {
    V <- categorical_intersection_df(Xcat, V, weight = 0.5, verbose = verbose)
  }
  
  list(V = V, nns = nns, pca_models = pca_models)
}





x2nn <- function(X, n_neighbors, metric, nn_method,
                 n_trees, search_k,
                 tmpdir = tempdir(),
                 n_threads, grain_size,
                 ret_model,
                 n_vertices = x2nv(X),
                 verbose = FALSE) {
  if (is.list(nn_method)) {
    # on first iteration n_neighbors is NULL
    # on subsequent iterations ensure n_neighbors is consistent for all data
    validate_nn(nn_method, n_vertices, n_neighbors = n_neighbors)
    nn <- nn_method
  }
  else {
    nn_method <- match.arg(tolower(nn_method), c("annoy", "fnn"))
    if (nn_method == "fnn" && metric != "euclidean") {
      stop(
        "nn_method = 'FNN' is only compatible with distance metric ",
        "'euclidean'"
      )
    }
    if (nn_method == "fnn" && ret_model) {
      stop("nn_method = 'FNN' is incompatible with ret_model = TRUE")
    }
    nn <- find_nn(X, n_neighbors,
                  method = nn_method, metric = metric,
                  n_trees = n_trees, search_k = search_k,
                  tmpdir = tmpdir,
                  n_threads = n_threads, grain_size = grain_size,
                  ret_index = ret_model, verbose = verbose
    )
  }
  nn
}




find_nn <- function(X, k, include_self = TRUE, method = "fnn",
                    metric = "euclidean",
                    n_trees = 50, search_k = 2 * k * n_trees,
                    tmpdir = tempdir(),
                    n_threads = NULL,
                    grain_size = 1,
                    ret_index = FALSE,
                    verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (methods::is(X, "dist")) {
    res <- dist_nn(X, k, include_self = include_self)
  }
  else if (methods::is(X, "sparseMatrix")) {
    # sparse distance matrix
    if (Matrix::isTriangular(X)) {
      res <- sparse_tri_nn(X, k, include_self = include_self)
    }
    else {
      res <- sparse_nn(X, k, include_self = include_self)
    }
  }
  else {
    # normal matrix
    if (method == "fnn") {
      res <- FNN_nn(X, k = k, include_self = include_self)
    }
    else {
      res <- annoy_nn(X,
                      k = k,
                      metric = metric,
                      n_trees = n_trees, search_k = search_k,
                      tmpdir = tmpdir,
                      n_threads = n_threads,
                      ret_index = ret_index,
                      verbose = verbose
      )
    }
  }
  
  res
}





dist_nn <- function(X, k, include_self = TRUE) {
  X <- as.matrix(X)
  
  if (!include_self) {
    k <- k + 1
  }
  
  nn_idx <- t(apply(X, 2, order))[, 1:k]
  nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
  for (i in seq_len(nrow(nn_idx))) {
    nn_dist[i, ] <- X[i, nn_idx[i, ]]
  }
  
  if (!include_self) {
    nn_idx <- nn_idx[, 2:ncol(nn_idx)]
    nn_dist <- nn_dist[, 2:ncol(nn_dist)]
  }
  
  attr(nn_idx, "dimnames") <- NULL
  attr(nn_dist, "dimnames") <- NULL
  
  list(idx = nn_idx, dist = nn_dist)
}





