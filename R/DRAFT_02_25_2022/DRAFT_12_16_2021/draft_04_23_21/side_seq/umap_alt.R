## Alternative umap package implementation

umap_alt <- 
  function (d, config = umap.defaults, method = c("naive", 
                                                  "umap-learn"), preserve.seed = TRUE, ...) 
  {
    method = config$method = match.arg(method)
    config = umap.check.config(config, ...)
    d = umap.prep.input(d, config)
    old.seed = get.global.seed()
    if (!is.na(config$random_state)) {
      set.seed(config$random_state)
    }
    if (nrow(d) <= 2) {
      result = umap.small(d, config)
    }
    else {
      implementations = c(naive = umap.naive, `umap-learn` = umap.learn)
      if (method %in% names(implementations)) {
        result = implementations[[method]](d, config)
      }
    }
    class(result) = "umap"
    if (preserve.seed) {
      set.global.seed(old.seed)
    }
    result
  }

