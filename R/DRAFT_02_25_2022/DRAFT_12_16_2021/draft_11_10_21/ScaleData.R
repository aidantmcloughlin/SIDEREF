### Directly Unplugging Seurat Data Scaling Function for use on Splatter.
SeuratScaleData <-
  function (mat, center = TRUE, scale = TRUE, scale_max = 10) 
  {
    if (center) {
      rm <- rowMeans2(x = mat, na.rm = TRUE)
    }
    if (scale) {
      if (center) {
        rsd <- rowSds(mat, center = rm)
      }
      else {
        rsd <- sqrt(x = rowSums2(x = mat^2)/(ncol(x = mat) - 1))
      }
    }
    if (center) {
      mat <- mat - rm
    }
    if (scale) {
      mat <- mat/rsd
    }
    if (scale_max != Inf) {
      mat[mat > scale_max] <- scale_max
    }
    return(mat)
  }