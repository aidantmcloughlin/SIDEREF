setSimilarityMeasure <- function(similarity_method, 
                                 index_method,
                                 do_sort) {
  
  if(similarity_method == "n_intersect") {
    if(index_method == "simultaneous") {
      similarityMeasure <- function(i, j, N, d_e_mat) {
        return(
          mean(sapply(setdiff(1:N, c(i,j)),
                     function(t) {
                       length(base::intersect(
                         d_e_mat[, getCombnPairIndex(i, t, N)],
                         d_e_mat[, getCombnPairIndex(j, t, N)]
                       ))}))
        )
      }
    } else if(index_method == "iterative") {
      similarityMeasure <- function(i, j, N, d_e_mat) {
        return(
          sum(sapply(1:(N-2),
                     function(t) {
                       length(base::intersect(
                         d_e_mat[, t],
                         d_e_mat[, (N-2+t)]
                       ))
                     }))
        )
      } 
    } else if(index_method == "SIDEREF") {
      if(do_sort) {
        similarityMeasure <- function(i, j, N,
                                      ref_set, 
                                      d_e_mat) {
          return(
            mean(sapply(seq_len(length(ref_set)),
                        function(t) {
                          ifelse(
                            ref_set[t] %in% c(i,j), NA,
                            sortedInterSum(
                              v1 = d_e_mat[, ((t-1) * N + i)],
                              v2 = d_e_mat[, ((t-1) * N + j)]))}), na.rm = TRUE)
            
          )
        }
      } else{
        similarityMeasure <- function(i, j, N,
                                      ref_set, 
                                      d_e_mat) {
          return(
            
            mean(sapply(seq_len(length(ref_set)),
                        function(t) {
                          ifelse(ref_set[t] %in% c(i,j), 
                                 NA,
                                 length(base::intersect(
                                   d_e_mat[, ((t-1) * N + i)],
                                   d_e_mat[, ((t-1) * N + j)])))
                        }), 
                 na.rm = TRUE)
          )
        } 
      }
    }
  }

  return(similarityMeasure)
}


sortedInterSum <- function(v1, v2) {
  c <- 0
  a <- 1
  b <- 1
  l <- length(v1)
  
  while(max(a, b) <= l) {
    if(identical(v1[a], v2[b])) {
      c <- c+1
      a <- a+1
      b <- b+1
    } else if(v1[a] < v2[b]) {
      a <- a + 1
    } else {
      b <- b + 1
    } 
  }
  
  return(c)
}

getCombnPairIndex <- function(i, j, N) {
  if (i < j) {
    return(i * (N - 0.5*(i + 1)) + j - N)
  }
  else {
    return(j * (N - 0.5*(j + 1)) + i - N)
  }
}


