setSimilarityMeasure <- function(similarity_method, 
                                 side_seq_method,
                                 rbo_p = NULL) {
  
  if(similarity_method == "n_intersect") {
    if(side_seq_method == "simultaneous") {
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
    } else if(side_seq_method == "iterative") {
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
    } else if(side_seq_method == "block") {
      ## TODO: once/if this is implemented.
      stop("Not implemented properly at the moment.")
    } else if(side_seq_method == "ref_set") {
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
  else if(similarity_method == "rbo" & side_seq_method == "ref_set") {
    if(is.null(rbo_p)) {rbo_p = 0.9}
    similarityMeasure <- function(i, j, N,
                                  ref_set, 
                                  d_e_mat) {
      return(
        mean(sapply(seq_len(length(ref_set)),
                    function(t) {
                      ifelse(ref_set[t] %in% c(i,j), 
                             NA,
                             rbo(d_e_mat[, ((t-1) * N + i)],
                                 d_e_mat[, ((t-1) * N + j)],
                                 p = rbo_p))
                    }), 
             na.rm = TRUE)
      )
    }
  }
  return(similarityMeasure)
}




getCombnPairIndex <- function(i, j, N) {
  if (i < j) {
    return(i * (N - 0.5*(i + 1)) + j - N)
  }
  else {
    return(j * (N - 0.5*(j + 1)) + i - N)
  }
}



