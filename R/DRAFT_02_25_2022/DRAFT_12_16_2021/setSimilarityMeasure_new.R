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
      similarityMeasure <- function(i, j, N, J,
                                    ref_set,
                                    d_e_mat,
                                    contribs_val) {

        gene_contrib_vec <- rep(0, J)


        sim_score <- list(0)
        if(!is.null(contribs_val)) {sim_score <- list(0, gene_contrib_vec)}
         
         
        R <- length(ref_set)
        R_trim <- length(setdiff(ref_set, c(i,j)))
        R_trim_inv <- 1 / R_trim

        for(r in seq_len(R)) {
          if(ref_set[r] %in% c(i,j)){
            next
          } else{
            inter <- base::intersect(
              d_e_mat[, ((r-1) * N + i)],
              d_e_mat[, ((r-1) * N + j)])

            ## update main similarity score.

            sim_score[[1]] <- sim_score[[1]] + length(inter) / R_trim
            # include gene contribs if specified.
            if(!is.null(contribs_val)) {
              sim_score[[2]][inter] <- sim_score[[2]][inter] + R_trim_inv
            }
          }


        }
        
        # finalizing contribution vector.
        if(!is.null(contribs_val)) {
          gene_contribs <-
            order(-sim_score[[2]], sample(seq_len(length(gene_contrib_vec))))[1:contribs_val]
          
          

          ## trim non-contributors from the list
          n_contributors <- sum(sim_score[[2]] > 1e-9)

          if(length(gene_contribs) > n_contributors) {
            gene_contribs <- gene_contribs[1:n_contributors]}
          if(n_contributors == 0) {gene_contribs <- 0}
        }
        
        ## keep scores as well.
        gene_contribs_scores <- sim_score[[2]][gene_contribs]
        
        ## return statement.
        if(!is.null(contribs_val)) {
          return(list(sim_score[[1]],
                      gene_contribs,
                      gene_contribs_scores))
        } else{return(sim_score[[1]])}
        

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



rbo <- function(x, y, p = 0.9) {
  assertthat::are_equal(length(x), length(y))
  return(
    (1-p) * sum(sapply(seq_len(length(x)),
                       function(d) {
                         a = 1/d * length(intersect(x[1:d], y[1:d]))
                         return(p^(d-1) * a)
                       }))
  )
}


