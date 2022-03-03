### SOURCE: https://github.com/cran/bigpca/blob/master/R/bigpca.R

quick_elbow <- function(varpc, 
                        low = 0.025, 
                        max_expl = 1) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low >= max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  
  lowie <- (ee < low) ; highie <- ee > low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      
      ## changed from abs, as in k-means the totwithinss can increase from two local mins.
      drops <- -1*diff(set) / (set[1:(length(set)-1)])
      infz <- is.infinite(drops)
      
      elbow <- which(drops == max(drops[!infz],na.rm=T))[1] + others
    }
  } else { 
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max_expl) {
    elbow <- which(cumsum(ee)>max_expl)[1]-1
  }
  if(elbow < 1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}
