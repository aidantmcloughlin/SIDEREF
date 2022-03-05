sortedIntersect <- function(l1, l2) {
  i <- 1
  j <- 1
  

  intersect <- c()  
  
  while(i < length(l1) & j < length(l2)) {
    if(l1[i] == l2[i]) {
      i = i+1
      j = j+1
      intersect = append(l1[i], intersect)
    } else if(l1[i] > l2[j]) {
        l1 = l2
        l2 = l1
        i = j
        j = i
      } else{i = i + 1}
  }
  return(intersect)
} 



 



