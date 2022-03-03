

#' SIMILARITY LEARNING FOR RAFSIL
#'
#' @param FeatMat matrix Feature matrix cells x genes
#' @param nrep integer The number of replicates to average over for RAFSIL1
#' @param oob logical Should out of bag proximity be used?
#' @param method string "RAFSIL1" or "RAFSIL2"
#' @return matrix of dissimilarities between cells
#' @importFrom randomForest randomForest
#' @importFrom fpc pamk
rafsilRF <- function(FeatMat,nrep=50,
                     oob=TRUE,
                     method="RAFSIL1",
                     parallelize = FALSE, 
                     n_cores = 1){
#==========================================================

  dis = NULL

  #- RAFSIL1
  #---------
  if(method=="RAFSIL1"){
    
    if(parallelize) {
      ## start clust
      cl <- makeCluster(n_cores) 
      
      ## export functions.
      clusterExport(cl, varlist=c("FeatMat", "oob"),
                    envir=environment())
      
      sims =
        parLapply(cl = cl, X = seq_len(nrep), 
                 fun = function(no_args) {
                   return(
                     randomForest::randomForest(x=FeatMat,y=NULL,
                                                proximity=TRUE,oob.prox=oob)$proximity
                   )})
      
      ## close clust
      stopCluster(cl)
      sim = Reduce("+", sims) / length(sims)
    } else{
      sims = replicate(nrep,randomForest(x=FeatMat,y=NULL,proximity=TRUE,oob.prox=oob)$proximity)
      sim  = apply(sims,c(1,2),mean)
    }
    
    dis  = (1-sim)

  #- RAFSIL 2
  #----------
  } else if( method == "RAFSIL2"){

    featFu <- function(i){
      cl   = fpc::pamk(FeatMat[,i])
      labs = as.factor(cl$pamobject$clustering)
      return(randomForest::randomForest(x=FeatMat[,-i], y=labs, proximity=TRUE,oob.prox=oob)$proximity)
    }
    
    
    if(parallelize) {
      ## start clust
      cl <- makeCluster(n_cores) 
      
      ## export functions.
      clusterExport(cl, varlist=c("FeatMat", "oob", "featFu"),
                    envir=environment())
      
      sims =
        parLapply(cl = cl, X = seq_len(ncol(FeatMat)), 
                  featFu)
      
      sims = array(unlist(sims), dim = c(nrow(FeatMat), nrow(FeatMat), length(sims)))
      
      ## close clust
      stopCluster(cl)
    } else{
      sims =  vapply(1:ncol(FeatMat),featFu,array(0,dim=c(nrow(FeatMat),nrow(FeatMat))))
    }
    
    sim  = apply(sims,c(1,2),mean)
    dis  = (1-sim)

  } else {
    stop("method not implemented")
  }

  return(dis)

}
