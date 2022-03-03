#' perform the RAFSIL dissimilarity learning and clustering algorithm
#'
#' @title RAFSIL
#'
#' @param data an (n x m) data matrix of gene expression measurements of single cells
#' @param NumC number of clusters to be extracted over data
#' @param method character "RAFISIL1" or "RAFSIL2"
#' @param gene_filter logical should gene filtering be performed?
#' @param nrep integer for RAFSIL1 the number of times similarities are learned
#' @return list of 3 elements describing the results:
#'		SF:  Spearman Feature Space,
#'  	D:   RAFSIL dissimilarity matrix,
#'    lab: inferred cluster labels
#'  	hc:  hierarchical clustering solution
#' @importFrom e1071 kurtosis
#' @importFrom dclone make.symmetric
#' @import stats
#' @export 
#'
RAFSIL<-function(data,NumC=NULL,nrep=50,method="RAFSIL1",gene_filter=TRUE,
                 parallelize = FALSE, n_cores = 1)

  {

  #- feature construction
  FE    <- rafsilFE(data,gene_filter=gene_filter,frq=0.06)

  #- random forest learning
  rfdis <- rafsilRF(FE$FeatMat, nrep=nrep, method = method,
                    parallelize = parallelize, n_cores = n_cores)

  dis<-rfdis

  if(method=="RAFSIL2")

  {

  #- "normalized" distance
  X            <- rfdis
  X            <- t(t(X) - colMeans(X))
  X            <- X / max(abs(X))
  rfdis2       <- make.symmetric(X)
  rfdis2       <- rfdis2 +  abs(min(rfdis2))

  #- normalize?
  t1<-abs(kurtosis(as.vector(rfdis)))
  t2<-abs(kurtosis(as.vector(rfdis2)))
  if (t1>5*t2){
    diag(rfdis2) <- 0
    dis<-rfdis2
  }else{
    dis<-rfdis}
}

hc = NULL
lab = NULL
  if(!is.null(NumC)){
	   hc  <- hclust(as.dist(dis),method="average")
  	lab <- cutree(hc,NumC)
  }

  res = list(  SF  = FE,
	             D   = rfdis, #- never normalize for now
	             lab = lab,
	             hc  = hc)

   return(res)
}
