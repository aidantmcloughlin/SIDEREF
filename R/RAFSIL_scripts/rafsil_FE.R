

#' FEATURE ENGINEERING FOR RAFSIL
#'
#' @param texpr matrix Expression values, cells x genes
#' @param frq numeric Frequency filtering fraction
#' @param gene_filter logical shold gene filtering be performed?
#' @return list with two components: feature matrix F and the result of gene clustering cl.gene
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom fpc pamk
rafsilFE <- function(texpr,gene_filter=TRUE,frq=0.06){

  #- GENE FILTERING
  if (gene_filter==TRUE){
    texpr  = texpr[,colSums(texpr>0)>floor(frq*nrow(texpr)) ]
  }
  
  #- PCA EMBEDDING OF GENES
  # pc.gen  = rpca(t(texpr), retx=TRUE, center=TRUE, scale=TRUE)
  # nc      = max(which(cumsum(pc.gen$sdev)/sum(pc.gen$sdev)< .9))
  # nc      = length(elbow_detection((pc.gen$sdev[1:nc])))+1 #- could use log
  # gen.dat = scale(pc.gen$x[,1:nc])
  gen.dat<-fast.pca(scale(t(texpr),scale=TRUE,center=TRUE))

  #- K-MEANS CLUSTERING
  clres   = sapply(2:10,  function(nc)
                              sum(KMeans_rcpp(gen.dat, max_iters = 50,
                              clusters = nc,initializer="kmeans++")$WCSS_per_cluster))
  nc      = max(elbow_detection(log(clres)))+1
  cl.gene = KMeans_rcpp(gen.dat, clusters =nc, num_init = 5, max_iters = 100, initializer = 'kmeans++')

  #- SPEARMAN FEATURE CONSTRUCTION
  sfc <- function(mat){
  #=====================
    #- specaial-treat low-variance cells
    ind       = apply(mat,1,sd) < 1E-6
    mat[ind,] = t(apply(mat[ind,],1,rank,ties.method="random"))
    #- spearman distance
    d1 = sqrt(1-cor(t(mat),method="spearman"))
    #- pca embedding
    pc.cell = prcomp(d1, retx=TRUE, center=TRUE, scale=TRUE)
    # nc      = max(which(cumsum(pc.cell$sdev)/sum(pc.cell$sdev)< .9))
    # nc      = length(elbow_detection((pc.cell$sdev[1:nc])))+1 #- could use log
    pcc.var<-(pc.cell$sdev)^2
    nc<-length(elbow_detection(pcc.var))  #- be donezz
    return(pc.cell$x[,1:nc])
  }

  FeatMat    = sapply(1:nc,function(cind) sfc(texpr[,cl.gene$clusters==cind]))
  Fmat = matrix(unlist(FeatMat),nrow=nrow(texpr))

  #- components per cluster in names
  ncomps         = unlist(lapply(FeatMat,ncol))
  names          = sapply(1:nc,function(ci) paste("C",rep(ci,ncomps[ci]),"-F",1:ncomps[ci],sep=""))
  names          = unlist(names)
  colnames(Fmat) = names
  rownames(Fmat) = rownames(texpr)
  #- return FEATURE MATRIX
  return(list(
    FeatMat = Fmat,
    geneClusters = cl.gene
  ))

}
