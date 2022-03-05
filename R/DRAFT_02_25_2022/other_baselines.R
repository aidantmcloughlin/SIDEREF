###############################################################################
### OTHER BASELINES
###############################################################################


### load modules
library(SIMLR)

## RAFSIL local repository source files
rafsil_direct <- here("R/RAFSIL_scripts")
for(f in list.files(rafsil_direct)) {
  source(rafsil_direct %p% f)
}



### 1: SIMLR
start_time_simlr <- Sys.time()
simlr_res <- SIMLR::SIMLR(expr_matrix, c=10)
end_time_simlr <- Sys.time()

end_time_simlr - start_time_simlr

## convert to distance
dist_simlr <- max(simlr_res$S) - simlr_res$S
diag(dist_simlr) <- 0




### 2: RAFSIL

## NOTE: input data is n x m, where rows are cells, columns genes.
## Repetitions is how many times to run the random forest.


raf_res1 <-
  RAFSIL(t(expr_matrix), NumC=NULL, nrep=1, 
         method="RAFSIL1", gene_filter = FALSE,
         parallelize=FALSE)


dist_rafsil <- raf_res1$D





### TODO: can also consider binarized matrix with combinatoric gene pairs + Hamming 
##          (relevant article finds performance advantages when decomposing to gene modules, and without doing so, 
##              can't do all combinations in decent time / lots of irrelevant rows.)




