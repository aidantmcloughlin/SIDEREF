## Load Libraries Used in Analysis
## SideSEQ implementation packages
library(utils)
library(parallel)
library(doParallel)
## Data cleaning, plotting
library(dplyr)
library(pryr)
library(tidyr)
library(ggplot2)
library(janitor)
library(stringr)
library(data.table)
library(here)
library(magrittr)
library(forcats)
library(GGally)
library(ggpubr)
library(readxl)
library(Matrix)
library(purrr)

## BioConductor packages
library(Seurat)
library(splatter)
library(scater)
library(SummarizedExperiment)
library(scRNAseq)
library(fgsea)

### other clustering / statistics packages
library(anocva)
library(diffusr)
library(statmod)
library(LICORS)
library(cluster)
library(mclust)
library(irlba)
library(uwot)
library(NetworkToolbox)

## RAFSIL dependencies:
library(pracma)
library(ClusterR)
library(randomForest)
library(fpc)

## Load scRNA seq dissimilarity alternatives
library(SIMLR)

## common paste helper function:
`%p%` <- function(x,y) {paste0(x,y)}

## RAFSIL local repository source files
rafsil_direct <- here("R/RAFSIL_scripts/")
for(f in list.files(rafsil_direct)) {
  source(rafsil_direct %p% f)
}



loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


## Load Functions:
source(here('R/computeDiffExprMat.R'))
source(here('R/setSimilarityMeasure.R'))
source(here('R/quick_elbow.R'))
source(here('R/SIDEseq.R'))
source(here('R/SIDEREF.R'))


## Constants:
## parallelization setting
N_CORES <- 3

## CONSTANTS 
USE_DROPLET_ONLY <- TRUE
VAR_FEATURES <- 3000
G <- 300
R <- 100

RAFSIL_NREP <- 20

TM_SAMP_SIZE <- 100

##ie, using elbow criteria:
N_CLUST <- NULL
N_PCS <- NULL
PC_ELBOW_CHG_THRES <- 0.025
KM_ELBOW_CHG_THRES <- 0.025

## params for supplementary analyses
## TODO: update these once the code functions fine.
CONVERGE_RES_REPS <- 3
SNN_K <- 15
SNN_REPS <- 5
SNN_SAMP_SIZE <- 1500
REF_SET_SIZES <- c(25, 50, 75, 100, 150, 200, 300)

COMPUTE_TIME_N_CORES <- N_CORES
N_CELLS_SIZES <- c(100, 250, 500, 750, 1000, 1500)
COMPUTE_TIME_REPS <- 1






