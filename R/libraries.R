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
library(ggnewscale)
library(gridExtra)
library(grid)
library(network)
library(gtable)
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
library(aricode)
library(irlba)
library(uwot)
library(NetworkToolbox)
library(aricode)

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
N_CORES <- 4

## CONSTANTS 
USE_DROPLET_ONLY <- TRUE
VAR_FEATURES <- 3000
G <- 300
R <- 100

RAFSIL_NREP <- 20

TM_SAMP_SIZE <- 100

## Distance setting constants:
SPLAT_SPECTRAL_DIMS <- c(5, 10, 25)
SPLAT_PCS <- c(25)

TM_SPECTRAL_DIMS <- c(3, 5, 10, 25)
TM_PCS <- c(25)


##ie, using elbow criteria:
N_CLUST <- NULL
N_PCS <- NULL
PC_ELBOW_CHG_THRES <- 0.025
KM_ELBOW_CHG_THRES <- 0.025

MIN_DIST <- 0.01
N_NEIGHBORS <- 15

## params for supplementary analyses
CONVERGE_RES_REPS <- 3
SNN_K <- 15
SNN_REPS <- 5
SNN_SAMP_SIZE <- 1500
REF_SET_SIZES <- c(25, 50, 75, 100, 150, 200, 300)

COMPUTE_TIME_N_CORES <- 8
N_CELLS_SIZES <- c(250, 500, 750, 1000, 1250, 1500, 2000)
COMPUTE_TIME_REPS <- 1

## params for computing performance metrics
RUN_SPLATTER_PERFORMANCE_METRICS <- TRUE
N_SPLATTER_SIMS_METRICS <- 10

## Plotting:
DPI <- 600

