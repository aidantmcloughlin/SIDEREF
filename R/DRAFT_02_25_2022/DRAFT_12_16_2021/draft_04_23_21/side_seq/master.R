## Load Libraries Used in Analysis
## SideSEQ implementation packages
library(utils)
library(parallel)

## Data cleaning, plotting
#library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)
library(data.table)
library(here)
#library(readxl)
library(magrittr)
library(forcats)
library(GGally)
library(ggpubr)
## BioConductor packages
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Mmusculus.v79)
library(Matrix)
library(GenomicRanges)
library(splatter)
library(scater)
library(SummarizedExperiment)
library(SeuratDisk)
library(gespeR)
## Statistics / clustering
library(Hmisc) 
library(diffusr)
library(statmod)
library(LICORS)
library(Spectrum)
library(cluster)
library(mclust)
library(uwot)
library(NetworkToolbox)
library(zinbwave)



## Load All sideSEQ functions
source(here('R/side_seq/side_seq.R'))
source(here('R/side_seq/side_seq_iter.R'))
source(here('R/side_seq/side_seq_block.R'))
source(here('R/side_seq/random_sampling_methods.R'))
source(here('R/side_seq/computeDiffExprMat.R'))
source(here('R/side_seq/setSimilarityMeasure.R'))


source(here("R/side_seq/side_seq_ref_set.R"))
source(here("R/side_seq/weightedNearestNeighbors.R"))
source(here("R/side_seq/evaluate.R"))

`%p%` <- function(x,y) {paste0(x,y)}


## TODO: just have it source everything in the directory once cleaned up..






