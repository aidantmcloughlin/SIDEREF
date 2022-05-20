set.seed(1)
library(here)
## packages + SIDEREF functions + Constants. ==================================
source(here('R/libraries.R'))


##### COMPUTATIONALLY EXPENSIVE DISTANCE SCRIPTS
source(here("R/heavy_computations.R"))


##### POST-HOC ANALYSIS
##### (Can be run after downloading output files from: 
#####   )
source(here("R/post_hoc_analysis.R"))


cat("Done!")

rm(list = ls())
