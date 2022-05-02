set.seed(1)
library(here)
## packages + SIDEREF functions + Constants. ==================================
source(here('R/libraries.R'))


## Simulation Analysis ========================================================

## track current obj
pre_objs <- ls()

cat("Running Splatter Sims...")
source(here('R/splatter_sim_create.R'))

## clear space
rm(list = setdiff(ls(), pre_objs))

## Sim Distance Computations  -----------------------------
## track current obj
pre_objs <- ls()


source(here("R/splatter_distance_computations.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))


## Sim Figures  -----------------------------
## track current obj
pre_objs <- ls()

source(here("R/splatter_post_hoc.R"))
source(here("R/splatter_perf_metrics.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))




## Real Data Analysis =========================================================

## Load and Bind Raw Files  -----------------------------
pre_objs <- ls()

cat("Loading Binding Raw Data...")
source(here("R/tm_combine_raw.R"))

## Preproc variable genes data  -----------------------------
source(here("R/tm_preproc.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))


## Run TM Distance Computations  -----------------------------
pre_objs <- ls()

cat("Real Data Distance Compute..,")
source(here("R/tm_distance_computations.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))


## Run TM Data Analysis and Figure Creation  -----------------------------
pre_objs <- ls()

source(here("R/tm_heatmaps.R"))


## clear space
rm(list = setdiff(ls(), pre_objs))

pre_objs <- ls()

source(here("R/tm_gsea.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))

pre_objs <- ls()

source(here("R/tm_perf_metrics.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))


## Supplemental Studies =======================================================
pre_objs <- ls()

cat("SIDEREF Relative Distance Convergence Study..")
source(here("R/SIDEREF_convergence.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))


pre_objs <- ls()

cat("SIDEREF Reference Set Size and Compute Time Study..")
source(here("R/snn_SIDEseq_SIDEREF.R"))

## clear space
rm(list = setdiff(ls(), pre_objs))

pre_objs <- ls()

source(here("R/compute_time_SIDEseq_SIDEREF.R"))



cat("Done!")

rm(list = ls())
