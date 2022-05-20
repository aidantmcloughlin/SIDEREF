library(here)
source(here('R/libraries.R'))
## Sim Figures  -----------------------------
## track current obj
pre_objs <- ls()

source(here("R/splatter_post_hoc.R"))
source(here("R/splatter_perf_metrics.R"))

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


## Supplemental Studies Figure Creation  -----------------------------

source(here("R/other_supplemental_figures.R"))


