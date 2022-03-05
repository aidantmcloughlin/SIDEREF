library(here)
source(here('R/libraries.R'))

## TODO: just for now
#source(here('R/genefishing.R'))
#source(here('R/hcclust.R'))

## Run analysis --------------------------------------------------------------- 

## Preprocessing:
pre_preproc_objs <- ls()

cat("Preprocessing data...")
source(here('R/PreProcessing.R'))

## clear preprocessing
rm(list = setdiff(ls(), pre_preproc_objs))


## Splatter:
pre_splat_objs <- ls()

cat("Running Splatter Sims...")
source(here('R/SplatterSimCreate.R'))
source(here('R/SplatterSIDEComputations.R'))
source(here('R/SplatterResMain.R'))
source(here('R/SplatterResOther.R'))

## clear splatter
rm(list = setdiff(ls(), pre_splat_objs))


## Snare-seq:

pre_snareatac_objs <- ls()
cat("Running SNARE ATAC results...")
source(here('R/SnareResATAC.R'))

## clear Snare ATAC
rm(list = setdiff(ls(), pre_snareatac_objs))


cat("Computing SNARE SIDEseq and SIDEref runs...")
source(here('R/SnareRNASIDEComputations.R'))

cat("Running SNARE RNA results...")
source(here('R/SnareResRNA.R'))


cat("Comparing SIDEseq to SIDEref...")
source(here("CompareSIDEseqSIDEREF.R"))

cat("Done!")


