### Create splatter simulation:
set.seed(128379)

library(here)
source(here("R/side_seq/master.R"))
source(here("R/side_seq/splatter_source/auxiliary_splatter_functions.R"))
source(here("R/side_seq/splatter_source/shared_genes_groups.R"))

params <- newSplatParams(batchCells = 250, nGenes = 1e4)
params@batchCells <- 250
params@nCells <- 250


## get multiplicative factor of 4
mean_mult <- 3

lnorm_mu1 <- log(mean_mult) -0.5 * 0.4^2

lnorm_mu2 <- log(mean_mult) -0.5 * 0.6^2

de_fac_scales <- c(0.4, 0.6, 0.4, 0.4, 0.4)
de_fac_locs <- c(lnorm_mu1, lnorm_mu2, lnorm_mu1, lnorm_mu1, lnorm_mu1)
de_downProbs <- c(0.5, 0.5, 0.5, 0.5, 0.5) 

de_vars <- c((exp(de_fac_scales^2) - 1) * exp(2*de_fac_locs + de_fac_scales^2))

## TODO: plug estimated dropout
params <- setParams(params, 
                    group.prob = c(0.05, 0.2, 0.45, 0.1, 0.2),
                    de.prob = c(0.02, 0.02, 0.025, 0.03, 0.04),
                    de.facLoc = de_fac_locs,
                    de.facScale = de_fac_scales,
                    de.downProb = de_downProbs)


sim1 <- splatSimulateGroups(params,
                            verbose = FALSE)



save(sim1, file = here("output/splatter_sideseq/splatter_sim1.RData"))



## TODO: SIM 2: rename below to sim 2, sim 3

###############################################################################
###   Simulation 3: Shared expression component within two subgroups
###############################################################################

set.seed(39817234)

params <- newSplatParams(batchCells = 250, nGenes = 1e4)
params@batchCells <- 250
params@nCells <- 250



## get multiplicative factor of 4
mean_mult <- 3



lnorm_mu1 <- log(mean_mult) -0.5 * 0.3^2

lnorm_mu2 <- log(mean_mult) -0.5 * 0.4^2

lnorm_mu3 <- log(2.5) -0.5*0.35^2

de_fac_scales <- c(0.3, 0.4, 0.3, 0.3, 0.3, 0.4)
de_fac_locs <- c(lnorm_mu1, lnorm_mu2, lnorm_mu1, lnorm_mu1, lnorm_mu1, lnorm_mu2)

de_vars <- c((exp(de_fac_scales^2) - 1) * exp(2*de_fac_locs + de_fac_scales^2))


### Group differential expression parameters
group_share_list <- list(c(1, 2, 3), c(4, 5))
group_share_probs <- c(0.02, 0.03)
group_share_fac_locs <- c(lnorm_mu3, lnorm_mu3)
group_share_fac_scales <- c(0.35, 0.35)
group_share_downProbs <- c(0.5, 0.5)

## mean
exp(lnorm_mu3 + 0.5*0.35^2)
#3 var
de_vars <- c((exp(0.35^2) - 1) * exp(2*lnorm_mu3 + 0.35^2))

## TODO: plug estimated dropout
params <- setParams(params, 
                    group.prob = c(0.1, 0.15, 0.3, 0.15, 0.1, 0.2),
                    de.prob = c(0.01, 0.01, 0.03, 0.03, 0.03, 0.04),
                    de.facLoc = de_fac_locs,
                    de.facScale = de_fac_scales,
                    de.downProb = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5))


sim3 <- splatSimulateAlt(params,
                         method = "groups",
                         group_share_list = group_share_list,
                         group_share_probs = group_share_probs,
                         group_share_facLocs = group_share_fac_locs,
                         group_share_facScales = group_share_fac_scales,
                         group_share_downProb = group_share_downProbs)


head(rowData(sim3), 15)




###############################################################################
###   Simulation 4: Augmented noise component within subgroups
###############################################################################

set.seed(523446)

params <- newSplatParams(batchCells = 250, nGenes = 1e4)
params@batchCells <- 250
params@nCells <- 250



## get multiplicative factor of 3
mean_mult <- 3



lnorm_mu1 <- log(mean_mult) -0.5 * 0.3^2

lnorm_mu2 <- log(mean_mult) -0.5 * 0.4^2



############################################################
### Augment noise inherent in all gene expression draws
############################################################




### up the noise levels
params@mean.rate <- 0.4
params@mean.shape <- 0.8



### set subgroup params

de_fac_scales <- c(0.3, 0.4, 0.3, 0.3, 0.3, 0.4)
de_fac_locs <- c(lnorm_mu1, lnorm_mu2, lnorm_mu1, lnorm_mu1, lnorm_mu1, lnorm_mu2)


### Group differential expression parameters
group_share_list <- list(c(1, 2, 3), c(4, 5))
group_share_probs <- c(0.02, 0.03)
group_share_fac_locs <- c(lnorm_mu3, lnorm_mu3)
group_share_fac_scales <- c(0.35, 0.35)
group_share_downProbs <- c(0.5, 0.5)


## TODO: plug estimated dropout
params <- setParams(params, 
                    group.prob = c(0.1, 0.15, 0.3, 0.15, 0.1, 0.2),
                    de.prob = c(0.01, 0.01, 0.03, 0.03, 0.03, 0.04),
                    de.facLoc = de_fac_locs,
                    de.facScale = de_fac_scales,
                    de.downProb = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5))


params@bcv.common <- 0.2

params@lib.scale <- 0.3

params@dropout.type = 'none'

sim4 <- splatSimulateAlt(params,
                         method = "groups",
                         group_share_list = group_share_list,
                         group_share_probs = group_share_probs,
                         group_share_facLocs = group_share_fac_locs,
                         group_share_facScales = group_share_fac_scales,
                         group_share_downProb = group_share_downProbs)


sim4 <- logNormCounts(sim4)
sim4_counts <- as.matrix(sim4@assays@data$logcounts)

## variance of True Counts
var(c(as.matrix(sim4@assays@data$TrueCounts)))

var(c(as.matrix(sim3@assays@data$TrueCounts)))

