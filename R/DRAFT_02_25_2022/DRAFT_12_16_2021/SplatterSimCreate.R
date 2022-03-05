### Create splatter simulation:
n_data_reps = 10

source(here("R/splatter_source/AuxiliarySplatterFunctions.R"))
source(here("R/splatter_source/splatSimulateAlt.R"))

## Primary simulation generator function:
generateSplatterSim <- function(
  cells = 250,
  n_genes = 1e4,
  group_probs,
  de_probs,
  mult_factors,
  fac_scales,
  group_share_params = NULL,
  save_filename = NULL,
  auxiliary_base_params = NULL,
  seed = 1) {
  
  set.seed(seed)
  
  ## constants
  n_groups <- length(group_probs)
  n_subgroups <- length(group_share_params[[1]])
  
  
  ##### Stop Checks -----------------------------------------------------------
  
  ## specifying the same number of groups for each parameter
  stopifnot(length(unique(
    c(length(group_probs),
      length(de_probs),
      length(mult_factors),
      length(fac_scales)))) == 1)
  
  ## names in group shared DE params
  stopifnot(setdiff(names(group_share_params),
                    c("shared_de_list", "shared_de_probs", 
                      "shared_de_mult_factors", "shared_de_fac_scales")) == "")
  
  ## appropriate length of group shared DE inputs if given
  if(!is.null(group_share_params)) {
    stopifnot(length(unique(c(
      length(group_share_params$shared_de_list),
      length(group_share_params$shared_de_probs),
      length(group_share_params$shared_de_mult_factors),
      length(group_share_params$shared_de_fac_scales)))) == 1)
    
    stopifnot(length(group_share_params) <= n_groups)
  }
  
  
  
  ## base splatter simulation
  base_params <- newSplatParams(batchCells = cells, nGenes = n_genes)
  base_params@nCells <- cells
  
  ## mean parameters for lognormal
  fac_locs <- 
    log(mult_factors) - 
    0.5 * fac_scales^2
  
  
  ## set base group DE parameters
  params <- splatter::setParams(base_params, 
                                group.prob  = group_probs,
                                de.prob     = de_probs,
                                de.facLoc   = fac_locs,
                                de.facScale = fac_scales,
                                de.downProb = rep(0.5, n_groups))
  
  
  ## Run Baseline Simulation before Noise adjustments
  
  if(is.null(group_share_params)) {

    sim_base <- splatSimulate(params)
  } else{
    shared_fac_locs <-
      log(group_share_params$shared_de_mult_factors) -
      0.5 * group_share_params$shared_de_fac_scales^2

    sim_base <- splatSimulateAlt(
      params,
      method = "groups",
      group_share_list      = group_share_params$shared_de_list,
      group_share_probs     = group_share_params$shared_de_probs,
      group_share_facLocs   = shared_fac_locs,
      group_share_facScales = group_share_params$shared_de_fac_scales,
      group_share_downProb  = rep(0.5, n_subgroups))
  }
  
  ## compute count variance
  var_sim_base = var(c(as.matrix(sim_base@assays@data$TrueCounts)))
  rm(sim_base)
  ## update base rate shape parameters if provided
  if(!is.null(auxiliary_base_params)) {
    params@mean.rate <- 
      ifelse(is.null(auxiliary_base_params$mean_rate),
             params@mean.rate, auxiliary_base_params$mean_rate)
    
    params@mean.shape <- 
      ifelse(is.null(auxiliary_base_params$mean_shape),
             params@mean.shape, auxiliary_base_params$mean_shape)
    
    params@bcv.common <- 
      ifelse(is.null(auxiliary_base_params$bcv_common),
             params@bcv.common, auxiliary_base_params$bcv_common)
    
    params@lib.scale <- 
      ifelse(is.null(auxiliary_base_params$lib_scale),
             params@lib.scale, auxiliary_base_params$lib_scale)
  }
  
  ## Run Final Simulation after Noise adjustments
  if(is.null(group_share_params)) {
    sim <- splatSimulate(params)
  } else{
    sim <- splatSimulateAlt(
      params,
      method = "groups",
      group_share_list      = group_share_params$shared_de_list,
      group_share_probs     = group_share_params$shared_de_probs,
      group_share_facLocs   = shared_fac_locs,
      group_share_facScales = group_share_params$shared_de_fac_scales,
      group_share_downProb  = rep(0.5, n_subgroups))
  }
  
  ## compute count variance
  var_sim = var(c(as.matrix(sim@assays@data$TrueCounts)))
  
  var_pct_change = (var_sim - var_sim_base) / var_sim_base
  
  ## compute DE variances for table
  de_vars <- 
    c((exp(fac_scales^2) - 1) * 
        exp(2*fac_locs + 
              fac_scales^2))
  
  if(is.null(group_share_params)) {
    ## Create table of simulation metadata.
    splat_params_table <-
      data.frame(cell_group     = paste0("group", seq_len(n_groups)),
                 de_prob        = de_probs,
                 de_factor_mean = mult_factors,
                 de_factor_var  = de_vars)
  } else{
    group_de_vars <-
      c((exp(group_share_params$shared_de_fac_scales^2) - 1) * 
          exp(2*shared_fac_locs + 
                group_share_params$shared_de_fac_scales^2))
    
    ## Create table of simulation metadata.
    splat_params_table <-
      data.frame(cell_group     = 
                   c(paste0("group", seq_len(n_groups)),
                     paste0("groups", 
                            sapply(group_share_params$shared_de_list, 
                                   function(x) paste(x, collapse="_")))),
                 de_prob        = c(de_probs, group_share_params$shared_de_probs),
                 de_factor_mean = c(mult_factors, group_share_params$shared_de_mult_factors),
                 de_factor_var  = c(de_vars, group_de_vars))
    
  }
  
  
  splat_res = list(sim = sim,
                   splat_params_table = splat_params_table,
                   var_pct_change = var_pct_change)
  
  if(!is.null(save_filename)) {
    save(splat_res, file = save_filename)
  }
  
  return(splat_res)
  
}



###############################################################################
##  Simulation 1: Fully isolated differentially expressed genes btw groups.
###############################################################################


generateSplatterSim(cells              = 250,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = c(0.05, 0.2, 0.45, 0.1, 0.2),
                    de_probs           = c(0.02, 0.02, 0.025, 0.03, 0.04),
                    mult_factors       = rep(3, 5),
                    fac_scales         = c(0.4, 0.35, 0.4, 0.4, 0.4),
                    ## group share params:
                    group_share_params = NULL,
                    save_filename = here("output/splatter_sim1.RData"),
                    auxiliary_base_params = NULL)



###############################################################################
###   Simulation 2: Shared expression component within two subgroups
###############################################################################


generateSplatterSim(cells              = 250,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = c(0.1, 0.15, 0.3, 0.15, 0.1, 0.2),
                    de_probs           = c(0.01, 0.01, 0.03, 0.03, 0.03, 0.04),
                    mult_factors       = rep(3, 6),
                    fac_scales         = c(0.3, 0.4, 0.3, 0.3, 0.3, 0.4),
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         = list(c(1, 2, 3), c(4, 5)),
                           shared_de_probs        = c(0.02, 0.03),
                           shared_de_mult_factors = c(3, 3),
                           shared_de_fac_scales   = c(0.3, 0.3)),
                    save_filename = here("output/splatter_sim2.RData"),
                    auxiliary_base_params = NULL)




###############################################################################
###   Simulation 2: Shared expression component within two subgroups
###############################################################################


generateSplatterSim(cells              = 250,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = c(0.1, 0.15, 0.3, 0.15, 0.1, 0.2),
                    de_probs           = c(0.01, 0.01, 0.03, 0.03, 0.03, 0.04),
                    mult_factors       = rep(3, 6),
                    fac_scales         = c(0.3, 0.4, 0.3, 0.3, 0.3, 0.4),
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         = list(c(1, 2, 3), c(4, 5)),
                           shared_de_probs        = c(0.02, 0.03),
                           shared_de_mult_factors = c(3, 3),
                           shared_de_fac_scales   = c(0.3, 0.3)),
                    save_filename = here("output/splatter_sim2.RData"),
                    auxiliary_base_params = NULL)


###############################################################################
###   Simulation 3: Augmented noise component within subgroups
###############################################################################

### Augment noise inherent in all gene expression draws

generateSplatterSim(cells              = 250,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = c(0.1, 0.15, 0.3, 0.15, 0.1, 0.2),
                    de_probs           = c(0.01, 0.01, 0.03, 0.03, 0.03, 0.04),
                    mult_factors       = rep(3, 6),
                    fac_scales         = c(0.3, 0.4, 0.3, 0.3, 0.3, 0.4),
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         = list(c(1, 2, 3), c(4, 5)),
                           shared_de_probs        = c(0.02, 0.03),
                           shared_de_mult_factors = c(3, 3),
                           shared_de_fac_scales   = c(0.3, 0.3)),
                    save_filename = here("output/splatter_sim3.RData"),
                    auxiliary_base_params = 
                      list(mean_rate  = 0.4,
                           mean_shape = 0.8,
                           bcv_common = 0.2,
                           lib_scale  = 0.3)
)




###############################################################################
###   Simulation 4: Shared expression component several subgroups
###############################################################################

## SHAREseq skin data shows following potential subgroup layout


### subgroup parameters ---------------------------------------------

## Indexed subgroup list
# skin_shared_de_list <-
#   c(1, 2, 2, 3, 3, 4, 4, 5, 
#     6, 6, 6, 6, 7, 7, 7, 7, 
#     8, 8, 8, 8, 8, 8)

skin_shared_de_list <-
  list(c(2, 3), c(4, 5), c(6, 7), 
       seq(9, 12, 1), 
       seq(13, 16, 1), 
       seq(17, 22, 1))
#length(skin_shared_de_list) == 6
## TODO: change params
skin_shared_de_probs   = c(0.03, 0.02, 0.01, 0.02, 0.02, 0.04)
skin_shared_mult_facs  = rep(3, 6)
skin_shared_fac_scales = c(rep(0.3, 4), 0.4, 0.3)


### group parameters ------------------------------------------------
skin_group_probs <-
  c(1, 
    1, 1, 
    2, 1, 
    0.5, 0.5, 
    0.5, 
    3, 3, 3, 3, 
    1, 1, 1, 1, 
    2.5, 2.5, 1, 3, 1, 1)


skin_de_probs <- 
  c(0.04, 
    0, 0.02, 
    0, 0.01,
    0.01, 0.02,
    0.02,
    0.0, 0.005, 0.01, 0.02,
    0.0, 0.005, 0.01, 0.02,
    0.0, 0.02, 0.03, 0.03, 0.035, 0.035)

skin_de_mult_facs <- rep(3, 22)
skin_de_fac_scales <-
  c(rep(0.3, 12), rep(0.45, 4), rep(0.3, 6))
  
skin_group_probs <- skin_group_probs / sum(skin_group_probs)

generateSplatterSim(cells              = 1000,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = skin_group_probs,
                    de_probs           = skin_de_probs,
                    mult_factors       = skin_de_mult_facs,
                    fac_scales         = skin_de_fac_scales,
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         = skin_shared_de_list,
                           shared_de_probs        = skin_shared_de_probs,
                           shared_de_mult_factors = skin_shared_mult_facs,
                           shared_de_fac_scales   = skin_shared_fac_scales),
                    save_filename = here("output/splatter_sim_skin.RData"),
                    auxiliary_base_params = 
                      list(mean_rate  = 0.35,
                           mean_shape = 0.65,
                           lib_scale  = 0.25)
)



###############################################################################
###   Simulation 5: More individual group DE
###############################################################################

## SHAREseq skin data shows following potential subgroup layout


### subgroup parameters ---------------------------------------------

## Indexed subgroup list
# skin_shared_de_list <-
#   c(1, 2, 2, 3, 3, 4, 4, 5, 
#     6, 6, 6, 6, 7, 7, 7, 7, 
#     8, 8, 8, 8, 8, 8)

skin_shared_de_list <-
  list(c(2, 3), c(4, 5), c(6, 7), 
       seq(9, 12, 1), 
       seq(13, 16, 1), 
       seq(17, 22, 1))

skin_shared_de_probs   = c(0.01, # 2, 3 
                           0.03, # 4, 5
                           0.03, # 6, 7
                           0.01, # 9, 10, 11, 12
                           0.01, # 13, 14, 15, 16
                           0.01) # 17, 18, 19, 20, 21, 22

skin_shared_mult_facs  = rep(3, 6)
skin_shared_fac_scales = c(rep(0.3, 4), 0.4, 0.4)


### group parameters ------------------------------------------------
skin_group_probs <-
  c(1, 
    1, 1,   # 2, 3
    2, 1,   # 4, 5
    0.5, 0.5, # 6, 7
    0.5, 
    3, 3, 3, 3, # 9, 10, 11, 12
    1, 1, 1, 1, # 13, 14, 15, 16
    2.5, 2.5, 1, 3, 1, 1) # 17, 18, 19, 20, 21, 22


skin_de_probs <- 
  c(0.04, 
    0.03, 0.03, # high ind (low shared)
    0.03, 0.03, # high ind (high shared)
    0.00, 0.01, # low ind (high shared)
    0.02,
    0.0, 0.005, 0.01, 0.02, # low ind (low shared)
    0.0, 0.005, 0.01, 0.02, # low ind (low shared, high variance)
    0.0, 0.02, 0.03, 0.03, 0.035, 0.035) # high ind (low shared, high variance)

skin_de_mult_facs <- rep(3, 22)
skin_de_fac_scales <-
  c(rep(0.3, 12), rep(0.4, 4), rep(0.4, 6))

skin_group_probs <- skin_group_probs / sum(skin_group_probs)

generateSplatterSim(cells              = 1000,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = skin_group_probs,
                    de_probs           = skin_de_probs,
                    mult_factors       = skin_de_mult_facs,
                    fac_scales         = skin_de_fac_scales,
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         = skin_shared_de_list,
                           shared_de_probs        = skin_shared_de_probs,
                           shared_de_mult_factors = skin_shared_mult_facs,
                           shared_de_fac_scales   = skin_shared_fac_scales),
                    save_filename = here("output/splatter_sim_skin_2_highvar.RData"),
                    auxiliary_base_params = 
                      list(mean_rate  = 0.35,
                           mean_shape = 0.70,
                           lib_scale  = 0.275)
)


generateSplatterSim(cells              = 1000,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        = skin_group_probs,
                    de_probs           = skin_de_probs,
                    mult_factors       = skin_de_mult_facs,
                    fac_scales         = skin_de_fac_scales,
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         = skin_shared_de_list,
                           shared_de_probs        = skin_shared_de_probs,
                           shared_de_mult_factors = skin_shared_mult_facs,
                           shared_de_fac_scales   = skin_shared_fac_scales),
                    save_filename = here("output/splatter_sim_skin_2_lowvar.RData"),
                    auxiliary_base_params = NULL)










###############################################################################
###   Simulation 6: More individual group DE
###############################################################################

## SHAREseq skin data shows following potential subgroup layout


### subgroup parameters ---------------------------------------------

## Indexed subgroup list
# skin_shared_de_list <-
#   c(1, 2, 2, 3, 3, 4, 4, 5, 
#     6, 6, 6, 6, 7, 7, 7, 7, 
#     8, 8, 8, 8, 8, 8)

skin_shared_de_list <-
  list(c(2, 3, 4), 
       c(5, 6, 7), 
       c(8, 9, 10), 
       c(12, 13, 14), 
       c(15, 16, 17), 
       c(18, 19, 20))

skin_shared_de_probs   = c(0.01, # 2, 3, 4
                           0.03, # 5, 6, 7
                           0.03, # 8, 9, 10
                           0.01, # 12, 13, 14
                           0.01, # 15, 16, 17
                           0.01) # 18, 19, 20

skin_shared_mult_facs  = rep(3, 6)
skin_shared_fac_scales = c(rep(0.3, 4), 0.4, 0.4)


### group parameters ------------------------------------------------
skin_group_probs <-
  c(1, 
    1, 1, 1, # 2, 3, 4
    1, 1, 1, # 5, 6, 7
    1, 1, 1, # 8, 9, 10
    1, 
    1, 1, 1, # 12, 13, 14
    1, 1, 1, # 15, 16, 17
    1, 1, 1) # 18, 19, 20


skin_de_probs <- 
  c(0.04, 
    0.02, 0.03, 0.04, # high ind (low shared)
    0.02, 0.03, 0.04, # high ind (high shared)
    0.005, 0.01, 0.015, # low ind (high shared)
    0.02,
    0.005, 0.01, 0.015, # low ind (low shared)
    0.005, 0.01, 0.015, # low ind (low shared, high variance)
    0.02, 0.03, 0.04) # high ind (low shared, high variance)

skin_de_mult_facs <- rep(3, 20)
skin_de_fac_scales <-
  c(rep(0.3, 14), rep(0.4, 3), rep(0.4, 3))

## floating point error with splatter package current version.
skin_group_probs <- round(skin_group_probs / sum(skin_group_probs), 24)


for(i in 1:n_data_reps) {
  
  generateSplatterSim(cells              = 1000,
                      n_genes            = 1e4,
                      ## ind group params:
                      group_probs        = skin_group_probs,
                      de_probs           = skin_de_probs,
                      mult_factors       = skin_de_mult_facs,
                      fac_scales         = skin_de_fac_scales,
                      ## group share params:
                      group_share_params = 
                        list(shared_de_list         = skin_shared_de_list,
                             shared_de_probs        = skin_shared_de_probs,
                             shared_de_mult_factors = skin_shared_mult_facs,
                             shared_de_fac_scales   = skin_shared_fac_scales),
                      save_filename = 
                        here("output/splatter_sim/splatter_sim_skin_eq_grps_highvar_rep" %p% i %p%
                               ".RData"),
                      auxiliary_base_params = 
                        list(mean_rate  = 0.35,
                             mean_shape = 0.70,
                             lib_scale  = 0.275)
  )
  
  
  generateSplatterSim(cells              = 1000,
                      n_genes            = 1e4,
                      ## ind group params:
                      group_probs        = skin_group_probs,
                      de_probs           = skin_de_probs,
                      mult_factors       = skin_de_mult_facs,
                      fac_scales         = skin_de_fac_scales,
                      ## group share params:
                      group_share_params = 
                        list(shared_de_list         = skin_shared_de_list,
                             shared_de_probs        = skin_shared_de_probs,
                             shared_de_mult_factors = skin_shared_mult_facs,
                             shared_de_fac_scales   = skin_shared_fac_scales),
                      save_filename = 
                        here("output/splatter_sim/splatter_sim_skin_eq_grps_lowvar_rep" %p% i %p% 
                               ".RData"),
                      auxiliary_base_params = NULL)
  
  
  
}
