### Create splatter simulation:
set.seed(1)

source(here("R/Splatter_scripts/AuxiliarySplatterFunctions.R"))
source(here("R/Splatter_scripts/splatSimulateAlt.R"))

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




### subgroup parameters ---------------------------------------------

## Indexed concat-group list
 shared_de_list <-
  list(c(3, 4, 5), 
       c(6, 7, 8), 
       c(9, 10, 11), 
       c(12, 13, 14), 
       c(15, 16, 17), 
       c(18, 19, 20))

shared_de_probs   = c(0.01, # 3, 4, 5
                      0.03, # 6, 7, 8
                      0.03, # 9, 10, 11
                      0.01, # 12, 13, 14
                      0.01, # 15, 16, 17
                      0.01) # 18, 19, 20

 shared_mult_facs  = rep(3, 6)
 shared_fac_scales = c(rep(0.3, 4), 0.4, 0.4)


### group parameters ------------------------------------------------
 group_probs <-
  c(1, 
    1, 1, 1, # 2, 3, 4
    1, 1, 1, # 5, 6, 7
    1, 1, 1, # 8, 9, 10
    1, 
    1, 1, 1, # 12, 13, 14
    1, 1, 1, # 15, 16, 17
    1, 1, 1) # 18, 19, 20


 de_probs <- 
  c(0.04, 
    0.02,
    0.02, 0.03, 0.04, # high ind (low shared)
    0.02, 0.03, 0.04, # high ind (high shared)
    0.005, 0.01, 0.015, # low ind (high shared)
    0.005, 0.01, 0.015, # low ind (low shared)
    0.005, 0.01, 0.015, # low ind (low shared, high variance)
    0.02, 0.03, 0.04) # high ind (low shared, high variance)

 de_mult_facs <- rep(3, 20)
 de_fac_scales <-
  c(rep(0.3, 14), rep(0.4, 3), rep(0.4, 3))

## floating point error with splatter package current version.
 group_probs <- round(group_probs / sum( group_probs), 24)


generateSplatterSim(cells              = 1000,
                    n_genes            = 1e4,
                    ## ind group params:
                    group_probs        =  group_probs,
                    de_probs           =  de_probs,
                    mult_factors       =  de_mult_facs,
                    fac_scales         =  de_fac_scales,
                    ## group share params:
                    group_share_params = 
                      list(shared_de_list         =  shared_de_list,
                           shared_de_probs        =  shared_de_probs,
                           shared_de_mult_factors =  shared_mult_facs,
                           shared_de_fac_scales   =  shared_fac_scales),
                    save_filename = 
                      here("output/splatter_sim/splatter_sim.RData"),
                    auxiliary_base_params = NULL)


