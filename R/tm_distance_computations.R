###############################################################################
### Tabula Muris Distance Computations.

set.seed(1)

library(here)
print(here())
source(here('R/libraries.R'))
source(here('R/relative_group_dist_comps.R'))
source(here("R/distance_compute_helpers.R"))
source(here('R/tm_sample.R'))


############################# CONSTANTS ####################################
## Data subset Constants
KEEP_DROPLET <- TRUE
if(USE_DROPLET_ONLY) {
  KEEP_FACS <- FALSE  
} else{
  KEEP_FACS <- TRUE
}


file_name <- "tab_muris_full"

SAMPLE_CELL_TYPES = TRUE
SAMPLE_PCT = FALSE

## Whether it's needed to run the preproc pipeline:
RUN_DIST <- TRUE


tm_data_list <- 
  tmLoadCut(file_name = file_name, data_pct = data_pct, 
            1, 1,
            sample_size = TM_SAMP_SIZE,
            sample_cell_types = SAMPLE_CELL_TYPES,
            keep_droplet = KEEP_DROPLET, keep_facs = KEEP_FACS,
            sample_pct = SAMPLE_PCT)

## extract objects
for(n in names(tm_data_list)) {
  assign(n, tm_data_list[[n]])
}

rm(tm_data_list)


for(x in data_store_list[1:length(data_store_list)]) {
    set.seed(data_store_list[[1]][[4]])

    tm_clust_by_file <- list()

    data_full <- x[[2]]
    meta_data_full <- x[[3]]

    if(RUN_DIST) {
      print("running full distance computation")
      ## Overall-distance matrices:
      full_dist_res <- 
        TabMurisDistResPipeline(data_full, 
                                meta_data_full, 
                                filename = NULL, 
                                spectral_dims = TM_SPECTRAL_DIMS, 
                                n_cores = N_CORES,
                                R = R,
                                G = G,
                                PCs = TM_PCS,
                                var_features = VAR_FEATURES,
                                save_loc = here("output/tab_muris_sc/dist_res/samp_ct") %p% 
                                  SAMPLE_CELL_TYPES %p% "_samp_pct" %p% SAMPLE_PCT %p% 
                                  "_samp_size" %p% TM_SAMP_SIZE %p% 
                                  "_use_droplet_only" %p% USE_DROPLET_ONLY,
                                n_clust = N_CLUST,
                                n_pcs = N_PCS,
                                min_dist = MIN_DIST,
                                n_neighbors = N_NEIGHBORS)
      
      ### save data 
      save(meta_data_full,
           file = here("output/tab_muris_sc/dist_res/meta_data" %p%
                         "_samp_ct" %p% SAMPLE_CELL_TYPES %p% 
                         "_samp_pct" %p% SAMPLE_PCT %p%
                         "_samp_size" %p% TM_SAMP_SIZE %p% ".RData")) 

    }
    
}
