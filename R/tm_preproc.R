## set seed
set.seed(10)

### load TM data and functions:
source(here("R/tm_helper_functions.R"))
### Preproc Functions
source(here("R/distance_compute_helpers.R"))

load(here("output/tab_muris_sc/preproc_data/tab_muris_full.RData"))

if(USE_DROPLET_ONLY) {
  droplet_cells <- which(grepl("droplet", meta_data_full$file_source))
  data_full <- data_full[, droplet_cells]
  meta_data_full <- meta_data_full[droplet_cells, ]
}

tm_preproc <-
  TabMurisProcSubset(full_data = data_full, 
                     meta_data = meta_data_full, 
                     var_features= VAR_FEATURES,
                     return_seurat = FALSE)


save(tm_preproc, 
     file = here("output/tab_muris_sc/tab_muris_full_scale_var_genes_" %p% 
                   VAR_FEATURES %p% 
                   ifelse(USE_DROPLET_ONLY, 
                          "_droplet_only", "") %p% 
                   ".RData"))
     
     
     
     