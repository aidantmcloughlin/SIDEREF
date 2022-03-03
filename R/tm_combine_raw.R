###############################################################################
### Binding Disparate Tabula Muris Data Sets
###############################################################################

library(here)
source(here("R/libraries.R"))
source(here("R/tm_helper_functions.R"))

## Path variables
tab_muris_rel_path = here("data/tab_muris_sc/seurat/")
tab_muris_out_path = here("output/tab_muris_sc/preproc_data")


sampleBindTabMuris <- 
  function(in_path, out_path,
           do_pct_sample = TRUE, sample_pct = 0.2) {
    
    tab_muris_files = list.files(in_path)    
    
    init = TRUE
    for(f in tab_muris_files) {
      print(f)
      load(tab_muris_rel_path %p% "/" %p% f,
           verbose = TRUE)
      
      data_f <- tiss@raw.data
      genes_data <- rownames(tiss@data)
      cells_data <- colnames(tiss@data)
      
      data_f <- data_f[genes_data, cells_data]
      
      
      n_genes <- dim(data_f)[1]
      if("old_genes" %in% ls()) {
        if(abs(old_genes-n_genes) > 1e-6) {
          print("Gene difference between files")
        }
      }
      print(n_genes)
      old_genes <- n_genes
      
      
      meta_data <- 
        data.frame(cell_id = rownames(tiss@meta.data),
                   tissue = tiss@meta.data$tissue,
                   mouse_id = tiss@meta.data$mouse.id,
                   mouse_sex = tiss@meta.data$mouse.sex,
                   free_anotation = tiss@meta.data$free_annotation,
                   cell_type = tiss@meta.data$cell_ontology_class, 
                   cluster = tiss@meta.data$cluster.ids,
                   file_source = f)
      
      if(do_pct_sample) {
        cell_n = nrow(meta_data)
        cell_sample = sample(cell_n, size = ceiling(cell_n * sample_pct), 
                             replace=FALSE)
        
        data_f <- data_f[, cell_sample]
        meta_data  <- meta_data[cell_sample, ]
        
        
      }
      
      
      if(init) {
        init <- FALSE
        data_full <- data_f
        meta_data_full <- meta_data
      } else{
        data_full <-
          cbind(data_full,
                data_f)
        meta_data_full <- 
          bind_rows(meta_data_full, 
                    meta_data)
      }
      
    }
    
    meta_data_full <- prepTMGroupLabels(
      meta_data = meta_data_full, 
      include_tissue = TRUE)
    
    save(data_full, meta_data_full, 
         file = out_path %p% "/" %p% "tab_muris_full" %p% 
           ifelse(do_pct_sample, "_" %p% sample_pct %p% "pct.RData",
                  ".RData"))
    
    return(list(data_full, meta_data_full))
  }


tab_muris_data_list <-
  sampleBindTabMuris(in_path  = tab_muris_rel_path,
                     out_path = tab_muris_out_path,
                     do_pct_sample = FALSE)


rm(tab_muris_data_list)


