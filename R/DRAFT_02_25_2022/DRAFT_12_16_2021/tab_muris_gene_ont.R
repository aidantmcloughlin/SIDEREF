###############################################################################
### Tab Muris Gene Ontology Analysis
###############################################################################


## Load and Process Gene Ontology

LoadProcMouseGO <- function() {
  
  mouse_go <- fread(here("data/go_terms/all_go/gene_association.mgi"), skip=3)
  
  names(mouse_go) <- c("db_designation", "mgi_marker_access_id", "marker_symbol", 
                       "not_designation", "go_term_id", "mgi_ref_accession_id", 
                       "go_evidence_code", "inferred_from", "ontology",
                       "mouse_marker_name", "mouse_marker_synonyms", "mouse_marker_type", 
                       "taxon", "modification_date", "assigned_by", "annotation_extension", 
                       "gene_product")
  
  ## not needed with next DF
  # ontology_class_cw <- 
  #   data.frame(ontology_name = c("biological process", "molecular function", "cellular component"),
  #              ontology_code = c("P", "F", "C"))
  
  ## Ontology term crosswalk:
  ontology_term_cw <- fread(here("data/go_terms/all_go/go_terms.mgi.txt"), header=FALSE,
                            col.names=c("ontology_name", "go_term_id", "go_term"))
  
  mouse_go <- 
    mouse_go %>% 
    left_join(ontology_term_cw)
  
  
  mouse_go <- 
    mouse_go %>% 
    dplyr::filter(!go_term %in% c("molecular_function", "cellular_component", "biological_process"))
  
  
  
  num_terms <- length(unique(mouse_go$go_term_id))
  
  term_freq <-
    mouse_go %>% 
    group_by(go_term_id) %>% 
    summarise(count=n())
  
  ## plot
  term_freq %>% 
    dplyr::filter(count < 500) %>% 
    ggplot(aes(x=count)) + geom_histogram()
  
  ## TODO: can join terms with similar names but that's a lot of work to do properly.
  
  
  high_term_freq <-
    term_freq$go_term_id[which(term_freq$count > 30)]
  
  
  
  mouse_go_high_freq <- 
    mouse_go %>% 
    dplyr::filter(go_term_id %in% high_term_freq)
  
  mouse_go_high_freq_selected <- 
    mouse_go_high_freq %>% 
    dplyr::select(marker_symbol, go_term_id, ontology_name, go_term) %>% distinct()
  
  return(mouse_go_high_freq_selected) 
}

mouse_go <- LoadProcMouseGO()

###############################################################################
### Peruse Tab Muris for good candidate GO tems.
###############################################################################

## TODO with small sample data for now.

## load cell types
tab_muris_out_path = here("output/tab_muris_sc/")
load(here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct.RData")
#load(here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct_scale_Var_genes_3000.RData")




### Track Gene Expression per GO Term

GO_term_gene_expression <- function(express_mat, go_df, cell_types, tissues = NULL) {
  go_term_ids <- unique(go_df$go_term_id)
  
  ## DFs for which we will populate information:
  module_df <- distinct(go_df, go_term_id, go_term)
  
  celltype_module_df <- module_df %>% 
    mutate(to_join=TRUE) %>% 
    left_join(data.frame(celltype=unique(cell_types), 
                         to_join=TRUE)) %>% 
    dplyr::select(-to_join)
  
  ## Generate Placeholder columns
  module_df <- module_df %>% 
    mutate(n_genes = NA)
  
  for(i in seq_len(length(go_term_ids))) {
    genes = go_df$marker_symbol[which(go_df$go_term_id == go_term_ids[i])]
    
    module_df$n_genes[i] = length(genes)
    
  }
  
  return(module_df)
}


module_df <- GO_term_gene_expression(data_full_variable, mouse_go,
                                     cell_types = meta_data_full$cell_type,
                                     tissues = meta_data_full$tissue)



immune_response_genes <- 
  mouse_go$marker_symbol[mouse_go$go_term == "immune response"]

###############################################################################
### Fish subgroup specific gene modules:
###############################################################################


library(GeneFishing)


# probe_res <-
#   GeneFishing::probeFishability(liver_mat,
#                                 potential_bait = bait_genes_chol_metab,
#                                 n_rounds = 50,
#                                 ncores = 5)


# doParallel::stopImplicitCluster()


# gf_res <-
#   GeneFishing::geneFishing(liver_mat,
#                            bait_genes = res$best_bait,
#                            fishing_rounds = 50,
#                            n_probing_rounds = 50)


