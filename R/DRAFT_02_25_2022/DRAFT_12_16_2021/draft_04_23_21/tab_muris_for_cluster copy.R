###############################################################################
### Tabula Muris Data
###############################################################################

`%p%` <- function(x,y) paste0(x,y)

CORES=16

## SIDEREF GENE Hyperparam:
G = c(150, 300, 450)

## Number of Cells Hyperparam:
R = 150

SPECTRAL_DIMS = 10


## Path variables
tab_muris_rel_path = here("data/tab_muris_sc/seurat/")
tab_muris_out_path = here("output/tab_muris_sc/")


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
        cell_sample = sample(cell_n, size = ceiling(cell_n * sample_pct), replace=FALSE)
        
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
    
    save(data_full, meta_data_full, 
         file = out_path %p% "/" %p% "tab_muris_full" %p% 
           ifelse(do_pct_sample, "_" %p% sample_pct %p% "pct.RData",
                  ".RData"))
    
    return(list(data_full, meta_data))
  }


# ## 2 Pct
tab_muris_data_list <-
  sampleBindTabMuris(in_path  = tab_muris_rel_path,
                     out_path = tab_muris_out_path,
                     do_pct_sample = TRUE,
                     sample_pct = 0.02)
# ## 2.5 Pct
# tab_muris_data_list <-
#   sampleBindTabMuris(in_path  = tab_muris_rel_path,
#                      out_path = tab_muris_out_path,
#                      do_pct_sample = TRUE,
#                      sample_pct = 0.025)
# ## 3.5 Pct
# tab_muris_data_list <-
#   sampleBindTabMuris(in_path  = tab_muris_rel_path,
#                      out_path = tab_muris_out_path,
#                      do_pct_sample = TRUE,
#                      sample_pct = 0.035)
# ## 5 Pct
# tab_muris_data_list <-
#   sampleBindTabMuris(in_path  = tab_muris_rel_path,
#                      out_path = tab_muris_out_path,
#                      do_pct_sample = TRUE,
#                      sample_pct = 0.05)
# ## 20 Pct
# tab_muris_data_list <-
#   sampleBindTabMuris(in_path  = tab_muris_rel_path, 
#                      out_path = tab_muris_out_path, 
#                      do_pct_sample = TRUE,
#                      sample_pct=0.2)
# 
# ## Full
# tab_muris_data_list <-
#   sampleBindTabMuris(in_path  = tab_muris_rel_path, 
#                      out_path = tab_muris_out_path, 
#                      do_pct_sample = TRUE,
#                      sample_pct = 0.5)


## Normalization and Variable Genes:

tabMurisVarGenes <- function(in_file, out_path, n_var_genes) {
  load(in_file) 
  data_full <- CreateSeuratObject(data_full)
  
  data_full <- NormalizeData(data_full)
  data_full <- FindVariableFeatures(data_full, nfeatures = n_var_genes)
  data_full <- ScaleData(data_full)
  
  data_full_variable <- data_full@assays$RNA@scale.data
  
  save(data_full_variable, file = out_path %p% "_scale_var_genes_" %p% 
         n_var_genes %p% ".RData")  
  
}

###############################################################################
### With SeuratV3 Batch Effects correction
###############################################################################



tabMurisCorrectBatch <- function(in_file, out_path, n_var_genes,
                                 max_dim = 15) {
  load(in_file) 
  meta_data_full <-
    meta_data_full %>% 
    mutate(source = gsub("\\.Robj|droplet_|facs_", "", file_source))
  data_full <- CreateSeuratObject(data_full)
  
  data_full <- SetIdent(data_full, 
                        value = meta_data_full$source) 
  
  data_full <- SplitObject(data_full, split.by = "ident")
  
  # normalize and identify variable features for each dataset independently
  data_full <- lapply(X = data_full, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", 
                              nfeatures = n_var_genes)
  })
  
  
  ## Integrate the data sets
  features <- SelectIntegrationFeatures(object.list = data_full)
  
  data_anchors <- 
    FindIntegrationAnchors(object.list = data_full, 
                           anchor.features = features,
                           dims = 1:max_dim,
                           n.trees = 25)
  
  data_full <- IntegrateData(anchorset = data_anchors)
  
  ## Scale and Save
  DefaultAssay(data_full) <- "integrated"
  data_full <- ScaleData(data_full, verbose = FALSE)
  data_full_variable <- data_full@assays$integrated@scale.data
  
  save(data_full_variable, file = out_path %p% "_batch_corrected_" %p% 
         n_var_genes %p% "_genes" %p% ".RData")
}



# ## 0.05 Pct Sample
# tabMurisVarGenes(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct.RData",
#                  out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct",
#                  n_var_genes = 3000)

tabMurisCorrectBatch(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.025pct.RData",
                     out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.025pct",
                     n_var_genes = 2000,
                     max_dim = 5)
# 
# tabMurisVarGenes(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct.RData",
#                  out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct",
#                  n_var_genes = 5000)
# 
# 
# ## 0.20 Pct Sample
# tabMurisVarGenes(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.2pct.RData", 
#                  out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.2pct",
#                  n_var_genes = 3000)
# 
# tabMurisVarGenes(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.2pct.RData", 
#                  out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.2pct",
#                  n_var_genes = 5000)

print("Variable Features on the Full Data")
## Full Sample
tabMurisVarGenes(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.5pct.RData", 
                 out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.5pct",
                 n_var_genes = 3000)

tabMurisVarGenes(in_file = here(tab_muris_out_path) %p% "/tab_muris_full_0.5pct.RData", 
                 out_path = here(tab_muris_out_path) %p% "/tab_muris_full_0.5pct",
                 n_var_genes = 5000)

## Run Distance Computations
source(here("R/subgroup_composition_gene_contrib.R"))
source(here("R/hcclust.R"))

## load cell types
load(here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct.RData")

rm(data_full)

## scaled data
load(here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct_scale_Var_genes_3000.RData")

cat("0.05 pct Dist Res...")
dist_res_tab_muris_0.05pct_3000g <-
  runDimReds(data_full_variable, meta_data_full$cell_type, 
             r = R, g=G, 
             n_clust = 20, n_cores = CORES)

save(dist_res_tab_muris_0.05pct_3000g, 
     file = here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.05pct_3000g.RData")
## Just SIDEREF:
save(dist_res_tab_muris_0.05pct_3000g$dist_list[1:3], 
     file = here(tab_muris_out_path) %p% "/side_ref_tab_muris_dist_res_0.05pct_3000g.RData")


cat("0.2 pct Dist Res...")
load(here(tab_muris_out_path) %p% "/tab_muris_full_0.2pct_scale_Var_genes_3000.RData")

dist_res_tab_muris_0.2pct_3000g <-
  runDimReds(data_full_variable, meta_data_full$cell_type, 
             r = R, g=G, 
             n_clust = 20, n_cores = CORES)

save(dist_res_tab_muris_0.2pct_3000g, 
     file = here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.2pct_3000g.RData")

save(dist_res_tab_muris_0.2pct_3000g$dist_list[1:3], 
     file = here(tab_muris_out_path) %p% "/side_ref_tab_muris_dist_res_0.2pct_3000g.RData")


cat("0.5 pct Dist Res...")
load(here(tab_muris_out_path) %p% "/tab_muris_full_0.5pct_scale_Var_genes_3000.RData")

dist_res_tab_muris_0.5pct_3000g <-
  runDimReds(data_full_variable, meta_data_full$cell_type, 
             r = R, g=G, 
             n_clust = 20, n_cores = CORES)

save(dist_res_tab_muris_0.5pct_3000g, 
     file = here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.5pct_3000g.RData")

save(dist_res_tab_muris_0.5pct_3000g$dist_list[1:3], 
     file = here(tab_muris_out_path) %p% "/side_ref_tab_muris_dist_res_0.5pct_3000g.RData")





load(here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct.RData")

load(here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.05pct_3000g_lowg.RData",
     verbose = TRUE)

dist_res_tab_muris_0.05pct_3000g$dist_list <- 
  dist_res_tab_muris_0.05pct_3000g$dist_list[2:3]

dist_res_tab_muris_0.05pct_3000g_smallg <- dist_res_tab_muris_0.05pct_3000g


load(here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.05pct_3000g.RData",
     verbose = TRUE)

appended_names <- 
  c(names(dist_res_tab_muris_0.05pct_3000g$dist_list), "side_ref_g50_dist", "side_ref_g100_dist")

## append new ngene values
n_dists <- length(dist_res_tab_muris_0.05pct_3000g$dist_list)
dist_res_tab_muris_0.05pct_3000g$dist_list[[n_dists+1]] <- 
  dist_res_tab_muris_0.05pct_3000g_smallg$dist_list[[1]]

dist_res_tab_muris_0.05pct_3000g$dist_list[[n_dists+2]] <- 
  dist_res_tab_muris_0.05pct_3000g_smallg$dist_list[[2]]

names(dist_res_tab_muris_0.05pct_3000g$dist_list) <- appended_names

save(dist_res_tab_muris_0.05pct_3000g,
     file = here('output/tab_muris_sc/dist_res/' %p%
                   'full_tab_muris_dist_res_0.05pct_3000g.RData'))


dist_res_tab_muris_0.05pct_3000g$plot_list$plot_list_splat_g3 + 
  theme(legend.position = "none")

dist_res_tab_muris_0.05pct_3000g$plot_list$p_splat_pca_umap_25 + 
  theme(legend.position = "none")


## TODO: environment frame issues for these constants:
min_dist = 0.01
n_neighbors = 15

splatUMAPPlot(data  = as.dist(dist_res_tab_muris_0.05pct_3000g$dist_list$side_ref_g1_dist),
              title = "SIDEREF (PCA Embed Samp) " %p% R %p% " Cells " %p% G[1] %p% " Top Genes",
              celltypes = meta_data_full$tissue,
              plot_sil_score = FALSE)



### ARI wrt different Granularities -------------------------------------------

## Generally reloading data to ensure the meta_data object has proper labels
source(here("R/hcclust.R"))
load(here(tab_muris_out_path) %p% "/tab_muris_full_0.05pct.RData")
meta_data_full <-
  meta_data_full %>% 
  mutate(source = gsub("\\.Robj|droplet_|facs_", "", file_source))

## TODO: function for this will use adaptive naming from load() method.
## TODO: potentially filter NA cell types in previous data cleaning.


## HClust Results (wrt source file)
has_celltype <- which(!is.na(meta_data_full$cell_type))
tab_muris_0.05pct_3000g_hclust_res_to_source_df <-
  lapply(dist_res_tab_muris_0.05pct_3000g$dist_list,
         function(x) {
           res <- 
             RunHClustPipeline(
               clusters = as.numeric(factor(meta_data_full$cell_type[has_celltype])), 
               subg_inds = as.numeric(factor(meta_data_full$source[has_celltype])),
               hclust_dissim_mat = x[has_celltype, has_celltype], 
               p_dissim_mat      = 
                 dist_res_tab_muris_0.05pct_3000g$dist_list$pca_dist_25[has_celltype, has_celltype])
         })


names(tab_muris_0.05pct_3000g_hclust_res_to_source_df) <- 
  names(dist_res_tab_muris_0.05pct_3000g$dist_list)

lapply(tab_muris_0.05pct_3000g_hclust_res_to_source_df, 
       function(x) x$adj_rand_index)


## HClust Results (wrt tissue)
tab_muris_0.05pct_3000g_hclust_res_to_tissue_df <-
  lapply(dist_res_tab_muris_0.05pct_3000g$dist_list,
         function(x) {
           res <- 
             RunHClustPipeline(
               clusters = as.numeric(factor(meta_data_full$cell_type[has_celltype])), 
               subg_inds = as.numeric(factor(meta_data_full$tissue[has_celltype])),
               hclust_dissim_mat = x[has_celltype, has_celltype], 
               p_dissim_mat      = 
                 dist_res_tab_muris_0.05pct_3000g$dist_list$pca_dist_25[has_celltype, has_celltype])
         })


names(tab_muris_0.05pct_3000g_hclust_res_to_tissue_df) <- 
  names(dist_res_tab_muris_0.05pct_3000g$dist_list)

lapply(tab_muris_0.05pct_3000g_hclust_res_to_tissue_df, 
       function(x) x$adj_rand_index)


###############################################################################
###  Appending Spectral Results
###############################################################################


for(dist_name in names(dist_res_tab_muris_0.05pct_3000g$dist_list[
  !grepl("pca", names(dist_res_tab_muris_0.05pct_3000g$dist_list))])) {
  print(dist_name)
  n <- ncol(dist_res_tab_muris_0.05pct_3000g$dist_list[[1]])
  
  l <- length(dist_res_tab_muris_0.05pct_3000g$dist_list)
  
  dissim <- dist_res_tab_muris_0.05pct_3000g$dist_list[[dist_name]]
  
  dissim <- dist(spectralEmbed(dissim, ndim = SPECTRAL_DIMS))
  
  mat <- matrix(0, nrow = n, ncol = n)
  
  mat[lower.tri(mat)]  <- c(dissim)
  
  mat[upper.tri(mat)] <-
    t(mat)[upper.tri(mat)]
  
  dist_res_tab_muris_0.05pct_3000g$dist_list[[l+1]] <- mat
  names(dist_res_tab_muris_0.05pct_3000g$dist_list)[l+1] <- 
    dist_name %p% "_spectral_d" %p% SPECTRAL_DIMS %p% "_dist"
}




## HClust Results (low cell type to high-level cell type)
has_lowtype_and_hightype <- 
  which(!is.na(meta_data_full$free_anotation) & !is.na(meta_data_full$cell_type))

tab_muris_0.05pct_3000g_hclust_res_to_celltype_df <-
  lapply(dist_res_tab_muris_0.05pct_3000g$dist_list,
         function(x) {
           res <- 
             RunHClustPipeline(
               clusters = as.numeric(factor(meta_data_full$free_anotation[has_lowtype_and_hightype])), 
               subg_inds = as.numeric(factor(meta_data_full$cell_type[has_lowtype_and_hightype])),
               hclust_dissim_mat = x[has_lowtype_and_hightype, has_lowtype_and_hightype], 
               p_dissim_mat      = 
                 dist_res_tab_muris_0.05pct_3000g$dist_list$pca_dist_25[has_lowtype_and_hightype, has_lowtype_and_hightype])
         })


names(tab_muris_0.05pct_3000g_hclust_res_to_celltype_df) <- 
  names(dist_res_tab_muris_0.05pct_3000g$dist_list)

lapply(tab_muris_0.05pct_3000g_hclust_res_to_celltype_df, 
       function(x) x$adj_rand_index)






###############################################################################
### TODO: Hclust animation plot  ##############################################
###############################################################################

load(here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.05pct_3000g_lowg.RData",
     verbose = TRUE)

sideref_50gdist <- dist_res_tab_muris_0.05pct_3000g$dist_list$side_ref_g2_dist
sideref_50gumap <- dist_res_tab_muris_0.05pct_3000g$umap_list$side_ref_g2_umap
sideref_50gplot <- dist_res_tab_muris_0.05pct_3000g$plot_list$plot_list_splat_g2
rm(dist_res_tab_muris_0.05pct_3000g)


load(here(tab_muris_out_path) %p% "/tab_muris_dist_res_0.05pct_3000g.RData",
     verbose = TRUE)



PlotHClustResTabMuris <- function(hclust_res_df,
                                  dist_mat_for_p,
                                  size = 0.5,
                                  n_neighbors = 15,
                                  min_dist = 0.01) {
  
  ## get HClust results
  
  ## Get Umap Embedding
  umap_res <-
    uwot::umap(as.dist(dist_mat_for_p),
               n_neighbors = n_neighbors,
               min_dist = min_dist) %>% data.frame()
  
  
  if(length(unique(hclust_res_df$cell_group_low)) == nrow(hclust_res_df)) {
    p <- 
      ggarrange( 
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cell_group_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("Actual Subgroups") +
          theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cluster_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("HClust Subgroups") +
          theme(title = element_text(size=8)),
        ncol = 2
      )
  } else{
    p <- 
      ggarrange( 
        # umap_res %>% 
        #   mutate(celltype = factor(hclust_res_df$cell_group_low)) %>%
        #   ggplot(aes(x = X1, y = X2, col = celltype)) + 
        #   geom_point(size = size) + 
        #   theme_bw() +
        #   labs(x="UMAP1", y="UMAP2") +
        #   theme(legend.position = "none") +
        #   ggtitle("Actual Low-Level Groups") +
        #   theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cell_group_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("Actual High-Level Groups") +
          theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cluster_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("HClust Subgroups") +
          theme(title = element_text(size=8)),
        ncol = 2
      )
  }
  
  return(p)
  
}








# dist_res_tab_muris_0.05pct_3000g$plot_list$plot_list_splat_g2 + 
#   theme(legend.position = "none")


### TODO: few gene pathways of interest.

### TODO: screening interesting pathways --- histogram of mean pathway expression per cell type.






