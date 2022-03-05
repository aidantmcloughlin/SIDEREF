### Real Data HClust Examples
`%p%` <- function(x,y) {paste0(x,y)}

## TODO: valid clustering pipeline for cell-line results? code dat up!

PRE_SCRIPT_OBJ <- ls()


###############################################################################
### SECTION 1: Load and Preprocess Data sets
###############################################################################


## Dataset 1: 68,000 PBMC Cells

main_loc <- here("data/pbmc68k/from10x")
data_folders <- list.files(main_loc)
data_folders <- data_folders[!grepl("tar\\.gz", data_folders)]
generic_data_loc <- "/filtered_matrices_mex/hg19/"

## Check gene lists
b_genes <- read.table(file = main_loc %p% "/" %p% data_folders[1] %p%
                        generic_data_loc %p% 
                        'genes.tsv', 
                      sep = '\t', header = FALSE)

for(i in seq_len(length(data_folders))) {
  new_genes <- read.table(file = main_loc %p% "/" %p% data_folders[i] %p%
                                       generic_data_loc %p% 
                                       'genes.tsv', 
                                     sep = '\t', header = FALSE)
  
  stopifnot(b_genes$V1 == new_genes$V1)
}
gene_names = new_genes
rm(b_genes, new_genes)

## all gene lists are the same.


## iterative load data and bind cols
for(i in seq_len(length(data_folders))) {
  print(i)
  curr_exp_matrix <- Matrix::readMM(main_loc %p% "/" %p% data_folders[i] %p%
                                      generic_data_loc %p% 
                                      "matrix.mtx")  
  
  curr_barcodes <- read.table(file = main_loc %p% "/" %p% data_folders[i] %p%
                                generic_data_loc %p% 
                                'barcodes.tsv', 
                              sep = '\t', header = FALSE)
  
  curr_celltypes <- rep(gsub("_filtered_gene_bc_matrices", "", data_folders[i]), nrow(curr_barcodes))
  
  if(i == 1) {
    exp_matrix     = curr_exp_matrix
    barcodes_full  = curr_barcodes 
    celltypes_full = curr_celltypes
  } else{
    exp_matrix = cbind(exp_matrix, curr_exp_matrix)
    barcodes_full = rbind(barcodes_full, curr_barcodes)
    celltypes_full = c(celltypes_full, curr_celltypes)
    
  }
}

## Store Full Object as Seurat

dimnames(exp_matrix) = list(gene_names$V1,barcodes_full$V1)
pbmc_seurat <- CreateSeuratObject(exp_matrix)

pbmc_seurat@meta.data$celltype = celltypes_full

## save Seurat object
save(pbmc_seurat, file = here("output/real_data/pbmc_seurat.RData"))

rm(exp_matrix, curr_exp_matrix, celltypes_full, curr_celltypes, 
   barcodes_full, curr_barcodes,
   gene_names)

pbmc_seurat <- FindVariableFeatures(pbmc_seurat, nfeatures = 2000)
pbmc_seurat <- NormalizeData(pbmc_seurat)
pbmc_seurat <- ScaleData(pbmc_seurat)

save(pbmc_seurat, file = here("output/real_data/pbmc_seurat_processed.RData"))

## can extract gene reduced matrix and sample?!

pbmc_var_features <- pbmc_seurat@assays$RNA@scale.data

# sample for sake of runtime
sample_size = 7500

sample_indices <- sample(ncol(pbmc_var_features), sample_size)

pbmc_var_features <- pbmc_var_features[,sample_indices]

## check cell distribution
sampled_barcodes <- colnames(pbmc_var_features)

table(pbmc_seurat@meta.data[sampled_barcodes,]$celltype)

table(pbmc_seurat@meta.data$celltype)

save(pbmc_var_features, file = here("output/real_data/pbmc_seurat_variable_sample.RData"))

## TODO: Dataset 2: 10X PBMC Cells



samp_cells <- colnames(full_res$dist_list$spearman_dist)
## Dataset 3: Pollen Glial Dataset 
#library(HGC)

#data(Pollen)

#tissue <- Pollen$Tissue
#cellline <- Pollen$CellLine


pollen_mRNA_seq <- readRDS(here("data/pollen/pollen.rds"))

pollen_counts <- pollen_mRNA_seq@assays$data$logcounts
#pollen_norm_counts <- pollen_mRNA_seq@assays$data$normcounts

pollen_cell_types <- pollen_mRNA_seq@colData$cell_type1
pollen_tissues <- pollen_mRNA_seq@colData$cell_type2

## Convert to Seurat to Find High Variability Features
pollen_Seurat <- CreateSeuratObject(pollen_counts)

pollen_Seurat <- FindVariableFeatures(pollen_Seurat, nfeatures = 3000)

head(pollen_Seurat@meta.data)

pollen_var_features <- 
  pollen_counts[rownames(pollen_counts) %in% 
                  pollen_Seurat@assays$RNA@var.features, ]

save(pollen_var_features, pollen_cell_types, pollen_tissues, 
     file = here("output/real_data/pollen_variable_genes_and_types.RData"))



POST_PROC_OBJ <- setdiff(ls(), PRE_SCRIPT_OBJ)

rm(POST_PROC_OBJ)


###############################################################################
### SECTION 2: Run Analytical Pipeline
###############################################################################


load(here("output/real_data/pbmc_dist_res.RData"))

load(here("output/real_data/pbmc_seurat_variable_sample.RData"))

load(here("output/real_data/pbmc_seurat.RData"))

pbmc_small = pbmc_seurat[, sample(ncol(pbmc_seurat),1000)]
#pbmc_small <- NormalizeData(pbmc_small)
#pbmc_small <- FindVariableFeatures(pbmc_small, nfeatures = 2000)
#pbmc_small <- ScaleData(pbmc_small)
#pbmc_small <- SCTransform(pbmc_small)

#pbmc_small_counts <- pbmc_small@assays$RNA@data[pbmc_small@assays$RNA@var.features, ]


#pbmc_small_counts <- pbmc_small@assays$SCT@scale.data
pbmc_small_counts <- pbmc_small@assays$RNA@counts
#pbmc_small_var = pbmc_var_features[, sample(ncol(pbmc_var_features), 1000)]


## remove low expression genes
pbmc_small_counts_filt <- pbmc_small_counts
pbmc_small_counts_filt <-
  pbmc_small_counts[which(rowSums(pbmc_small_counts) >= quantile(rowSums(pbmc_small_counts), probs=0.95)), ]


umap_res = 
  t(as.matrix(pbmc_small_counts_filt)) %>%
  uwot::umap(n_neighbors=15, pca=25)

library(Rtsne)
tsne_res = Rtsne(t(as.matrix(pbmc_small_counts_filt)))




cells_small = colnames(pbmc_small_counts_filt)


cell_types_small <- pbmc_seurat@meta.data[cells_small, ]$celltype

# umap_res <- 
#   uwot::umap(as.dist(full_res$dist_list$pca_dist_25),
#              n_neighbors = 5)
  

umap_res %>%
  data.frame() %>%
  mutate(celltype = cell_types_small) %>%
  ggplot(aes(x=X1, y=X2, color = as.factor(clust$cluster))) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "K-Means Cell Type")

umap_res %>%
  data.frame() %>%
  mutate(celltype = cell_types_small) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "Marked Cell Type")

tsne_res$Y %>%
  data.frame() %>%
  mutate(celltype = cell_types_small) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="TSNE1", y="TSNE2", color = "Marked Cell Type")

###############################################################################
## Kmeans cluster from UMAP results 
###############################################################################


diss <- stats::dist(umap_res)
k.max = 15
n.reps = 25
v <- rep(0, k.max)
for(i in 2:k.max) {
  print(i)
  for(r in 1:n.reps) {
    clust = kmeanspp(umap_res, i)
    v[i] <- v[i] + factoextra:::.get_ave_sil_width(diss, clust$cluster) / n.reps
  }
}

plot_df <- data.frame(clusters = as.factor(1:k.max), y = v, 
                      stringsAsFactors = TRUE)

p <- ggpubr::ggline(plot_df, x = "clusters", y = "y", 
                    group = 1, color = "steelblue", ylab = "Avg. Sil. Width", xlab = "Number of clusters k", 
                    main = "Optimal number of clusters (25 Runs)") +
  theme(title = element_text(size = 18))






###############################################################################
### Load Pollen Data
###############################################################################


load(here("output/real_data/pollen_variable_genes_and_types.RData"))

umap_res = 
  t(as.matrix(pollen_var_features)) %>%
  uwot::umap(n_neighbors=15, pca=15)


pollen_clust <- 
  kmeanspp(umap_res, 9)

umap_res %>%
  data.frame() %>%
  mutate(celltype = pollen_cell_types) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() +  
  labs(x="UMAP1", y="UMAP2", color = "Marked Cell Type") 


umap_res %>%
  data.frame() %>%
  mutate(celltype = pollen_tissues) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "Tissue Source")


umap_res %>%
  data.frame() %>%
  mutate(celltype = as.factor(pollen_clust$cluster)) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "K-Means Res")


test_hclust_res <- 
  pollen_hclust_res_df$pear_dist$hclust_res$clust_res_df

umap_res %>%
  data.frame() %>%
  mutate(celltype = as.factor(test_hclust_res$cluster_high)) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "HClust-Res")


source(here("R/subgroup_composition_gene_contrib.R"))
source(here("R/hcclust.R"))

pollen_dist_res <-
  runDimReds(pollen_var_features, pollen_cell_types, 
             r = 50, g=c(50,150,300), 
             n_clust = 9, n_cores = detectCores()-3)


### Hierarchical Clustering Results
pollen_hclust_res_df2 <-
  lapply(pollen_dist_res$dist_list,
         function(x) {
           res <- 
             RunHClustPipeline(
               clusters = as.numeric(pollen_clust$cluster), 
               subg_inds = sample(3, 301, replace=TRUE), 
               hclust_dissim_mat = x, 
               p_dissim_mat      = pollen_dist_res$dist_list$pca_dist_25)
         })


names(pollen_hclust_res_df) <- names(pollen_dist_res$dist_list)


lapply(pollen_hclust_res_df, function(x) x$adj_rand_index)
ggdendrogram(
  pollen_hclust_res_df$side_ref_g2_dist$hclust_res$hclust_res, 
  rotate = TRUE, theme_dendro = FALSE)

ggdendrogram(pollen_hclust_res_df$side_ref_g2_dist$hclust_res$hclust_res) +
  theme(ylim= c(100, 1000))

# HGC::HGC.PlotDendrogram(
#   pollen_hclust_res_df$side_ref_g2_dist$hclust_res$hclust_res, 
#   k = 9, plot.label = TRUE, labels = 1:9)

# sideres <- RunHClustPipeline(
#   clusters = as.numeric(splat_group_labels), 
#   subg_inds = splat_subgs_inds_w_singletons, 
#   hclust_dissim_mat = full_res$dist_list$side_ref_g150_dist, 
#   p_dissim_mat      = full_res$dist_list$pca_dist_25)
# 
# sideres$p


###############################################################################
### Load Other PBMC Data from Seurat Package

## How to access data:
#library(SeuratData)
#InstallData("ifnb")
#data("ifnb")
#ifnb.list <- SplitObject(ifnb, split.by = "stim")
#save(ifnb.list, file=here("output/ifnb_data.RData"))



## Load and preprocess
load(here("output/ifnb_data.RData"))


ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
})

ifnb <- ifnb.list$CTRL@assays$RNA@scale.data

ifnb_meta <- ifnb.list$CTRL@meta.data

### UMAP run
umap_res <-
  t(as.matrix(ifnb)) %>%
  uwot::umap(n_neighbors=15, pca=15)


umap_res %>%
  data.frame() %>%
  mutate(celltype = ifnb_meta$seurat_annotations) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "Cell Type")


ifnb_sample_ind <- sample(nrow(ifnb_meta), 1000)

ifnb_sample <- ifnb[,ifnb_sample_ind]

ifnb_meta_sample <- ifnb_meta[ifnb_sample_ind, ]


save(ifnb, ifnb_meta, ifnb_sample, ifnb_meta_sample, 
     file=here("output/ifnb_processed.RData"))



source(here("R/subgroup_composition_gene_contrib.R"))
source(here("R/hcclust.R"))

load(here("output/ifnb_processed.RData"))

ifnb_sample_dist_res <-
  runDimReds(ifnb_sample, ifnb_meta_sample$seurat_annotations, 
             r = 100, g=c(50,150,300), 
             n_clust = 15, n_cores = detectCores()-2)



### Hierarchical Clustering Results
ifnb_sample_hclust_res <-
  lapply(ifnb_sample_dist_res$dist_list,
         function(x) {
           res <- 
             RunHClustPipeline(
               clusters = as.numeric(as.factor(ifnb_meta_sample$seurat_annotations)), 
               hclust_dissim_mat = x, 
               subg_inds = sample(3, 301, replace=TRUE), 
               p_dissim_mat      = ifnb_sample_dist_res$dist_list$pca_dist_25)
         })


ifnb_meta_sample <- 
  ifnb_meta_sample %>% 
  mutate(subgs = 
           case_when(grepl("CD4", seurat_annotations) ~ "CD4 T",
                     grepl("CD8", seurat_annotations) ~ "CD8 T",
                     grepl("Mono", seurat_annotations) ~ "Monocyte",
                     seurat_annotations %in% c("B", "B Activated") ~ "B",
                     grepl("DC", seurat_annotations) ~ "Dendritic"))


## Not a clear set of annotations here.


# Get Res DF:
#names(pollen_hclust_res_df) <- names(pollen_dist_res$dist_list)


lapply(pollen_hclust_res_df, function(x) x$adj_rand_index)
ggdendrogram(
  pollen_hclust_res_df$side_ref_g2_dist$hclust_res$hclust_res, 
  rotate = TRUE, theme_dendro = FALSE)

ggdendrogram(pollen_hclust_res_df$side_ref_g2_dist$hclust_res$hclust_res) +
  theme(ylim= c(100, 1000))






### Allen Institute Human Brain Cells
dendro_info <- 
  readRDS(here("data/tab_muris_sc/10x_brain_human/human_dendrogram.rds"))




###############################################################################
### Tabula Muris Data
###############################################################################

`%p%` <- function(x,y) paste0(x,y)


sampleBindTabMuris <- 
  function(in_path, out_path,
           do_pct_sample = TRUE, sample_pct = 0.2) {
    
    tab_muris_files = list.files(in_path)    
    
    for(f in tab_muris_files) {
      print(f)
      load(tab_muris_rel_path %p% "/" %p% f,
           verbose = TRUE)
      
      scale_data <- tiss@scale.data
      n_genes <- dim(scale_data)[1]
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
        
        scale_data <- scale_data[, cell_sample]
        meta_data  <- meta_data[cell_sample, ]
        
        
      }
      
      
      rm(tiss)
      if(!"scale_data_full" %in% ls()) {
        scale_data_full <- scale_data
        meta_data_full <- meta_data
      } else{
        scale_data_full <-
          cbind(scale_data_full,
                scale_data)
        meta_data_full <- 
          bind_rows(meta_data_full, 
                    meta_data)
      }
      
    }
    
    save(scale_data_full, meta_data_full, 
         file = out_path %p% "/" %p% "tab_muris_full" %p% 
           ifelse(do_pct_sample, "_" %p% sample_pct %p% "pct.RData",
                  ".RData"))
  }


tab_muris_rel_path = here("data/tab_muris_sc/seurat/")
tab_muris_out_path = here("output/tab_muris_sc/")

## 20 Pct
sampleBindTabMuris(in_path  = tab_muris_rel_path, 
                   out_path = tab_muris_out_path, 
                   do_pct_sample = TRUE, sample_pct = 0.05)




## Looks Good
umap_res <- 
  t(as.matrix(scale_data_var_genes)) %>%
  uwot::umap(n_neighbors=15, pca=15)


umap_res %>%
  data.frame() %>%
  mutate(celltype = meta_data$cell_type) %>%
  ggplot(aes(x=X1, y=X2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", color = "Cell Type")



## TODO: dynamic sizing of geom_points would be a nice touch.



