### Real Data Examples --------------------------------------------------------


## Set environment
library(here)
source(here("R/side_seq/master.R"))

## Set constants
n_neighbors <- 15L
min_dist <- 0.01

### 1. SNAREseq --------------------------------------------------------------- 

## load
load(here("output/snareseq_final.RData"))

## extract data of interest

## ATAC
DefaultAssay(snare) <- "ATAC"

## get variable chromatin regions
load(here("output/snareseq_atac_peaks_variable.RData"))
atac_peaks_variable <- snare_atac_variable %>% as.matrix()


## get variable chromatin genes 
atac_genes_variable <- 
  snare@assays$ATAC_gene@scale.data

## RNAseq
gex <- snare@assays$RNA@scale.data %>% 
  as.matrix()

#########################################
## A: SVD on Variable Peaks 
#########################################

### SVD
atac_peaks_svd <- RunSVD(snare_atac_variable)

## standard embedding
umap_atac_peaks_svd <-
  uwot::umap(data.frame(atac_peaks_svd@cell.embeddings)[,2:30],
             n_neighbors = n_neighbors,
             min_dist = min_dist)


p_umap_atac_pca <-
  data.frame(umap_atac_peaks_svd) %>% 
  ggplot(aes(x = X1, y = X2)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 SVD approx Data Matrix") + 
  theme(legend.position = "none")

p_umap_atac_pca



#########################################
## B: DEFAULT
#########################################

snare <- RunUMAP(snare, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')

p2 <- DimPlot(snare, reduction = 'umap.atac') + NoLegend() + ggtitle("ATAC UMAP")

p2

#########################################
## B2: PCA on full mat
#########################################

DefaultAssay(snare) <- "ATAC"

snare <- FindVariableFeatures(snare, nfeatures = 35000)
#brain <- FindTopFeatures(brain, min.cutoff = "q5")
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 50, reduction.key = "pcs", reduction.name = "pcs")


snare_pca <- snare@reductions$pcs@cell.embeddings %>% as.data.frame()

umap_snare_atac_pca <-
  uwot::umap(snare_pca,
             n_neighbors = n_neighbors,
             min_dist = min_dist)


p_umap_snare_atac_pca <-
  data.frame(umap_snare_atac_pca) %>% 
  ggplot(aes(x = X1, y = X2)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 PCA Normalized ATAC")
#theme(legend.position = "none")


p_umap_snare_atac_pca



#########################################
## C: PCA on Variable Gene ATAC Activity
#########################################



## standard embedding
umap_atac_genes_pca <-
  uwot::umap(data.frame(t(atac_genes_variable)),
             n_neighbors = n_neighbors,
             min_dist = min_dist,
             pca = 50)


p_umap_atac_genes_pca <-
  data.frame(umap_atac_genes_pca) %>% 
  ggplot(aes(x = X1, y = X2)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 PCA Variable Genes ATAC") + 
  theme(legend.position = "none")

p_umap_atac_genes_pca


#########################################
## D: PCA on Variable Gene Expression
#########################################


## standard embedding
umap_rna_genes_pca <-
  uwot::umap(data.frame(t(gex)),
             n_neighbors = n_neighbors,
             min_dist = min_dist,
             pca = 50)


p_umap_rna_genes_pca <-
  data.frame(umap_rna_genes_pca) %>% 
  ggplot(aes(x = X1, y = X2)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 PCA Variable Genes RNA") + 
  theme(legend.position = "none")

p_umap_rna_genes_pca






### 1. SHAREseq --------------------------------------------------------------- 


load(here("data/share_seq/brain/brain.rda"))

### Cell Labels


barcodes.brain <- 
  row.names(brain@meta.data) %>% 
  data.frame()
names(barcodes.brain) <- "barcode"

ct <- fread(here('data/share_seq/brain/celltype_brain.txt'))

ct <- ct %>% 
  mutate_at(vars(atac.bc, rna.bc), 
            ~ gsub("[[:punct:]]", "\\.", .))

celltype_df <- ct %>% 
  dplyr::select(barcode = rna.bc, celltype) %>% 
  dplyr::filter(barcode %in% names(data.frame(as.matrix(brain@assays$RNA@counts))))

barcodes.brain <- 
  barcodes.brain %>% 
  left_join(celltype_df)


############################################################
### A: Standard preprocessing approaches:
############################################################

share_atac <- brain@assays$ATAC@counts %>% as.matrix()

share_atac_norm <-brain@assays$ATAC@data %>% as.matrix()

DefaultAssay(brain) <- "ATAC"

brain <- FindVariableFeatures(brain, nfeatures = 35000)
#brain <- FindTopFeatures(brain, min.cutoff = "q5")
brain <- ScaleData(brain)
brain <- RunPCA(brain, npcs = 30, reduction.key = "pcs", reduction.name = "pcs")


share_pca <- brain@reductions$pcs@cell.embeddings %>% as.data.frame()

umap_share_atac_pca <-
  uwot::umap(share_pca,
             n_neighbors = n_neighbors,
             min_dist = min_dist)


p_umap_share_atac_pca <-
  data.frame(umap_share_atac_pca) %>% 
  mutate(celltype = barcodes.brain$celltype) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", col = "celltype") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 30 PCA Normalized ATAC")
  #theme(legend.position = "none")


p_umap_share_atac_pca


###################################
### A2: ON GENE ACITIVTY MATRIX
###################################

DefaultAssay(brain) <- "ACTIVITY"

brain <- NormalizeData(brain)
brain <- ScaleData(brain)
brain <- RunPCA(brain, npcs = 30, reduction.key = "pcact_", reduction.name = "pcact")


share_pca <- brain@reductions$pcact@cell.embeddings %>% as.data.frame()

umap_share_atac_pca <-
  uwot::umap(share_pca,
             n_neighbors = n_neighbors,
             min_dist = min_dist)


p_umap_share_atac_pca_act <-
  data.frame(umap_share_atac_pca) %>% 
  mutate(celltype = barcodes.brain$celltype) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", col = "celltype") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 30 PCA Normalized ATAC")
#theme(legend.position = "none")


p_umap_share_atac_pca_act




############################################################
### B: LSI Reduction
############################################################

DefaultAssay(brain) <- "ATAC"
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 10)
brain <- RunSVD(brain, reduction.name = "newlsi", reduction.key = "newlsi_")


brain <- RunUMAP(brain, reduction = 'newlsi', 
                 dims = 2:30, reduction.name = 'newumap.atac')

p1 <- DimPlot(brain, reduction = 'newumap.atac') + NoLegend() + ggtitle("ATAC UMAP")

p1

load(here("output/snareseq_intermed.RData"))

snare <- RunUMAP(snare, reduction = 'lsi', 
                 dims = 2:30, reduction.name = 'newumap.atac')

p2 <- DimPlot(snare, reduction = 'newumap.atac') + NoLegend() + ggtitle("ATAC UMAP")

p2

############################################################
### C: SVD on Variable ATAC Peaks
############################################################

## get variable chromatin genes 

DefaultAssay(brain) <- "ATAC"
brain <- FindTopFeatures(brain, min.cutoff = 10)
share_atac_variable <-
  RunTFIDF(brain@assays$ATAC@counts[brain@assays$ATAC@var.features, ])

atac_peaks_variable <- FindVariableFeatures(share_atac_variable, 
                                            selection.method = "disp")

atac_peaks_variable <- 
  row.names(atac_peaks_variable[order(atac_peaks_variable$mvp.dispersion, 
                                      decreasing = TRUE), 
                                , 
                                drop = FALSE][1:3000, ])

share_atac_variable <- share_atac_variable[atac_peaks_variable, ]


### SVD
atac_peaks_svd <- RunSVD(share_atac_variable)

## standard embedding
umap_atac_peaks_svd <-
  uwot::umap(data.frame(atac_peaks_svd@cell.embeddings)[,2:10],
             n_neighbors = n_neighbors,
             min_dist = min_dist)


p_umap_atac_pca <-
  data.frame(umap_atac_peaks_svd) %>% 
  mutate(celltype = barcodes.brain$celltype) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 2-30 SVD approx Data Matrix") #+ 
  #theme(legend.position = "none")

p_umap_atac_pca



############################################################
### D: PCA on Variable ATAC Genes
############################################################

DefaultAssay(brain) <- "ACTIVITY"

brain@assays$ACTIVITY

#brain <- FindVariableFeatures(brain, nfeatures = 3000)
brain <- NormalizeData(brain)
brain <- ScaleData(brain)

## extract variable chromatin genes 
atac_genes_variable <- 
  brain@assays$ACTIVITY@scale.data



## standard embedding
umap_atac_genes_pca <-
  uwot::umap(data.frame(t(atac_genes_variable)),
             n_neighbors = n_neighbors,
             min_dist = min_dist,
             pca = 30)


p_umap_atac_genes_pca <-
  data.frame(umap_atac_genes_pca) %>% 
  mutate(celltype = barcodes.brain$celltype) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 PCA Variable Genes ATAC") + 
  theme(legend.position = "none")

p_umap_atac_genes_pca






### ---------------------------------------------------------------------------
### SIDEseq vs PCA of ATAC ----------------------------------------------------
### ---------------------------------------------------------------------------

##############################################################
### A: SIDEseq on normalized chromatin accessibility values
##############################################################

ref_size = 25
select_method = "cell_embed_sample"


## Get scaled select chromatin peak regions
DefaultAssay(brain) <- "ATAC"

brain <- FindVariableFeatures(brain, nfeatures = 35000)
brain <- NormalizeData(brain)
brain <- ScaleData(brain)

atac_matrix <- brain@assays$ATAC@scale.data



brain_sample <- subset(brain, cells = sample(Cells(brain), 500))

atac_matrix_samp <- brain_sample@assays$ATAC@scale.data

cells_samp <- colnames(atac_matrix_samp)

cell_types_samp <- celltype_df$celltype[sapply(cells_samp, function(x){
  which(x == celltype_df$barcode)})]

ptm <- proc.time()
test_atac_sideseq_ref <- 
  SIDEseqRefSet(expr_matrix = atac_matrix_samp, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = select_method,
                size_ref_set = ref_size,
                n_clust = 25,
                B = 3,
                D = 5,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)
ct_elap <- (proc.time() - ptm)[3]



dissim_atac_sideseq_ref <- stats::as.dist(test_atac_sideseq_ref$dissim_final)


umap_atac_full_sideseq_ref <-
  uwot::umap(dissim_atac_sideseq_ref,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_atac_genes_sideseq <-
  data.frame(umap_atac_full_sideseq_ref) %>% 
  mutate(celltype = cell_types_samp) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("SIDEseq ATAC") + 
  theme(legend.position = "none")

p_umap_atac_genes_sideseq


## Versus PCA Approach:
DefaultAssay(brain_sample) <- "ATAC"
brain_sample <- RunPCA(brain_sample, npcs = 30, 
                       reduction.key = "pcs", reduction.name = "pcs")


share_samp_pca <- brain_sample@reductions$pcs@cell.embeddings %>% as.data.frame()

umap_share_samp_atac_pca <-
  uwot::umap(share_samp_pca,
             n_neighbors = n_neighbors,
             min_dist = min_dist)


p_umap_share_samp_atac_pca <-
  data.frame(umap_share_samp_atac_pca) %>% 
  mutate(celltype = cell_types_samp) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", col = "celltype") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 30 PCA Normalized ATAC")
#theme(legend.position = "none")


p_umap_share_samp_atac_pca







##############################################################
### B: SIDEseq on SNAREseq vs PCA->UMAP
##############################################################

ref_size = 25
select_method = "cell_embed_sample"


## Get scaled select chromatin peak regions
load(here("output/snareseq_final.RData"))
DefaultAssay(snare) <- "ATAC"




### Preprocess
snare_sample <- FindVariableFeatures(snare_sample, nfeatures = 5000)




## Tracking what I've saved so far:

save(snare_samp_pca, file = here("output/embeds/snare_full_pca_1500"))
save(snare_samp_pca, file = here("output/embeds/snare_150Kpeaks_pca_1500"))
save(snare_samp_pca, file = here("output/embeds/snare_100Kpeaks_pca_1500"))
save(snare_samp_pca, file = here("output/embeds/snare_50Kpeaks_pca_1500"))
save(snare_samp_pca, file = here("output/embeds/snare_5Kpeaks_pca_1500"))
save(snare_samp_pca, file = here("output/embeds/snare_full_pca_500"))

p_umap_snare_samp_atac_pca <-
  data.frame(umap_snare_samp_atac_pca) %>% 
  mutate(celltype = snare_sample@meta.data$celltype) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 30 PCA Normalized ATAC")


p_umap_snare_samp_atac_pca


###############################################################################
### PCA for differing amounts of peaks
###############################################################################

UMAPPlot <- function(data, title, celltypes) {
  umap_snare_samp_atac_pca <-
    uwot::umap(data,
               n_neighbors = n_neighbors,
               min_dist = min_dist)
  
  p <-
    data.frame(umap_snare_samp_atac_pca) %>% 
    mutate(celltype = celltypes) %>%
    ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
    labs(x="UMAP1", y="UMAP2", col = "Cell Type") + 
    theme_bw() + 
    ggtitle(title)
  return(p)
}


## TODO, put in loop
plot_list_num_peaks <- list()

load(here("output/embeds/snare_full_pca_1500"))
plot_list_num_peaks[[1]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "All 244,544 Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_150Kpeaks_pca_1500"))
plot_list_num_peaks[[2]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "150,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_100Kpeaks_pca_1500"))
plot_list_num_peaks[[3]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "100,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_50Kpeaks_pca_1500"))
plot_list_num_peaks[[4]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "50,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_5Kpeaks_pca_1500"))
plot_list_num_peaks[[5]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "5,000 Most Variable Peaks",
           celltypes = snare_sample@meta.data$celltype)

ggarrange(plotlist = plot_list_num_peaks, nrow = 5,
          common.legend = TRUE, legend = "bottom")


ggsave(here("output/figures/snare_atac_num_peaks.png"), 
       width = 5, height = 25)

save(plot_list_num_peaks, 
     file = here("output/figures/plot_objects/atac_num_peaks_pca.RData"))


### For RNA-seq ---------------------------------------------------------------
### TODO: clean pipeline

DefaultAssay(snare) <- "RNA"

gene_counts <- c(dim(snare@assays$RNA@counts)[1], 
                 25000, 10000, 2500, 500)


## TODO: test if we really need to resample each time.
for(g in gene_counts) {
  print(g)
  set.seed(209827)
  snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))
  DefaultAssay(snare) <- "RNA"
  ### Preprocess
  snare_sample <- FindVariableFeatures(snare_sample, nfeatures = g)
  snare_sample <- NormalizeData(snare_sample)
  snare_sample <- ScaleData(snare_sample)
  
  snare_sample <- RunPCA(snare_sample, npcs = 50, 
                         reduction.key = "pcs_", reduction.name = "pcs")
  
  
  snare_samp_pca <- snare_sample@reductions$pcs@cell.embeddings %>% as.data.frame()
  save(snare_samp_pca, file = here(paste0("output/embeds/snare_", g, "_genes_pca_1500.RData")))
}


plot_list_n_genes <- list()

## TODO loop
load(here("output/embeds/snare_33160_genes_pca_1500.RData"))
plot_list_n_genes[[1]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "All 33,160 Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_25000_genes_pca_1500.RData"))
plot_list_n_genes[[2]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "25,000 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_10000_genes_pca_1500.RData"))
plot_list_n_genes[[3]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "10,000 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_2500_genes_pca_1500.RData"))
plot_list_n_genes[[4]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "2,500 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

load(here("output/embeds/snare_500_genes_pca_1500.RData"))
plot_list_n_genes[[5]] <- 
  UMAPPlot(data  = snare_samp_pca, 
           title = "500 Most Variable Genes",
           celltypes = snare_sample@meta.data$celltype)

ggarrange(plotlist = list(plot_list_n_genes[[1]], plot_list_num_peaks[[1]], 
                          plot_list_n_genes[[2]], plot_list_num_peaks[[2]], 
                          plot_list_n_genes[[3]], plot_list_num_peaks[[3]], 
                          plot_list_n_genes[[4]], plot_list_num_peaks[[4]], 
                          plot_list_n_genes[[5]], plot_list_num_peaks[[5]]),
          nrow = 5, ncol = 2,
          common.legend = TRUE, legend = "bottom")





###############################################################################
### Testing hyperparameters of UMAP
###############################################################################

neighbors_vec <- c(5, 15, 25, 50)
min_dist_vec <- c(1e-3, 0.01, 0.1, 0.5)
plot_list <- list()

avg_sil_vec <- c()

i <- 1
for(n in neighbors_vec) {
  n_neighbors = n
  for(m in min_dist_vec) {
    print(i)
    min_dist = m
    umap_snare_samp_atac_pca <-
      uwot::umap(snare_samp_pca,
                 n_neighbors = n_neighbors,
                 min_dist = min_dist)
    
    
    ## compute Silhouette Scores
    sil_scores <- cluster::silhouette(as.integer(factor(snare_sample@meta.data$celltype)),
                                      dist = stats::dist(umap_snare_samp_atac_pca, 
                                                         method = "euclidean"))[,3]
    avg_sil_vec[i] <- mean(sil_scores)
    
    p <-
      data.frame(umap_snare_samp_atac_pca) %>% 
      mutate(celltype = snare_sample@meta.data$celltype) %>%
      ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
      labs(x="UMAP1", y="UMAP2", col = "Cell Type") + 
      theme_bw() + 
      ggtitle(paste(n, "Neighbors,", m, "Minimum Dist,", 
                    round(avg_sil_vec[i], 2), "Sil. Score"))
    
    
    
    plot_list[[i]] <- p
    
    i <- i+1
    
  }
}


ggarrange(plotlist = plot_list, ncol = 4, nrow = 4,
          common.legend = TRUE, legend = "right")


ggsave(here("output/figures/snare_atac_pca_hypers.png"), 
       width = 20, height = 14)


save(plot_list, 
     file = here("output/figures/plot_objects/atac_pca_hyper_grid.RData"))




###############################################################################
###############################################################################
###############################################################################
###
### 3: 10X Human PBMC Cells 
###
###############################################################################
###############################################################################
###############################################################################




##############################################################
### B: SIDEseq on SNAREseq vs PCA->UMAP
##############################################################

ref_size = 25
select_method = "cell_embed_sample"


## Get scaled select chromatin peak regions
load(here("output/snareseq_final.RData"))
DefaultAssay(snare) <- "ATAC"



set.seed(209827)
snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))

### Preprocess
snare_sample <- FindVariableFeatures(snare_sample, nfeatures = 50000)
snare_sample <- NormalizeData(snare_sample)
snare_sample <- ScaleData(snare_sample)


atac_matrix_samp <- snare_sample@assays$ATAC@scale.data

cells_samp <- colnames(atac_matrix_samp)

## TODO: get SNAREseq samples
# cell_types_samp <- celltype_df$celltype[sapply(cells_samp, function(x){
#   which(x == celltype_df$barcode)})]

ptm <- proc.time()
test_atac_sideseq_ref <- 
  SIDEseqRefSet(expr_matrix = atac_matrix_samp, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = select_method,
                size_ref_set = ref_size,
                n_clust = 25,
                B = 3,
                D = 5,
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)
ct_elap <- (proc.time() - ptm)[3]



dissim_atac_sideseq_ref <- stats::as.dist(test_atac_sideseq_ref$dissim_final)


umap_atac_full_sideseq_ref <-
  uwot::umap(dissim_atac_sideseq_ref,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_atac_genes_sideseq <-
  data.frame(umap_atac_full_sideseq_ref) %>% 
  mutate(celltype = cell_types_samp) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("SIDEseq ATAC") + 
  theme(legend.position = "none")

p_umap_atac_genes_sideseq


#####################################
## Versus PCA Approach:
DefaultAssay(snare_sample) <- "ATAC"
snare_sample <- RunPCA(snare_sample, npcs = 50, 
                       reduction.key = "pcs_", reduction.name = "pcs")


snare_samp_pca <- snare_sample@reductions$pcs@cell.embeddings %>% as.data.frame()

umap_snare_samp_atac_pca <-
  uwot::umap(snare_samp_pca,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

## Tracking what I've saved so far:

#save(snare_samp_pca, file = here("output/embeds/snare_full_pca_1500"))
#save(snare_samp_pca, file = here("output/embeds/snare_150Kpeaks_pca_1500"))
#save(snare_samp_pca, file = here("output/embeds/snare_100Kpeaks_pca_1500"))
#save(snare_samp_pca, file = here("output/embeds/snare_50Kpeaks_pca_1500"))
#save(snare_samp_pca, file = here("output/embeds/snare_full_pca_500"))

p_umap_snare_samp_atac_pca <-
  data.frame(umap_snare_samp_atac_pca) %>% 
  mutate(celltype = snare_sample@meta.data$celltype) %>%
  ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 30 PCA Normalized ATAC")


p_umap_snare_samp_atac_pca


###############################################################################
### Testing hyperparameters of UMAP
###############################################################################

neighbors_vec <- c(5, 15, 25, 50)
min_dist_vec <- c(1e-3, 0.01, 0.1, 0.5)
plot_list <- list()
i <- 1
for(n in neighbors_vec) {
  n_neighbors = n
  for(m in min_dist_vec) {
    print(i)
    min_dist = m
    umap_snare_samp_atac_pca <-
      uwot::umap(snare_samp_pca,
                 n_neighbors = n_neighbors,
                 min_dist = min_dist)
    
    p <-
      data.frame(umap_snare_samp_atac_pca) %>% 
      mutate(celltype = snare_sample@meta.data$celltype) %>%
      ggplot(aes(x = X1, y = X2, col = celltype)) + geom_point() + 
      labs(x="UMAP1", y="UMAP2") + 
      theme_bw() + 
      ggtitle(paste(n, "neighbors", m, "min_dist"))
    
    plot_list[[i]] <- p
    i <- i+1
    
  }
}


ggarrange(plotlist = plot_list, ncol = 4, nrow = 4)















