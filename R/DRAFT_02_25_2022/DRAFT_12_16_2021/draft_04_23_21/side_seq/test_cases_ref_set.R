## Testing the cell reference set approach to SIDEseq 
source(here("R/side_seq/side_seq_ref_set.R"))
source(here("R/side_seq/evaluate.R"))

`%p%` <- function(x,y) {paste0(x,y)}

## Add cell labels
load(here("data/share_seq/brain/brain.rda"))


## Access normalized count data
umi_counts_norm <- 
  as.matrix(brain@assays$SCT@data) %>% 
  as.data.frame()

## add annotated cell types
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
  dplyr::filter(barcode %in% names(umi_counts_norm))

barcodes.brain <- 
  barcodes.brain %>% 
  left_join(celltype_df)



brain@meta.data <-cbind(brain@meta.data, celltype = barcodes.brain$celltype)


## 
cell_labels <- barcodes.brain


## Trim full data to higher mean expression genes.
GENE_PCTILE_CUTOFF <- 0.90

umi_counts_norm$mean_express <- 
  apply(umi_counts_norm, 1, mean)

umi_counts_high_expr <- 
  umi_counts_norm %>%
  dplyr::filter(mean_express >= quantile(umi_counts_norm$mean_express, GENE_PCTILE_CUTOFF))



test_expr_matrix <- umi_counts_high_expr[, sample(ncol(umi_counts_high_expr), 300)]


test_ref_set_size = 40



set.seed(2020)
test_run_upd_mat <- 
  SIDEseqRefSet(expr_matrix = test_expr_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 300,
                ## cell reference set selection parameters
                selection_method = "cell_embed_sample",
                size_ref_set = test_ref_set_size,
                n_clust = 10,
                B = 10,
                prior_mat_method = "new", ## c("new", "average", "both")
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)

set.seed(2020)
test_run_avg_mat <- 
  SIDEseqRefSet(expr_matrix = test_expr_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 300,
                ## cell reference set selection parameters
                selection_method = "cell_embed_sample",
                size_ref_set = 10,
                n_clust = 10,
                B = 1,
                prior_mat_method = "average", ## c("new", "average", "both")
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)



## Let's look at our results!
test_run_upd_mat$convergence_res

test_run_avg_mat$convergence_res



## Plot results:

spect_clust <- spectralCluster(test_run_upd_mat$sim_matrix,
                               n_clust = n_clust)

spect_clust_data <- data.frame(spect_clust$G) %>% 
  mutate(cluster = spect_clust$clusters)

## get cell types
test_cells <- names(test_expr_matrix)

test_celltypes <- cell_labels$celltype[sapply(test_cells, function(x){
  which(x == cell_labels$barcode)})]

spect_clust_data <- spect_clust_data %>% 
  mutate(celltype = test_celltypes)

p1 <-
  spect_clust_data %>% 
  ggplot() + geom_point(aes(x=X1, y = X2, color = factor(cluster)))

p2 <-
  spect_clust_data %>% 
  ggplot() + geom_point(aes(x=X1, y = X2, color = factor(celltype)))


ggarrange(p1, p2, ncol = 2)#, widths = c(.5,.02,.5))


## TODO: warm-starting mechanism


### Differential expression matrix size plot ----------------------------------




## TODO: examine improvement when filtering to highly variable genes before any SIDEseq comp.


## Run evaluation loop:
## get cell types


### Subset cell types

set.seed(2020)
test_n_cells <- 300

test_expr_matrix <- 
  umi_counts_high_expr[, sample(ncol(umi_counts_high_expr), 
                                test_n_cells)]

test_cells <- names(test_expr_matrix)

test_celltypes <- cell_labels$celltype[sapply(test_cells, function(x){
  which(x == cell_labels$barcode)})]


## Expand grid of options:
select_methods <- c("random", 
                    "total_reads",
                    "n_genes_detected",
                    "n_genes_low_expr",
                    "n_highly_variable_genes",
                    "cell_embed_sample")

ref_set_sizes <- c(15, seq(25, 300, by = 25))

B_opts <- c(1, 5)

prior_mat_methods <- c("new", "average")


eval_opt_df <- 
  expand.grid(selection_method = select_methods, 
              size_ref_set     = ref_set_sizes, 
              B                = B_opts, 
              prior_mat_method = prior_mat_methods,
              stringsAsFactors = FALSE)


## filtering options we don't need to run
eval_opt_df <- 
  eval_opt_df %>% 
  dplyr::filter(!(selection_method %in% c("random") & 
                    prior_mat_method == "new")) %>% 
  dplyr::filter(!(selection_method %in% c("total_reads",
                                          "n_genes_detected",
                                          "n_genes_low_expr",
                                          "n_highly_variable_genes") & 
                    (prior_mat_method == "average" | B > 1))) %>% 
  dplyr::filter(!(prior_mat_method == "average" & B == 1)) %>% 
  dplyr::filter(!B>1) %>%
  arrange(size_ref_set)



### Initialize
## TODO cleaner structure.
eval_df <- 
  data.frame(diff_expr_method = NA, 
             similarity_method = NA,
             n_top_genes = NA, 
             selection_method = NA,
             size_ref_set = NA,
             n_clust = NA,
             B = NA,
             prior_mat_method = NA,
             parallelize = NA,
             n_cores = NA,
             run_time = NA,
             avg_sil_width = NA,
             ## TODO: change to adj_rand_index
             rand_index = rep(NA, nrow(eval_opt_df)),
             stringsAsFactors = FALSE)



for(r in seq_len(nrow(eval_opt_df))) {
  cat(r, "of", nrow(eval_opt_df), "...")
  eval_df[r,] <-
    evaluateSIDEseqRun(cell_types = test_celltypes,
                       n_runs = 1,
                       expr_matrix = test_expr_matrix,
                       diff_expr_method = "diff_expr_norm",
                       similarity_method = "n_intersect",
                       n_top_genes = 300,
                       ## Selection method changes::
                       selection_method = eval_opt_df$selection_method[r],
                       ## Ref set size changes::
                       size_ref_set = eval_opt_df$size_ref_set[r],
                       n_clust = 10,
                       ## B changes::
                       B = eval_opt_df$B[r],
                       prior_mat_method = eval_opt_df$prior_mat_method[r],
                       parallelize = TRUE,
                       n_cores = detectCores()-1,
                       verbose = TRUE)
  
  saveRDS(eval_df, here("output/SIDE_seq_ref_set_test_runs_300_genes.RDS"))
  
}

full_side_seq_res <-
  evaluateSIDEseqRun(cell_types = test_celltypes,
                     n_runs = 1,
                     expr_matrix = test_expr_matrix,
                     diff_expr_method = "diff_expr_norm",
                     similarity_method = "n_intersect",
                     n_top_genes = 300,
                     ## Selection method changes::
                     selection_method = "random",
                     ## Ref set size changes::
                     size_ref_set = 300,
                     n_clust = 10,
                     ## B changes::
                     B = 1,
                     prior_mat_method = "new",
                     parallelize = TRUE,
                     n_cores = detectCores()-1,
                     verbose = TRUE)


## Plot computation time
# test_eval_df <-
#   eval_df %>% 
#   mutate(selection_method = case_when(
#     selection_method == 1 ~ "random",
#     selection_method == 2 ~ "total_reads",
#     selection_method == 3 ~ "n_genes_detected",
#     selection_method == 4 ~ "n_genes_low_expr",
#     selection_method == 5 ~ "n_highly_variable_genes",
#     selection_method == 6 ~ "cell_embed_sample"
#   ))



plot_eval_df <- 
  eval_df %>% 
  bind_rows(full_side_seq_res %>% 
              mutate(selection_method = "All Cells")) %>%
  dplyr::filter(!is.na(B)) %>%
  mutate(B = paste("B = ", B)) %>%
  mutate(prior_mat_method = ifelse(prior_mat_method == "average",
                                   "Mean", "Last Iter's Sim Matrix"))

### Computation time:
plot_eval_df %>%
  ggplot(aes(x = size_ref_set, y = run_time, 
             color = selection_method)) + 
  geom_point() + geom_line() + 
  facet_wrap(~prior_mat_method + B, scales = "free_y") + 
  theme_bw() + 
  labs(x = "Size of Cell Ref Set", y = "Run Time", 
       color = "Ref Set Select Method")
  

## Clustering Metrics
plot_eval_df %>%
  ggplot(aes(x = size_ref_set, y = avg_sil_width, 
             color = selection_method)) + 
  geom_point() + geom_line() + 
  facet_wrap(~prior_mat_method + B, scales = "free_y") + 
  theme_bw() + 
  labs(x = "Size of Cell Ref Set", y = "Average Sil Width", 
       color = "Ref Set Select Method")


plot_eval_df %>%
  ggplot(aes(x = size_ref_set, y = rand_index, 
             color = selection_method)) + 
  geom_point() + geom_line() + 
  facet_wrap(~prior_mat_method + B, scales = "free_y") + 
  theme_bw() + 
  labs(x = "Size of Cell Ref Set", y = "Adj Rand Index", 
       color = "Ref Set Select Method")


#saveRDS(eval_df, here("output/SIDE_seq_ref_set_test_runs.RDS"))



### Examining UMAP Embeddings -------------------------------------------------


n_neighbors <- 15L
min_dist <- 0.01


## Pearson correlation
sim_matrix <- rcorr(as.matrix(test_expr_matrix),
                    type = "pearson")$r

diag(sim_matrix) <- 0

dissim_matrix <- max(sim_matrix) - sim_matrix
diag(dissim_matrix) <- 0

dissim_matrix <- stats::as.dist(dissim_matrix)

set.seed(2021)
umap_pearson <-
  uwot::umap(dissim_matrix,
             n_neighbors = 20,
             min_dist = min_dist)

p_umap_pearson <-
  data.frame(umap_pearson) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Pearson Corr.") + 
  theme(legend.position = "none")


### spearman corr
sim_matrix <- rcorr(as.matrix(test_expr_matrix),
                    type = "spearman")$r

diag(sim_matrix) <- 0

dissim_matrix <- max(sim_matrix) - sim_matrix
diag(dissim_matrix) <- 0

dissim_matrix <- stats::as.dist(dissim_matrix)

set.seed(2021)
umap_spearman <-
  uwot::umap(dissim_matrix,
             n_neighbors = 20,
             min_dist = .01)

p_umap_spearman <-
  data.frame(umap_spearman) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Spearman Corr.")


## vanilla
set.seed(2021)
umap_vanilla <-
  uwot::umap(t(test_expr_matrix),
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_vanilla <-
  data.frame(umap_vanilla) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Euclid Dist on Original Data Matrix") + 
  theme(legend.position = "none")


umap_pca <-
  uwot::umap(t(test_expr_matrix), pca = 50,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_pca <-
  data.frame(umap_pca) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 PCs Data Matrix") + 
  theme(legend.position = "none")


### testing manual distance matrix, euclidean

umap_euclid_dist <- stats::dist(t(test_expr_matrix),
                                method = "euclidean")

set.seed(2021)
umap_man_dist <-
  uwot::umap(umap_euclid_dist,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_man_dist <-
  data.frame(umap_man_dist) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw()



set.seed(2021)
umap_auto_pearson <-
  uwot::umap(t(test_expr_matrix),
             n_neighbors = n_neighbors,
             min_dist = min_dist,
             nn_method = "annoy",
             metric = "correlation")

p_umap_auto_pearson <-
  data.frame(umap_auto_pearson) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Pearson Corr.") + 
  theme(legend.position = "none")


ggarrange(p_umap_vanilla, p_umap_spearman, p_umap_pca,
          p_umap_man_dist, p_umap_pearson,
          nrow = 2,
          ncol = 3)




# ggarrange(p_umap_pearson, p_umap_auto_pearson,
#           ncol = 2)



### SIDEseq result ------------------------------------------------------------

test_sideseq_res <-
  SIDEseqSimult(expr_matrix = test_expr_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                parallelize = TRUE)

dissim_sideseq <- stats::as.dist(test_sideseq_res)


umap_full_sideseq <-
  uwot::umap(dissim_sideseq,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_full_sideseq <-
  data.frame(umap_full_sideseq) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("SIDEseq Full") + 
  theme(legend.position = "none")


## TODO: create n_top_genes constant for these plots


### SIDEseq refset (1 run) ----------------------------------------------------
ref_size = 75
select_method = "random"

test_sideseq_ref <- 
  SIDEseqRefSet(expr_matrix = test_expr_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = select_method,
                size_ref_set = ref_size,
                n_clust = 10,
                B = 1,
                prior_mat_method = "new", ## c("new", "average", "both")
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)


dissim_sideseq_ref <- as.dist(test_sideseq_ref$dissim_matrix)

umap_full_sideseq_ref <-
  uwot::umap(dissim_sideseq_ref,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_full_sideseq_ref <-
  data.frame(umap_full_sideseq_ref) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() +
  ggtitle("SIDEseq Ref Set " %p% 
            ref_size %p% " Cells " %p% select_method %p% " Selection") +
  theme(legend.position = "none")



### SIDEseq refset (10 runs) --------------------------------------------------

select_method = "cell_embed_sample"

test_sideseq_ref_iter <- 
  SIDEseqRefSet(expr_matrix = test_expr_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = select_method,
                size_ref_set = ref_size,
                n_clust = 10,
                B = 10,
                clust_repeats = 3,
                purity_tol = 0.95,
                prior_mat_method = "new", ## c("new", "average", "both")
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)


dissim_sideseq_ref_iter <- as.dist(test_sideseq_ref_iter$dissim_matrix)

umap_full_sideseq_ref_iter <-
  uwot::umap(dissim_sideseq_ref_iter,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_full_sideseq_ref_iter <-
  data.frame(umap_full_sideseq_ref_iter) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() +
  ggtitle("SIDEseq Ref Set " %p% 
            ref_size %p% " Cells " %p% "Iterative Stratified Sampling of Clusters") +
  theme(legend.position = "none")





## Plot results RNAseq --------------------------------------------------------
ggarrange(p_umap_vanilla, p_umap_pca, p_umap_spearman,
          p_umap_full_sideseq, p_umap_full_sideseq_ref, p_umap_pearson,
          nrow = 2,
          ncol = 3,
          common.legend = TRUE, legend = "right")

ggsave(here("output/figures/sideseqRNA300cells.png"),
       width = 12, height = 6)











###############################################################################
## ATAC-seq clustering attempts -----------------------------------------------
###############################################################################
###############################################################################

## Get normalized atac-seq data
#atac_seq_peak_counts <- as.matrix(brain@assays$ATAC@counts)
peak.count= readMM(here("data/share_seq/brain/GSM4156599_brain.counts.txt.gz"))
peak.bed=read.delim(here('data/share_seq/brain/GSM4156599_brain.peaks.bed.gz'), 
                    header=FALSE)
peak.granges=GRanges(seqnames=peak.bed$V1, 
                     ranges=IRanges(st=peak.bed$V2, end=peak.bed$V3))
rm(peak.bed)
rownames(peak.count)=paste(peak.granges)
colnames(peak.count)=ct$rna.bc # we will use the rna barcode from here and on

grange.use <- seqnames(peak.granges) %in% standardChromosomes(peak.granges)
peak.count <- peak.count[as.vector(grange.use), ]
peak.granges <- peak.granges[as.vector(grange.use)]
rm(grange.use)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Yuriko mentioned that there is issue loading the fragment file
# Leave the fragment file out for now.
# frag.file <- "GSM4156599_brain.atac.fragments.bed.gz"
brain_test <- CreateChromatinAssay(
  counts = peak.count,
  sep = c(":", "-"),
  genome = 'mm10',
  # fragments = frag.file,
  min.cells = 10,
  annotation = annotations)


# brain_test <- CreateSeuratObject(counts = brain@assays$RNA@counts)


## TODO: Retrying Try the given data preprocessing
DefaultAssay(brain) <- "ATAC"

brain <- FindTopFeatures(brain, min.cutoff = 10)
brain <- RunTFIDF(brain)



atac_gene_peak_counts_scaled <- as.matrix(brain@assays$ACTIVITY@scale.data)

atac_gene_peak_counts_unscaled <- as.matrix(brain@assays$ACTIVITY@data)

## Identify high variance ATAC genes 
atac_sd <- apply(atac_gene_peak_counts_scaled, 1, sd)

GENE_PCTILE_CUTOFF <- 0.90

atac_peak_counts_high_sd <- 
  atac_gene_peak_counts_scaled[atac_sd > quantile(atac_sd, GENE_PCTILE_CUTOFF), ]



test_atac_matrix <- 
  atac_peak_counts_high_sd[, which(names(umi_counts_high_expr) %in% 
                                     test_cells)]


### Now produce clustering results for the ATACseq data:


### Examining UMAP Embeddings -------------------------------------------------

## TODO: need to tidy this up to be functional.


## Pearson correlation
sim_matrix <- rcorr(as.matrix(test_atac_matrix),
                    type = "pearson")$r

diag(sim_matrix) <- 0

dissim_matrix <- max(sim_matrix) - sim_matrix
diag(dissim_matrix) <- 0

dissim_matrix <- stats::as.dist(dissim_matrix)

set.seed(2021)

umap_atac_pearson <-
  uwot::umap(dissim_matrix,
             n_neighbors = 20,
             min_dist = min_dist)

p_umap_atac_pearson <-
  data.frame(umap_atac_pearson) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Pearson Corr.") + 
  theme(legend.position = "none")


### spearman corr
sim_matrix <- rcorr(as.matrix(test_atac_matrix),
                    type = "spearman")$r

diag(sim_matrix) <- 0

dissim_matrix <- max(sim_matrix) - sim_matrix
diag(dissim_matrix) <- 0

dissim_matrix <- stats::as.dist(dissim_matrix)

set.seed(2021)
umap_atac_spearman <-
  uwot::umap(dissim_matrix,
             n_neighbors = 20,
             min_dist = .01)

p_umap_atac_spearman <-
  data.frame(umap_atac_spearman) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Spearman Corr.")


## vanilla
set.seed(2021)
umap_atac_vanilla <-
  uwot::umap(t(test_atac_matrix),
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_atac_vanilla <-
  data.frame(umap_atac_vanilla) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Euclid Dist on Original Data Matrix") + 
  theme(legend.position = "none")


umap_atac_pca <-
  uwot::umap(t(test_atac_matrix), pca = 50,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_atac_pca <-
  data.frame(umap_atac_pca) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("Euclid Dist w 50 PCs Data Matrix") + 
  theme(legend.position = "none")





### SIDEseq result ------------------------------------------------------------

test_atac_sideseq_res <-
  SIDEseqSimult(expr_matrix = test_atac_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                parallelize = TRUE)

dissim_atac_sideseq <- stats::as.dist(test_atac_sideseq_res)


umap_atac_full_sideseq <-
  uwot::umap(dissim_atac_sideseq,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_umap_atac_full_sideseq <-
  data.frame(umap_atac_full_sideseq) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() + 
  ggtitle("SIDEseq Full") + 
  theme(legend.position = "none")


## TODO: create n_top_genes constant for these plots


### SIDEseq refset (1 run) ----------------------------------------------------
ref_size = 75
select_method = "random"

test_atac_sideseq_ref <- 
  SIDEseqRefSet(expr_matrix = test_atac_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = select_method,
                size_ref_set = ref_size,
                n_clust = 10,
                B = 1,
                prior_mat_method = "new", ## c("new", "average", "both")
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)


dissim_atac_sideseq_ref <- max(test_atac_sideseq_ref$sim_matrix) - test_atac_sideseq_ref$sim_matrix
diag(dissim_atac_sideseq_ref) <- 0

dissim_atac_sideseq_ref <- stats::as.dist(dissim_atac_sideseq_ref)


umap_atac_full_sideseq_ref <-
  uwot::umap(dissim_atac_sideseq_ref,
             n_neighbors = n_neighbors,
             min_dist = min_dist)

p_atac_umap_full_sideseq_ref <-
  data.frame(umap_atac_full_sideseq_ref) %>% 
  mutate(cluster = factor(test_celltypes)) %>%
  ggplot(aes(x = X1, y = X2, color = cluster)) + geom_point() + 
  labs(x="UMAP1", y="UMAP2", color="Cell Type") + 
  theme_bw() +
  ggtitle("SIDEseq Ref Set " %p% 
            ref_size %p% " Cells " %p% select_method %p% " Selection") +
  theme(legend.position = "none")



### SIDEseq refset (10 iters) -------------------------------------------------






## Plot results RNAseq --------------------------------------------------------
ggarrange(p_umap_vanilla, p_umap_pca, p_umap_spearman,
          p_umap_full_sideseq, p_umap_full_sideseq_ref, p_umap_pearson,
          nrow = 2,
          ncol = 3,
          common.legend = TRUE, legend = "right")

ggsave(here("output/figures/sideseqRNA300cells.png"),
       width = 12, height = 6)






### Weighted nearest neighbor approach ----------------------------------------  



## TODO: unnormalized distance metric

## TODO: examine cluster convergence


## TODO: Umap fair hyperparameter search, Ask Phil!


## TODO: build simulation from Side-seq package then extend to realistic multinomial proj.
  
## TODO: how do results change when subsetting to highly variable genes?



### Compute time for similarity scores ----------------------------------------
ref_size = 50
ptm <- proc.time()

test_sim_time_stand <- 
  SIDEseqRefSet(expr_matrix = test_expr_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "n_intersect",
                n_top_genes = 150,
                ## cell reference set selection parameters
                selection_method = "random",
                size_ref_set = ref_size,
                n_clust = 10,
                B = 1,
                prior_mat_method = "new", ## c("new", "average", "both")
                ## parallelization parameters
                parallelize = TRUE,
                n_cores = detectCores()-1,
                ## other
                verbose = TRUE)


ct_elap1 <- (proc.time() - ptm)[3]


ptm <- proc.time()

## Too slow:
# test_sim_time_rbo <- 
#   SIDEseqRefSet(expr_matrix = test_expr_matrix, 
#                 diff_expr_method = "diff_expr_norm",
#                 similarity_method = "rbo",
#                 n_top_genes = 150,
#                 ## cell reference set selection parameters
#                 selection_method = "random",
#                 size_ref_set = ref_size,
#                 n_clust = 10,
#                 B = 1,
#                 prior_mat_method = "new", ## c("new", "average", "both")
#                 ## parallelization parameters
#                 parallelize = TRUE,
#                 n_cores = detectCores()-1,
#                 ## other
#                 verbose = TRUE)
# 
# 
# ct_elap2 <- (proc.time() - ptm)[3]

