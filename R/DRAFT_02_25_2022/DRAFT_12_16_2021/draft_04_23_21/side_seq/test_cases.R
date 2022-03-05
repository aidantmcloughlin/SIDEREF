### Test cases -----------------------------------------------------------------

## SideSEQ implementation packages
library(utils)
library(parallel)
library(bigmemory)
## Data cleaning, plotting
library(tidyverse)
library(janitor)
library(data.table)
library(here)
library(readxl)
library(magrittr)
library(forcats)
library(GGally)
library(ggpubr)
## BioConductor packages
library(Seurat)
library(Signac)
## Statistics / clustering
library(Hmisc) 
library(diffusr)
library(statmod)
library(LICORS)
library(Spectrum)
library(cluster)
library(mclust)
library(uwot)
library(gespeR)

source(here('R/side_seq/side_seq.R'))
source(here('R/side_seq/side_seq_iter.R'))
source(here('R/side_seq/side_seq_block.R'))
source(here('R/side_seq/random_sampling_methods.R'))
source(here('R/side_seq/computeDiffExprMat.R'))
source(here('R/side_seq/setSimilarityMeasure.R'))
load(here("data/share_seq/brain/brain.rda"))


## Access normalized count data
umi_counts_norm <- 
  as.matrix(brain@assays$SCT@data) %>% 
  as.data.frame()

# Just a check that the counts are just log transformed.
# umi_counts_counts <- 
#   as.matrix(brain@assays$RNA@data) %>% 
#   as.data.frame()
# 
# umi_counts_filtered <- 
#   umi_counts_counts[row.names(umi_counts_norm),]
# 
# row.names(umi_counts_filtered) <- NULL
# 
# check <- umi_counts_filtered[,1:100] %>% as.matrix() %>% c()
# 
# only_1 <- check[check==1]
# 
# check2 <- umi_counts_norm[,1:100] %>% as.matrix() %>% c()


# 
# example_matrix <- 
#   umi_counts_subset


## Trim full data to higher mean expression genes.
GENE_PCTILE_CUTOFF <- 0.90

umi_counts_norm$mean_express <- 
  apply(umi_counts_norm, 1, mean)

umi_counts_high_expr <- 
  umi_counts_norm %>%
  dplyr::filter(mean_express >= quantile(umi_counts_norm$mean_express, GENE_PCTILE_CUTOFF))

expr_matrix <- umi_counts_high_expr[,1:25]


small_matrix <- 
  expr_matrix[,1:6]


set.seed(200)
small_no_par <-
  SIDEseqSimult(expr_matrix = small_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "top_n_gene_intersect",
                n_top_genes = 20,
                parallelize = FALSE)

set.seed(200)
small_par <-
  SIDEseqSimult(expr_matrix = small_matrix, 
                diff_expr_method = "diff_expr_norm",
                similarity_method = "top_n_gene_intersect",
                n_top_genes = 20,
                parallelize = TRUE)

set.seed(200)
small_iter_no_par <-
  SIDEseqIterated(expr_matrix = small_matrix, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "top_n_gene_intersect",
                  n_top_genes = 20,
                  parallelize = FALSE)


set.seed(200)
small_iter_no_par <-
  SIDEseqIterated(expr_matrix = small_matrix, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "top_n_gene_intersect",
                  n_top_genes = 20,
                  parallelize = TRUE)


set.seed(200)
small_iter_par_rand <-
  SIDEseqIterated(expr_matrix = small_matrix, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "top_n_gene_intersect",
                  n_top_genes = 20,
                  parallelize = TRUE)



set.seed(200)
small_block_par <-
  SIDEseqBlockSimult(expr_matrix = small_matrix, 
                     diff_expr_method = "diff_expr_norm",
                     similarity_method = "top_n_gene_intersect",
                     n_top_genes = 20,
                     block_memory = 3e3,
                     parallelize = TRUE)

## reproducible:
set.seed(200)
small_block_par_rep <-
  SIDEseqBlockSimult(expr_matrix = small_matrix, 
                     diff_expr_method = "diff_expr_norm",
                     similarity_method = "top_n_gene_intersect",
                     n_top_genes = 20,
                     block_memory = 3e3,==
                     parallelize = TRUE)

small_block_par
small_block_par_rep

## test usage of multiple blocks:
set.seed(200)
small_block_par_multiblock <-
  SIDEseqBlockSimult(expr_matrix = small_matrix, 
                     diff_expr_method = "diff_expr_norm",
                     similarity_method = "top_n_gene_intersect",
                     n_top_genes = 20,
                     block_memory = 0.001,
                     parallelize = TRUE)


## Plot computation time as umi.counts.matrix expands -------------------------


num_cells <- seq(30, 400, by = 10)

## initialize the time benchmark dataset
compute_time_df <-
  data.frame(num_cells,
             #simult_no_para = 0,
             simult_para = 0,
             #iterate_no_para = 0,
             iterate_para = 0)

counter <- 1

set.seed(303030)

for(c in num_cells) {
  print(c)
  ## subset to number of cells of interest
  counts_subset <- umi_counts_high_expr[,sample(ncol(umi_counts_high_expr), c)]
  
  ## Track run time of each function option and store in DF
  
  ## Simultaneous compute, no parallelization
  # ptm <- proc.time()
  # side_seq_res <-
  #   SIDEseqSimult(expr_matrix = counts_subset, 
  #                 diff_expr_method = "diff_expr_norm",
  #                 similarity_method = "top_n_gene_intersect",
  #                 n_top_genes = 1000,
  #                 parallelize = FALSE)
  # ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  # 
  # compute_time_df$simult_no_para[counter] <- ct_elap
  
  ## Simultaneous compute, parallelization
  ptm <- proc.time()
  side_seq_res <-
    SIDEseqSimult(expr_matrix = counts_subset, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "top_n_gene_intersect",
                  n_top_genes = 1000,
                  parallelize = TRUE,
                  n_cores = 6)
  ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  
  compute_time_df$simult_para[counter] <- ct_elap
  
  ## Iterated, no parallelization
  # ptm <- proc.time()
  # side_seq_res <-
  #   SIDEseqIterated(expr_matrix = counts_subset, 
  #                   diff_expr_method = "diff_expr_norm",
  #                   similarity_method = "top_n_gene_intersect",
  #                   n_top_genes = 1000,
  #                   parallelize = FALSE)
  # ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  # 
  # 
  # compute_time_df$iterate_no_para[counter] <- ct_elap
  
  ## Iterated, parallelization
  # ptm <- proc.time()
  # side_seq_res <-
  #   SIDEseqIterated(expr_matrix = counts_subset, 
  #                   diff_expr_method = "diff_expr_norm",
  #                   similarity_method = "top_n_gene_intersect",
  #                   n_top_genes = 1000,
  #                   parallelize = TRUE,
  #                   n_cores = 6)
  # ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  # 
  # compute_time_df$iterate_para[counter] <- ct_elap
  
  ## updated counter 
  counter <- counter+1
  
}


## remove NAs, save
compute_time_df <- 
  compute_time_df %>% 
  mutate_all(.funs = ~ (ifelse(.==0, NA, .)))

#saveRDS(compute_time_df, here("output/SIDE_seq_test_runs.RDS"))


## Include Block computation results
#compute_time_df <- readRDS(here("output/SIDE_seq_test_runs.RDS"))

compute_time_df$block_para <- NA
  
counter <- 1
## compute block computation results
for(c in compute_time_df$num_cells) {
  print(c)
  ## subset to number of cells of interest
  counts_subset <- umi_counts_high_expr[,sample(ncol(umi_counts_high_expr), c)]
  
  ## track compute time
  ptm <- proc.time()
  side_seq_res <-
    SIDEseqBlockSimult(expr_matrix = counts_subset, 
                       diff_expr_method = "diff_expr_norm",
                       similarity_method = "top_n_gene_intersect",
                       n_top_genes = 1000,
                       block_memory = 50,
                       parallelize = TRUE,
                       n_cores = 6)
  
  ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  
  compute_time_df$block_para[counter] <- ct_elap
  
  counter <- counter + 1
}



### Plot results
compute_time_df %>%
  gather(method, time_secs, simult_para, iterate_para,
         iterate_random_subset_para) %>%
  filter_all(~ !is.na(.)) %>%
  ggplot(aes(x = num_cells, y = time_secs, color = method)) + 
  geom_point() + geom_line() + 
  theme_bw() + 
  labs(x = "Number of cells", y = "Run Time (sec)") + 
  ggtitle("SIDE-seq runtimes, number of cells\n Includes iterative random cell subset approach")


### Test runs of random sampling strategies -----------------------------------
#compute_time_df <- readRDS(here("output/SIDE_seq_test_runs.RDS"))

num_cells <- compute_time_df$num_cells
rand_iter_times <- rep(0, length(num_cells))
for(c in 9:length(num_cells)) {
  print(c)
  counts_subset <- 
    umi_counts_high_expr[,sample(ncol(umi_counts_high_expr), 
                                 num_cells[c])]
  
  ## set timer, run function
  ptm <- proc.time()
  small_iter_par_rand <-
    SIDEseqIteratedRand(expr_matrix = counts_subset, 
                        diff_expr_method = "diff_expr_norm",
                        similarity_method = "top_n_gene_intersect",
                        n_top_genes = 1000,
                        parallelize = TRUE)
  ## elapsed time
  ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  rand_iter_times[c] <- ct_elap
}

iter_times <- rep(0, length(num_cells))
for(c in seq_len(length(num_cells))) {
  print(c)
  counts_subset <- 
    umi_counts_high_expr[,sample(ncol(umi_counts_high_expr), 
                                 num_cells[c])]
  ptm <- proc.time()
  small_iter_par_rand <-
    SIDEseqIterated(expr_matrix = small_matrix, 
                        diff_expr_method = "diff_expr_norm",
                        similarity_method = "top_n_gene_intersect",
                        n_top_genes = 1000,
                        parallelize = TRUE)
  ct_elap <- (proc.time() - ptm)[3]  # compute time elapsed
  iter_times[c] <- ct_elap
}

### TODO: comparison of simulatenous against two code improvements ------------

### Object size comparisons ---------------------------------------------------


top_gene_options <- seq(100, 1000, 100)

num_cell_options <- seq(500, 3000, 250)
stop_gene_options <- c()

object_size_df <- 
  expand.grid(n_genes = top_gene_options, 
              n_cells = num_cell_options) %>% 
  mutate(object_size = 0)


## NOTE: this still needs to sample ranks, not runif for correct storage size.
for(i in 38:nrow(object_size_df)) {
  print(i)
  gc()
  if(!object_size_df[i,1] %in% stop_gene_options){
  #ERROR HANDLING
  possibleError <- tryCatch(
    matrix(runif(object_size_df[i,1] * dim(combn(object_size_df[i,2], 2))[2]), 
           nrow = object_size_df[i,1], 
           ncol = dim(combn(object_size_df[i,2], 2))[2]),
    error=function(e) e
  )
  if(!inherits(possibleError, "error")){
    rm(possibleError)
    gc()
    mat <- matrix(runif(object_size_df[i,1] * dim(combn(object_size_df[i,2], 2))[2]), 
                  nrow = object_size_df[i,1], 
                  ncol = dim(combn(object_size_df[i,2], 2))[2])
    ## megabytes
    object_size_df[i,3] = object.size(mat) * 1e-6
    rm(mat)
  } else{stop_gene_options <- 
    append(stop_gene_options, object_size_df[i,1])
  print(object_size_df[i,1])}
  
  }
} 

object_size_df <- 
  object_size_df %>% 
  mutate_all(.funs = ~ (ifelse(.==0, NA, .)))

#saveRDS(object_size_df, here("output/matrix_size_df.RDS"))

object_size_df %>% 
  filter_all(~ !is.na(.)) %>%
  ggplot(aes(x = n_cells, y = object_size, color = factor(n_genes))) + 
  geom_point() + geom_line() + 
  theme_bw() + 
  labs(x = "Number of cells (cols)", y = "Object Size (mb)", color = "Top Genes (rows)") + 
  ggtitle("Top Gene Matrix Storage Size for all Cell Pairings")



