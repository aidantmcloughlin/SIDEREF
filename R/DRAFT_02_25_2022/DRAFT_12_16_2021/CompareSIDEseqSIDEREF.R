###############################################################################
### Compare full side seq to SIDEREF runs.
###############################################################################

set.seed(1)

## Commonly used function:

getKNNList <- function(dist_mat, k) {
  knn_list <- 
    lapply(seq_len(dim(dist_mat)[1]),
           function(i) {
             return(order(dist_mat[,i])[-1][1:k])
           })
}



#### --------------------------------------------------------------------------
#### Analysis on SNAREseq Data
#### --------------------------------------------------------------------------


load(here("output/snare_sideseq/side_seq_snare_1500.RData"))


knn_full_snare <- getKNNList(full_side_seq, k = 15)


rand_res <- list()
embed_res <- list()
pca_embed_res <- list()


methods <- c("Cell Subgroup Strata (Iterative SIDEREF)", 
             "Cell Subgroup Strata (PCA, UMAP)",
             "Random Cells")


ref_set_size <- c(25, 50, 75, 100, 150, 200, 300, 400)



rand_table      <- list()
embed_table     <- list()
pca_embed_table <- list()

nn_consistency <- 
  as.data.frame(
    expand.grid(method = methods, 
                ref_set_size = ref_set_size),
    stringsAsFactors = FALSE) %>% 
  mutate(avg_nn_consistency_run1 = NA,
         avg_nn_consistency_run2 = NA,
         avg_nn_consistency_run3 = NA,
         avg_nn_consistency_run4 = NA,
         avg_nn_consistency_run5 = NA,
         compute_time_run1 = NA,
         compute_time_run2 = NA,
         compute_time_run3 = NA,
         compute_time_run4 = NA,
         compute_time_run5 = NA)

i <- 0
for(r in ref_set_size) {
  print(r)
  for(s in 1:5) {
    load(here("output/snare_sideseq/run" %p% s %p% "_" %p% r %p% "cell_embed_ref_set_snare_1500.RData"))
    load(here("output/snare_sideseq/run" %p% s %p% "_" %p% r %p% "cell_embed_pca_umap_ref_set_snare_1500.RData"))
    load(here("output/snare_sideseq/run" %p% s %p% "_" %p% r %p% "rand_ref_set_snare_1500.RData"))
    
    
    
    knn_embed_snare     <- getKNNList(side_seq_ref_cell_embed$dissim_final, k = 15)
    knn_pca_embed_snare <- getKNNList(side_seq_ref_cell_embed_pca_umap$dissim_final, k = 15)
    knn_rand_snare      <- getKNNList(side_seq_ref_rand$dissim_final,  k = 15)
    
    
    neighbor_intersections_embed <- 
      sapply(seq_len(nrow(full_side_seq)), 
             function(i){length(intersect(knn_embed_snare[[i]],
                                          knn_full_snare[[i]]))})
    
    neighbor_intersections_pca_embed <- 
      sapply(seq_len(nrow(full_side_seq)), 
             function(i){length(intersect(knn_pca_embed_snare[[i]],
                                          knn_full_snare[[i]]))})
    
    neighbor_intersections_rand <- 
      sapply(seq_len(nrow(full_side_seq)), 
             function(i){length(intersect(knn_rand_snare[[i]],
                                          knn_full_snare[[i]]))})
    
    table_embed     <- table(neighbor_intersections_embed)
    table_pca_embed <- table(neighbor_intersections_pca_embed)
    table_rand      <- table(neighbor_intersections_rand)
    
    avg_embed     <- mean(neighbor_intersections_embed)
    avg_pca_embed <- mean(neighbor_intersections_pca_embed)
    avg_rand      <- mean(neighbor_intersections_rand)
    
    embed_table[[i+1]]     <- table_embed
    pca_embed_table[[i+1]] <- table_pca_embed
    rand_table[[i+1]]      <- table_rand
    
    
    nn_consistency[(3*i+1), (2+s)] = avg_embed
    nn_consistency[(3*i+2), (2+s)] = avg_pca_embed
    nn_consistency[(3*i+3), (2+s)] = avg_rand
    
    nn_consistency[(3*i+1), (7+s)] = side_seq_ref_cell_embed[[3]]
    nn_consistency[(3*i+2), (7+s)] = side_seq_ref_cell_embed_pca_umap[[3]]
    nn_consistency[(3*i+3), (7+s)] = side_seq_ref_rand[[3]]
    
  }
  i <- i + 1 
}


## summarizing:
nn_consistency <-
  nn_consistency %>%
  mutate(
    ## NN Consistency
    avg_nn_consistency = 
      0.2*(avg_nn_consistency_run1 +
             avg_nn_consistency_run2 +
             avg_nn_consistency_run3 +
             avg_nn_consistency_run4 +
             avg_nn_consistency_run5),
    min_nn_consistency = 
      pmin(avg_nn_consistency_run1,
           avg_nn_consistency_run2,
           avg_nn_consistency_run3,
           avg_nn_consistency_run4,
           avg_nn_consistency_run5),
    max_nn_consistency = 
      pmax(avg_nn_consistency_run1,
           avg_nn_consistency_run2,
           avg_nn_consistency_run3,
           avg_nn_consistency_run4,
           avg_nn_consistency_run5),
    ## Compute Time
    avg_compute_time = 
      0.2*(compute_time_run1 +
             compute_time_run2 +
             compute_time_run3 +
             compute_time_run4 +
             compute_time_run5),
    min_compute_time = 
      pmin(compute_time_run1,
           compute_time_run2,
           compute_time_run3,
           compute_time_run4,
           compute_time_run5),
    max_compute_time = 
      pmax(compute_time_run1,
           compute_time_run2,
           compute_time_run3,
           compute_time_run4,
           compute_time_run5))


nn_consistency_long <- 
  nn_consistency %>% 
  dplyr::select(method, ref_set_size,
                avg_nn_consistency, min_nn_consistency, max_nn_consistency) %>%
  tidyr::gather(measure, nn_consistency,
                avg_nn_consistency, min_nn_consistency, max_nn_consistency) %>% 
  mutate(method = 
           factor(case_when(measure == "avg_nn_consistency" ~ method %p% " (Avg)",
                            measure == "min_nn_consistency" ~ method %p% " (Min)",
                            measure == "max_nn_consistency" ~ method %p% " (Max)")))



nn_consistency_long$method <-
  forcats::fct_relevel(
    nn_consistency_long$method,
    c("Cell Subgroup Strata (Iterative SIDEREF) (Max)", 
      "Cell Subgroup Strata (Iterative SIDEREF) (Avg)", 
      "Cell Subgroup Strata (Iterative SIDEREF) (Min)",
      "Cell Subgroup Strata (PCA, UMAP) (Max)", 
      "Cell Subgroup Strata (PCA, UMAP) (Avg)", 
      "Cell Subgroup Strata (PCA, UMAP) (Min)",
      "Random Cells (Max)", "Random Cells (Avg)", "Random Cells (Min)"))

## plotting it all
nn_consistency_long %>%
  ggplot(aes(x = ref_set_size, fill = method, y = nn_consistency)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black", 
           width=16, alpha = 0.75) + 
  coord_cartesian(ylim=c(12.5,15)) +
  theme_bw() + 
  scale_x_continuous(name = "Reference Set Size", 
                     breaks=c(25, 50, 75, 100,
                              150, 200, 300, 400),
                     minor_breaks = NULL) +
  scale_fill_manual(values = c("dark blue", "blue", "light blue",
                               "dark red", "red", "coral",
                               "forest green", "green", "light green")) +
  labs(x = "Reference Set Size",
       y = "Mean Shared Nearest Neighbors (of 15)",
       fill = "Method (5 Runs Statistic)")


ggsave(file = here("output/figures/snare_sideref_consistency.png"),
       height = 6, width = 12)





## compute time:
compute_time_long <- 
  nn_consistency %>% 
  dplyr::select(method, ref_set_size,
                avg_compute_time, min_compute_time, max_compute_time) %>%
  tidyr::gather(measure, compute_time,
                avg_compute_time, min_compute_time, max_compute_time) %>% 
  mutate(method = 
           factor(case_when(measure == "avg_compute_time" ~ method %p% " (Avg)",
                            measure == "min_compute_time" ~ method %p% " (Min)",
                            measure == "max_compute_time" ~ method %p% " (Max)")))



compute_time_long$method <-
  forcats::fct_relevel(
    nn_consistency_long$method,
    c("Cell Subgroup Strata (Iterative SIDEREF) (Max)", 
      "Cell Subgroup Strata (Iterative SIDEREF) (Avg)", 
      "Cell Subgroup Strata (Iterative SIDEREF) (Min)",
      "Cell Subgroup Strata (PCA, UMAP) (Max)", 
      "Cell Subgroup Strata (PCA, UMAP) (Avg)", 
      "Cell Subgroup Strata (PCA, UMAP) (Min)",
      "Random Cells (Max)", "Random Cells (Avg)", "Random Cells (Min)"))

## linear fit?
summary_subgroup_glm <-
  summary(
    glm(compute_time ~ ref_set_size,
        data = compute_time_long %>% dplyr::filter(method == "Cell Subgroup Strata (Iterative SIDEREF) (Avg)"))
    )

summary_subgroup_glm$coefficients[2,1]

summary_rand_glm <-
  summary(
    glm(compute_time ~ ref_set_size,
        data = compute_time_long %>% dplyr::filter(method == "Random Cells (Avg)"))
    )

## plotting it all
compute_time_long %>%
  ggplot(aes(x = ref_set_size, fill = method, y = compute_time/60)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black", 
           width=16, alpha = 0.75) + 
  geom_abline(slope = summary_rand_glm$coefficients[2,1] / 60, 
              intercept = summary_rand_glm$coefficients[1,1] / 60, 
              col = "black", alpha = 0.5, lty = 2) +
  geom_abline(slope = summary_subgroup_glm$coefficients[2,1] / 60, 
              intercept = summary_subgroup_glm$coefficients[1,1] / 60, 
              col = "black", alpha = 0.5, lty = 2) +
  theme_bw() + 
  scale_x_continuous(name = "Reference Set Size", 
                     breaks=c(25, 50, 75, 100,
                              150, 200, 300, 400),
                     minor_breaks = NULL) +
  scale_fill_manual(values = c("dark blue", "blue", "light blue",
                               "dark red", "red", "coral",
                               "forest green", "green", "light green")) +
  labs(x = "Reference Set Size",
       y = "Compute Time (Mins)",
       fill = "Method (5 Runs Statistic)")


ggsave(file = here("output/figures/snare_sideref_compute_time.png"),
       height = 6, width = 12)




