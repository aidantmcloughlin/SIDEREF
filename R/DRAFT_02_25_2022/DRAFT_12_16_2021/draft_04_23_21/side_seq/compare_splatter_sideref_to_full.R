### Compare 1500 full side seq to subset SIDEseq runs.

getKNNList <- function(dist_mat, k) {
  knn_list <- 
    lapply(seq_len(dim(dist_mat)[1]),
           function(i) {
             return(order(dist_mat[,i])[-1][1:k])
           })
}


load(here("output/splatter_sideseq/splatter_sideseq_full.RData"))


knn_full_splatter <- getKNNList(splatter_sideseq_full_res, k = 15)


rand_res <- list()
embed_5_res <- list()
embed_10_res <- list()

methods <- c("Cell Subgroups (5)",
             "Cell Subgroups (10)",
             "Random Cells")



ref_set_size <- c(5, 10, 25, 50, 100, 150)


nn_consistency <- 
  as.data.frame(
    expand.grid(method = methods, 
                ref_set_size = ref_set_size),
    stringsAsFactors = FALSE) %>% 
  mutate(avg_nn_consistency_run1 = NA,
         avg_nn_consistency_run2 = NA,
         avg_nn_consistency_run3 = NA,
         avg_nn_consistency_run4 = NA,
         avg_nn_consistency_run5 = NA)

i <- 0
for(r in ref_set_size) {
  print(r)
  for(s in 1:5) {
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "cell_embed_ref_set_splatter.RData"))
    side_seq_ref_cell_embed_5  <- side_seq_ref_cell_embed
    if(r >= 10) {
      load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                  "cell_embed_ref_set_splatter_10_clust.RData"))
      side_seq_ref_cell_embed_10 <- side_seq_ref_cell_embed  
    }
    
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "rand_ref_set_splatter.RData"))
    
    
    
    knn_embed_5_splatter  <- getKNNList(side_seq_ref_cell_embed_5$dissim_final,  k = 15)
    
    if(r >= 10) {
      knn_embed_10_splatter <- getKNNList(side_seq_ref_cell_embed_10$dissim_final, k = 15)}
    
    knn_rand_splatter     <- getKNNList(side_seq_ref_rand$dissim_final,  k = 15)
    
    
    neighbor_intersections_embed_5 <- 
      sapply(seq_len(nrow(splatter_sideseq_full_res)), 
             function(i){length(intersect(knn_embed_5_splatter[[i]],
                                          knn_full_splatter[[i]]))})
    if(r >= 10) {
      neighbor_intersections_embed_10 <- 
        sapply(seq_len(nrow(splatter_sideseq_full_res)), 
               function(i){length(intersect(knn_embed_10_splatter[[i]],
                                            knn_full_splatter[[i]]))}) 
    }
    
    neighbor_intersections_rand <- 
      sapply(seq_len(nrow(splatter_sideseq_full_res)), 
             function(i){length(intersect(knn_rand_splatter[[i]],
                                          knn_full_splatter[[i]]))})
    

    avg_embed_5  <- mean(neighbor_intersections_embed_5)
    if(r >= 10) {
      avg_embed_10 <- mean(neighbor_intersections_embed_10)} else{avg_embed_10 <- NA}
    avg_rand     <- mean(neighbor_intersections_rand)

    
    nn_consistency[(3*i+1), (2+s)] = avg_embed_5
    nn_consistency[(3*i+2), (2+s)] = avg_embed_10
    nn_consistency[(3*i+3), (2+s)] = avg_rand
    

    
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
           avg_nn_consistency_run5))


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


## reorder levels
nn_consistency_long$method <-
  forcats::fct_relevel(
    nn_consistency_long$method,
    c("Cell Subgroups (5) (Max)", "Cell Subgroups (5) (Avg)", "Cell Subgroups (5) (Min)",
      "Cell Subgroups (10) (Max)", "Cell Subgroups (10) (Avg)", "Cell Subgroups (10) (Min)",
      "Random Cells (Max)", "Random Cells (Avg)", "Random Cells (Min)"))

## plotting it all
nn_consistency_long %>%
  ggplot(aes(x = factor(ref_set_size), fill = method, y = nn_consistency)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black", 
           width=0.5, alpha = 0.75) + 
  theme_bw() + 
  scale_fill_manual(values = c("dark blue", "blue", "light blue",
                               "dark red", "red", "coral",
                               "darkgreen", "forestgreen", "darkolivegreen1")) +
  labs(x = "Reference Set Size",
       y = "Mean Nearest Neighbor Intersection (of 15)",
       fill = "Method (5 Runs Statistic)")

ggsave(file = here("output/figures/splatter_sideref_consistency.png"),
       height = 6, width = 9)




###############################################################################
### Embedded method purity trends
###############################################################################

purity_data_5 <- data.frame(run_info = rep("", length(ref_set_size) * 5),
                            purity_2 = rep(0, length(ref_set_size) * 5),
                            purity_3 = rep(0, length(ref_set_size) * 5))

purity_data_10 <- data.frame(run_info = rep("", length(ref_set_size) * 5),
                             purity_2 = rep(0, length(ref_set_size) * 5),
                             purity_3 = rep(0, length(ref_set_size) * 5))

i <- 1
for(r in ref_set_size) {
  print(r)
  for(s in 1:5) {
    load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                "cell_embed_ref_set_splatter.RData"))
    
    purity_data_5[i,] <- c(r %p% "_" %p% s, side_seq_ref_cell_embed$purity_scores[-1])
    
    
    if(r >= 10) {
      load(here("output/splatter_sideseq/run" %p% s %p% "_" %p% r %p% 
                  "cell_embed_ref_set_splatter_10_clust.RData"))
      purity_data_10[i,] <- c(r %p% "_" %p% s, side_seq_ref_cell_embed$purity_scores[-1])
    }
    
    
    
    i <- i + 1
  }
}


purity_data_clean_5 <- purity_data_5 %>% 
  tidyr::separate(run_info, into = c("ref_set_size", "run"), sep = "_")
purity_data_clean_10 <- purity_data_10 %>% 
  tidyr::separate(run_info, into = c("ref_set_size", "run"), sep = "_")


purity_data_clean_5 %>% 
  tidyr::gather(iteration, purity, purity_2, purity_3) %>%
  mutate(purity = as.numeric(purity),
         iteration = ifelse(iteration == "purity_2", 2, 3)) %>%
  dplyr::filter(run == 1) %>% 
  ggplot(aes(color = factor(ref_set_size, levels = c("25", "50", "100", "200")), 
             lty = run, 
             x = iteration, y = purity)) + 
  geom_line() + 
  coord_cartesian(xlim = c(2, 3)) +
  theme_bw() + scale_x_continuous(name = "Iteration", breaks=c(2, 3),
                                  minor_breaks = NULL) + 
  labs(y = "Purity of Current Clusters with Prior Clusters",
       lty = "Run Number", col = "Ref Set Size") + 
  scale_color_manual(values = c("red", "blue", "orange", "dark green"))

ggsave(here("output/figures/purity_sideref_splatter_5_clust.png"),
       width = 8, height = 5)




purity_data_clean_10 %>% 
  tidyr::gather(iteration, purity, purity_2, purity_3) %>%
  mutate(purity = as.numeric(purity),
         iteration = ifelse(iteration == "purity_2", 2, 3)) %>%
  dplyr::filter(ref_set_size %in% c(25, 50, 100, 200)) %>%
  ggplot(aes(color = factor(ref_set_size, levels = c("25", "50", "100", "200")), 
             lty = run, 
             x = iteration, y = purity)) + 
  geom_line() + 
  coord_cartesian(xlim = c(2, 3)) +
  theme_bw() + scale_x_continuous(name = "Iteration", breaks=c(2, 3),
                                  minor_breaks = NULL) + 
  labs(y = "Purity of Current Clusters with Prior Clusters",
       lty = "Run Number", col = "Ref Set Size") + 
  scale_color_manual(values = c("red", "blue", "orange", "dark green"))

ggsave(here("output/figures/purity_sideref_splatter_10_clust.png"),
       width = 8, height = 5)



### Extended purity for many iterations
for(r in c(10, 25, 75)) {
  print(r)
  test_purity_10_clust <- 
    SIDEseqRefSet(expr_matrix = sim1_counts, 
                  diff_expr_method = "diff_expr_norm",
                  similarity_method = "n_intersect",
                  n_top_genes = def_n_top_genes,
                  ## cell reference set selection parameters
                  selection_method = "cell_embed_sample",
                  size_ref_set = r,
                  n_clust = 10,
                  B = 10,
                  D = 1,
                  purity_tol = 1,
                  ## parallelization parameters
                  parallelize = TRUE,
                  n_cores = 5,
                  ## other
                  verbose = TRUE)
  
  save(test_purity_10_clust, file = here("output/splatter_sideseq/test_purity_" %p% r %p% ".RData"))
}






