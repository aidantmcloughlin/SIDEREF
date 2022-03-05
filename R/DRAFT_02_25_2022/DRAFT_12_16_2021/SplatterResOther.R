################################################################################
### Running Splatter simulations in functional manner
################################################################################

set.seed(291827)

n_neighbors <- 15
min_dist <- 0.01

source(here("R/PlotFunctions.R"))


###############################################################################
### Load assorted simulations and run plots.
###############################################################################



#### --------------------------------------------------------------------------
#### Simulation 2 -- Shared DE Groups
#### --------------------------------------------------------------------------

## load simulation data
load(here("output/splatter_sim2.RData"))

sim2 <- logNormCounts(sim2)
sim2_counts <- as.matrix(sim2@assays@data$logcounts)
## remove genes with no expression
sim2_counts <- sim2_counts[rowSums(sim2_counts) > 0, ]

splat_group_labels <- factor(sim2@colData@listData$Group)


levels(splat_group_labels) <- 
  c("1: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 2, 3)", 
    "2: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 1, 3)",
    "3: D.E. Prob 3.0%\nShared D.E. Prob 2.0% (with Grp. 1, 2)",
    "4: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 5)",
    "5: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 4)",
    "6: D.E. Prob 4.0%\nNo Shared D.E.")


plot_list_shared_de_123_56 <-
  runDimRedsAndPlot(sim_counts = sim2_counts, 
                    group_labels = splat_group_labels)


ggarrange(plot_list_shared_de_123_56[[1]], plot_list_shared_de_123_56[[4]], 
          plot_list_shared_de_123_56[[2]], plot_list_shared_de_123_56[[5]], 
          plot_list_shared_de_123_56[[3]], plot_list_shared_de_123_56[[6]],
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/splatter_shared_de.png"),
       width = 8.5, height = 11)



#### --------------------------------------------------------------------------
#### Simulation 3 - Shared DE Groups, Increased Noise
#### --------------------------------------------------------------------------

## load simulation data
load(here("output/splatter_sim3.RData"))


sim3 <- logNormCounts(sim3)
sim3_counts <- as.matrix(sim3@assays@data$logcounts)
## remove genes with no expression
sim3_counts <- sim3_counts[rowSums(sim3_counts) > 0, ]

splat_group_labels <- factor(sim3@colData@listData$Group)


levels(splat_group_labels) <- 
  c("1: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 2, 3)", 
    "2: D.E. Prob 1.0%\nShared D.E. Prob 2.0% (with Grp. 1, 3)",
    "3: D.E. Prob 3.0%\nShared D.E. Prob 2.0% (with Grp. 1, 2)",
    "4: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 5)",
    "5: D.E. Prob 3.0%\nShared D.E. Prob 3.0% (with Grp. 4)",
    "6: D.E. Prob 4.0%\nNo Shared D.E.")


plot_list_shared_de_123_56_noisy <-
  runDimRedsAndPlot(sim_counts = sim3_counts, 
                    group_labels = splat_group_labels)


ggarrange(plot_list_shared_de_123_56_noisy[[1]], plot_list_shared_de_123_56_noisy[[4]], 
          plot_list_shared_de_123_56_noisy[[2]], plot_list_shared_de_123_56_noisy[[5]], 
          plot_list_shared_de_123_56_noisy[[3]], plot_list_shared_de_123_56_noisy[[6]],
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "bottom")

ggsave(here("output/figures/splatter_noisy.png"),
       width = 8.5, height = 11)


