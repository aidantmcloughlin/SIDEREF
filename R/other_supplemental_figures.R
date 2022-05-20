library(here)
source(here("R/libraries.R"))



##### SNN Results -------------------------------------------------------------


SIDEREF_snn_plot <- 
  function(nn_consistency_df,
           snn_k = SNN_K) {
    
    ## get run time details from DF
    n_reps <- max(as.numeric(
      str_sub(names(nn_consistency_df), -1, -1)), na.rm=T)
    
    ref_set_sizes <- unique(nn_consistency_df$ref_set_size)
    
    nn_consistency_df_long <-
      nn_consistency_df %>% 
      dplyr::select(method, ref_set_size,
                    avg_snn, 
                    min_snn, max_snn) %>%
      tidyr::gather(measure, snn_stat,
                    avg_snn, 
                    min_snn, max_snn) %>% 
      mutate(method = 
               factor(case_when(measure == "avg_snn" ~ 
                                  method %p% " (Avg)",
                                measure == "min_snn" ~ 
                                  method %p% " (Min)",
                                measure == "max_snn" ~ 
                                  method %p% " (Max)")))
    
    nn_consistency_df_long$method <-
      forcats::fct_relevel(
        nn_consistency_df_long$method,
        c("Cell Type Stratified Sample (Max)", 
          "Cell Type Stratified Sample (Avg)", 
          "Cell Type Stratified Sample (Min)", 
          "Random Sample (Max)", "Random Sample (Avg)", "Random Sample (Min)"))
    
    p <-
      nn_consistency_df_long %>%
      ggplot(aes(x = ref_set_size, fill = method, y = snn_stat)) + 
      geom_bar(stat = "identity", position = "dodge", col = "black", 
               width=16) + #, alpha = 0.75) + 
      coord_cartesian(ylim = c(floor(min(nn_consistency_df_long$snn_stat)), 
                               snn_k)) +
      theme_bw() + 
      scale_x_continuous(name = "Reference Set Size", 
                         breaks = ref_set_sizes,
                         minor_breaks = NULL) +
      scale_fill_manual(values = c("dark blue", "blue", "light blue",
                                   "dark red", "red", "coral",
                                   "forest green", "green", "light green")) +
      labs(x = "Reference Set Size",
           y = "Mean Shared Nearest Neighbors (of " %p% snn_k %p% ")",
           fill = "Method (Statistic over " %p% n_reps %p% " Runs)")
    
    return(p)
    
  }


### load, plot

load(here("output/computations/" %p%
            "snn_SIDEREF_SIDEseq.RData"),
     verbose = TRUE)

SIDEREF_snn_plot(nn_consistency,
                 snn_k = SNN_K) + 
  theme(axis.title = element_text(size=13),
        legend.title = element_text(size=13))

ggsave(here("manuscript_files/FigureS1.eps"),
       plot = last_plot(),
       width = 12, height = 7,
       device='eps', dpi=DPI)


##### Compute Time ------------------------------------------------------------


SIDEREF_SIDEseq_compute_time_plot <- 
  function(compute_time_df,
           n_reps,
           log_scale=FALSE) {
    n_cells <- unique(compute_time_df$n_cells)
    
    compute_time_df <-
      compute_time_df %>% 
      gather(method, run_time, 
             SIDEREF_comp_time, SIDEseq_comp_time) %>% 
      mutate(method = gsub("_comp_time", "", method),
             method = ifelse(method == "SIDEREF",
                             method %p% " (" %p% R %p% " Ref. Set Size)",
                             method))
    
    compute_time_df <-
      compute_time_df %>% 
      mutate(log_n_cells = log(n_cells),
             log_run_time = log(run_time))
    
    
    
    p1 <- ggplot(data = compute_time_df) + 
      geom_line(aes(x = n_cells, 
                    y = run_time,
                    colour = method)) +
      geom_point(aes(x = n_cells, 
                     y = run_time,
                     colour = method)) +
      theme_bw()
    
    p1 <- p1 +
      scale_x_continuous(name = "Number of Cells", 
                         breaks= n_cells,
                         minor_breaks = NULL) + 
      labs(x = "Number of Cells",
           y = "Compute Time, Minutes",
           colour = "Method")
    
    ## LM for power rule.
    m1 <- lm(log_run_time ~ log_n_cells, data=compute_time_df %>% 
               dplyr::filter(method == "SIDEREF (100 Ref. Set Size)"))$coefficients %>% unname()
    m2 <- lm(log_run_time ~ log_n_cells, data=compute_time_df %>% 
               dplyr::filter(method == "SIDEseq"))$coefficient %>% unname()
    
    m1[1] <- exp(m1[1])
    m2[1] <- exp(m2[1])
    
    p2 <- ggplot(data = compute_time_df) + 
      geom_line(aes(x = log_n_cells, 
                    y = log_run_time,
                    colour = method)) +
      geom_point(aes(x = log_n_cells, 
                     y = log_run_time,
                     colour = method)) +
      theme_bw() + 
      labs(x = "ln(Number of Cells)",
           y = "ln(Compute Time, Minutes)")
    
    label_data <- 
      data.frame(
        y = c(6, 5.5),
        x = 6.25,
        slope = "OLS Slope: " %p% as.character(round(c(m2[2], m1[2]), 3)),
        method = c("SIDEseq", "SIDEREF (100 Ref. Set Size)"))
    
    p2 <- p2 + 
      geom_text(data = label_data, aes(x=x,y=y,label=slope,colour=method)) + 
      labs(colour = "Method")
    
    p <- ggarrange(p1, p2, ncol=2, common.legend = TRUE)
    
    
    return(p)
  }



### load, plot
load(here("output/computations/" %p%
            "compute_time_SIDEREF_SIDEseq.RData"),
     verbose = TRUE)

SIDEREF_SIDEseq_compute_time_plot(compute_time_df,
                                  n_reps = COMPUTE_TIME_REPS,
                                  log_scale=TRUE)


ggsave(here("manuscript_files/FigureS2.eps"),
       plot = last_plot(),
       width = 10, height = 5,
       device='eps', dpi=DPI)





##### Converge Results --------------------------------------------------------

load(here("output/computations/sideref_converge_df_list.RData"),
     verbose=TRUE)

### average over runs
AVG_RUNS <- TRUE

if(AVG_RUNS) {
  nonzero_reps <- sum(sapply(res_df_list, function(x) x$norm_diff[nrow(x)] > 0))
  
  avg_res_df <- Reduce("+", res_df_list[1:nonzero_reps]) / nonzero_reps  
} else{
  avg_res_df <- res_df_list[[1]]
}

last_entry <- max(which(avg_res_df$norm_diff > 0))
avg_res_df <- avg_res_df[1:last_entry, ]

### Plot Res
p_converge <-
  avg_res_df %>%
  mutate(log_10_Fnorm = log(norm_diff, base=10) )%>% 
  gather(cat, val, quantile_90:log_10_Fnorm) %>%
  mutate(cat = case_when(
    cat == "log_10_Fnorm" ~ "Log (Base 10) of Frobenius Norm",
    cat == "max_diff" ~ "Maximum",
    cat == "quantile_90" ~ "90th Quantile",
    cat == "quantile_95" ~ "95th Quantile",
    cat == "quantile_97.5" ~ "97.5th Quantile",
    cat == "quantile_99" ~ "99th Quantile",
    TRUE ~ cat)
  ) %>%
  mutate(cat = factor(
    cat, 
    levels = c("Log (Base 10) of Frobenius Norm",
               "Maximum",
               "90th Quantile",
               "95th Quantile",
               "97.5th Quantile",
               "99th Quantile"))) %>%
  ggplot(aes(x=n_prev, y=val, color=cat)) +
  geom_line() + 
  geom_point(size=0.5) +
  scale_x_continuous(name = "Sample Size per Cell Type (n1, n2)", 
                     breaks= unique(avg_res_df$n_prev),
                     labels= '(' %p% unique(avg_res_df$n_prev) %p% ',' %p% 
                       unique(avg_res_df$n) %p% ')',
                     minor_breaks = NULL) + 
  theme_bw() + 
  labs(y="",color=expression(paste("Stat of: |", D[n2]-D[n1], "|"))) + 
  theme(axis.title = element_text(size=13),
        legend.title = element_text(size=13))


ggsave(here("manuscript_files/FigureS5.eps"),
       plot = p_converge,
       width = 15, height = 6,
       device='eps', dpi=DPI)



