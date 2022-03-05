## scratchpad




load(file=here("output/hclust/side_50_hclust.RData"))

side_50_full_hclust_res$adj_rand_index

load(file=here("output/hclust/side_150_hclust.RData"))
side_150_full_hclust_res$adj_rand_index

load(file=here("output/hclust/pca_25_hclust.RData"))
pca_25_full_hclust_res$adj_rand_index


load(file= here("output/hclust/pca_3_hclust.RData"))
pca_3_full_hclust_res$adj_rand_index

pca_3_full_hclust_res$p


load(file=here("output/hclust/hclust_full_splat.RData"))


hclust_res_df$side_ref_g50_dist$p

ari_df <-
  data.frame(dist_metric = names(hclust_res_df),
             ari = sapply(hclust_res_df, function(x) x$adj_rand_index))


rownames(ari_df)=NULL

ari_df %>% 
  fwrite(file=here("output/hclust/hclust_ari_low_var_splat.csv"))



### no low clustering.
load(file=here("output/hclust/hclust_full_splat_no_low_clust.RData"))

ari_df <-
  data.frame(dist_metric = names(hclust_res_df),
             ari = sapply(hclust_res_df, function(x) x$adj_rand_index))


rownames(ari_df)=NULL

ari_df %>% 
  fwrite(file=here("output/hclust/hclust_ari_low_var_splat_no_low_clust.csv"))





###############################################################################
###############################################################################
higvar_name = "splatter_sim_skin_eq_grps_highvar"

load(file=here("output/hclust/hclust_" %p% higvar_name %p% ".RData"))

ari_df <-
  data.frame(dist_metric = names(hclust_res_df),
             ari = sapply(hclust_res_df, function(x) x$adj_rand_index))


rownames(ari_df)=NULL

ari_df %>% 
  fwrite(file=here("output/hclust/hclust_ari_high_var_splat.csv"))



### no low clustering.
load(file=here("output/hclust/hclust_" %p% higvar_name %p% "_no_low_clust.RData"))

ari_df <-
  data.frame(dist_metric = names(hclust_res_df),
             ari = sapply(hclust_res_df, function(x) x$adj_rand_index))


rownames(ari_df)=NULL

ari_df %>% 
  fwrite(file=here("output/hclust/hclust_ari_high_var_splat_no_low_clust.csv"))


###############################################################################
###############################################################################

## Function for measuring group-wise clust performance


load(file=here("output/hclust/hclust_splatter_sim_skin_eq_grps_highvar_no_low_clust.RData"))



subgroup_metric_df <-
  do.call(rbind,
          lapply(seq_len(length(hclust_res_df)), function(x) {
            df = hclust_res_df[[x]]$hclust_res$clust_res_df
            
            ## relevel cell_group_high
            df <- 
              df %>% 
              group_by(cell_group_high) %>% 
              mutate(min_cell_group_low = min(cell_group_low),
                     n_cell_groups = length(unique(cell_group_low))) %>%
              ungroup() %>%
              dplyr::filter(n_cell_groups > 1)
            
            min_cells = data.frame(min_cell_group_low=unique(df$min_cell_group_low)) %>% 
              mutate(new_subgroup = rank(min_cell_group_low))
            
            df <- 
              df %>%
              left_join(min_cells) %>% 
              mutate(cell_group_high = new_subgroup) %>%
              group_by(cell_group_high) %>%
              summarise(
                group_variance = var(cluster_high)) %>% 
              ungroup() %>%
              mutate(method = names(hclust_res_df)[x])
            
          })) %>% 
  left_join(name_xwalk) %>% 
  mutate(ggplot_name = factor(ggplot_name, levels = method_levels))

subgroup_metric_df %>% 
  mutate(cell_group_high = as.character(cell_group_high)) %>%
  dplyr::filter(!grepl("spectral", ggplot_name, ignore.case = TRUE)) %>%
  ggplot(aes(x = cell_group_high, 
             y = group_variance, 
             fill = ggplot_name)) + 
  geom_bar(stat = "identity", position="dodge") + 
  labs(x="Actual Subgroup", y="Group Variance", fill = "Method") + 
  theme_bw() + 
  theme(legend.position="bottom")
  

subgroup_metric_df %>% 
  mutate(cell_group_high = as.character(cell_group_high)) %>%
  dplyr::filter(grepl("spectral", ggplot_name, ignore.case = TRUE)) %>%
  ggplot(aes(x = cell_group_high, 
             y = group_variance, 
             fill = ggplot_name)) + 
  geom_bar(stat = "identity", position="dodge") + 
  labs(x="Actual Subgroup", y="Group Variance", fill = "Method") + 
  theme_bw() + 
  theme(legend.position="bottom")



