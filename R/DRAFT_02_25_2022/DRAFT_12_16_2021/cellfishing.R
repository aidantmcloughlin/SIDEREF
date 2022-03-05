###############################################################################
###
### CellFishing Experiments
### Aidan McLoughlin
### 
###############################################################################


## Constants
SPECTRAL_DIMS <- 15
celltype_keyword <- "microglial cell"
file_name <- "tab_muris_full_0.05pct"

## Packages, Modules
library(here)
source(here('R/libraries.R'))
library(GeneFishing)


### load Tabula Muris data
load(here('output/tab_muris_sc/' %p% file_name %p% '_scale_var_genes_3000.RData'))
load(here('output/tab_muris_sc/' %p% file_name %p% '.RData'))



## Example 1: select cell type
table(meta_data_full$cell_type)


cells_of_int <- 
  colnames(data_full_variable)[which(meta_data_full$cell_type == celltype_keyword)]



probe_res <-
  GeneFishing::probeFishability(data_full_variable %>% t(),
                                potential_bait = cells_of_int,
                                n_rounds = 50,
                                ncores = 5)


doParallel::stopImplicitCluster()


gf_res <-
  GeneFishing::geneFishing(data_full_variable %>% t(),
                           bait_genes = probe_res$best_bait,
                           fishing_rounds = 50,
                           n_probing_rounds = 50)

doParallel::stopImplicitCluster()

save(probe_res, gf_res, file = here("output/tab_muris_sc/gf_res" %p% "microglial.RData"))

### Plot cells WRT genefishing status
# (we stick with Spearman => Spectral For now)
spearman_dist <- 1 - abs(cor(data_full_variable, method = "spearman"))

spearman_dist_spectral <- 
  dist(spectralEmbed(spearman_dist, ndim = SPECTRAL_DIMS))

save(spearman_dist_spectral, spearman_dist, 
     file = here("output/tab_muris_sc/dist_res/spearman_dist_0.05pct.RData"))

umap_embed <- 
  uwot::umap(as.dist(spearman_dist_spectral))


umap_embed_df <- 
  umap_embed %>% 
  data.frame() %>% 
  mutate(cell_names = colnames(data_full_variable)) %>% 
  mutate(cell_type = meta_data_full$cell_type,
         tissue = meta_data_full$tissue) %>% 
  mutate(fish_status = 
           case_when(cell_names %in% probe_res$best_bait ~ "Bait Cells",
                     cell_names %in% gf_res$fished_genes ~ "Fished Cells",
                     TRUE ~ "Other Cells"))



p1 <- 
  umap_embed_df %>% 
  mutate(fish_status=factor(fish_status)) %>% 
  ggplot(aes(x=X1, y=X2, color = fish_status)) + 
  geom_point(alpha = 1) +
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", col="Fished Cell") + 
  ggtitle(str_to_title(celltype_keyword) %p% " Bait") + 
  scale_color_manual(values = c("blue", "orange", "dark red")) + 
  theme(legend.position="bottom")

p1


p2 <- 
  umap_embed_df %>% 
  mutate(fish_status=factor(fish_status)) %>% 
  dplyr::filter(fish_status != "Other Cells") %>% 
  ggplot(aes(x=X1, y=X2, color = fish_status)) + 
  geom_point(alpha = 1) +
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", col="Fished Cell") + 
  ggtitle("Fished Cells Only") + 
  scale_color_manual(values = c("blue", "orange", "dark red")) + 
  theme(legend.position="bottom")

p2

p3 <- 
  umap_embed_df %>% 
  mutate(fish_status=factor(fish_status)) %>% 
  group_by(cell_type) %>% 
  mutate(num_cells = n()) %>% 
  ungroup() %>% 
  dplyr::filter(fish_status == "Fished Cells") %>% 
  group_by(cell_type) %>% 
  mutate(cell_type = cell_type %p% " (" %p% n() %p% " of " %p% num_cells %p% ")") %>% 
  ungroup() %>% 
  ggplot(aes(x=X1, y=X2, color = cell_type)) + 
  geom_point(alpha = 1) +
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", col="Cell Type") + 
  ggtitle("Fished Cells Cell Types") + 
  theme(legend.position="bottom")

p3

p4 <- 
  umap_embed_df %>% 
  mutate(fish_status=factor(fish_status)) %>% 
  group_by(tissue) %>% 
  mutate(num_cells = n()) %>% 
  ungroup() %>% 
  dplyr::filter(fish_status == "Fished Cells") %>% 
  group_by(tissue) %>% 
  mutate(tissue = tissue %p% " (" %p% n() %p% " of " %p% num_cells %p% ")") %>% 
  ungroup() %>% 
  ggplot(aes(x=X1, y=X2, color = tissue)) + 
  geom_point(alpha = 1) +
  theme_bw() + 
  labs(x="UMAP1", y="UMAP2", col="") + 
  ggtitle("Fished Cells Tissues") + 
  theme(legend.position="bottom")

p4


fish_res_plot_h = 6
fish_res_plot_w = 6


ggsave(here("output/tab_muris_sc/microglial_fished_all.png"), 
       plot = p1,
       height = fish_res_plot_h-1, 
       width = fish_res_plot_w-1)

ggsave(here("output/tab_muris_sc/microglial_fished_only.png"), 
       plot = p2,
       height = fish_res_plot_h-1,
       width = fish_res_plot_w-1)

ggsave(here("output/tab_muris_sc/microglial_fished_celltypes.png"), 
       plot = p3 + theme(legend.position = "right"),
       height = fish_res_plot_h,
       width = fish_res_plot_w+4)

ggsave(here("output/tab_muris_sc/microglial_fished_tissues.png"), 
       plot = p4 + theme(legend.position = "right"),
       height = fish_res_plot_h-1,
       width = fish_res_plot_w+2)




###############################################################################
###   Sample Spectral Embedding
###############################################################################

set.seed(72621)

p_list = vector(mode = "list", length=4)


for(p in 1:length(p_list)) {
  fished_cells <- gf_res$fished_genes
  
  bait_cells <- probe_res$best_bait
  
  pool_cells <- sample(setdiff(colnames(data_full_variable), bait_cells), 
                       length(bait_cells) * 5)
  
  
  
  spearman_dist_pool <- 1 - abs(cor(data_full_variable[, c(bait_cells, pool_cells)], 
                                    method = "spearman"))
  
  spearman_dist_spectral_pool <- 
    dist(spectralEmbed(spearman_dist_pool, ndim = SPECTRAL_DIMS))
  
  umap_embed_pool <- 
    uwot::umap(as.dist(spearman_dist_spectral_pool))
  
  p_list[[p]] <-
    umap_embed_pool %>% 
    data.frame() %>% 
    mutate(cell = c(bait_cells, pool_cells)) %>% 
    mutate(cell_status = 
             case_when(cell %in% bait_cells ~ "Bait",
                       cell %in% cells_of_int ~ "Other Microbial",
                       cell %in% fished_cells ~ "Fished", 
                       TRUE ~ "Other")) %>% 
    ggplot(aes(x=X1, y=X2)) + 
    geom_point(aes(col=cell_status)) + 
    labs(x="UMAP1", y="UMAP2") + 
    theme_bw() + 
    scale_color_manual(values=c("red", "blue", "orange", "forest green"))
  
}

ggarrange(p_list[[1]], p_list[[2]],
          p_list[[3]], p_list[[4]],
          legend = "right", common.legend = TRUE)
  
  
  
  

###############################################################################
###  How many Gene Modules across all the cells?
###############################################################################




spearman_dist_gene <- 1 - abs(cor(t(data_full_variable), 
                                  method = "spearman"))

spearman_dist_spectral_gene <- 
  dist(spectralEmbed(spearman_dist_gene, ndim = SPECTRAL_DIMS))

umap_embed_gene <- 
  uwot::umap(as.dist(spearman_dist_spectral_gene))

umap_embed_gene %>% 
  data.frame() %>% 
  ggplot(aes(x=X1, y=X2)) + 
  geom_point(col="blue") + 
  labs(x="UMAP1", y="UMAP2") + 
  theme_bw()




###############################################################################
###  Attempt GeneFish on the Data
###############################################################################

data_nonzero <- 
  data_full[rowSums(data_full) > 50, ]



probe_res_gene <-
  GeneFishing::probeFishability(as.matrix(data_nonzero),
                                potential_bait = immune_response_genes[immune_response_genes %in% 
                                                                         rownames(data_nonzero)],
                                n_rounds = 50,
                                ncores = 5)


doParallel::stopImplicitCluster()


gf_res_gene <-
  GeneFishing::geneFishing(as.matrix(data_nonzero),
                           bait_genes = probe_res_gene$best_bait,
                           fishing_rounds = 50,
                           n_probing_rounds = 50)

doParallel::stopImplicitCluster()

save(probe_res_gene, gf_res_gene, 
     file = here("output/tab_muris_sc/gf_res" %p% "immune_genes.RData"))


load(here("output/tab_muris_sc/gf_res" %p% "immune_genes.RData"))




















