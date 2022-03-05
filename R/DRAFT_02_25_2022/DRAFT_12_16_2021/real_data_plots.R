library(here)
source(here("R/libraries.R"))

core_names <- 
  c("pca_3", "pca_10", "pca_25",
    "euclid", "pear", "spearman",
    "side_ref_g50", "side_ref_g150", "side_ref_g300")

plot_names <- 
  c("PCA (3)", "PCA (10)", "PCA (25)", 
    "Euclidean", "Pearson", "Spearman", 
    "SIDE REF (50)", "SIDE REF (150)", "SIDE REF (300)")

names_map <- 
  data.frame(core_name = core_names,
             plot_name = plot_names)


loadAvgFolds <- function(files_to_load) {

  N = length(files_to_load)
  
  store_res <- NULL
  for(f in files_to_load) {
    print(f)
    load(f)
    results_list_sub <- 
      results_list[c("hclust_ari_celltypes",
                     "hclust_ari_celllevel")]
    if(is.null(store_res)){
      store_res <- results_list_sub
    } else{
      for(i in seq_len(length(store_res[[1]]))) {
        store_res[[1]][[i]] <- 
          store_res[[1]][[i]] + 
          results_list_sub[[1]][[i]]
        
        store_res[[2]][[i]] <- 
          store_res[[2]][[i]] + 
          results_list_sub[[2]][[i]]
      }
    }
  }
  ## average the results:
  for(i in seq_len(length(store_res[[1]]))) {
    store_res[[1]][[i]] <- store_res[[1]][[i]] / N
    store_res[[2]][[i]] <- store_res[[2]][[i]] / N
  }
  
  return(store_res)
}


cleanARIRes <- function(ari_list, use_wtd_pca=TRUE) {
  all_names <- gsub("_dist", "", names(ari_list))
  core_names <- unique(gsub("wtd", "", all_names))
  core_names <- unique(gsub("_spectral_d10|_spectral__d10", "", core_names))
  
  all_df <- data.frame(full_dist_name = all_names, 
                       name_wout_wtd = gsub("wtd", "", all_names),
                       ari = sapply(ari_list, function(x)x),
                       core_name = "")
  rownames(all_df) <- NULL
  
  all_df <- all_df[!grepl("pca_[[:digit:]]", all_df$full_dist_name),]
  
  all_df$full_dist_name <- gsub("pca_wtd", "pca", all_df$full_dist_name)
  
  for(n in core_names) {
    all_df$core_name[which(grepl(n, all_df$name_wout_wtd))] <- n
  }           
  
  all_df$y = ifelse(grepl("wtd", all_df$full_dist_name), 
                    "weighted",
                    ifelse(grepl("spectral", all_df$full_dist_name),
                           "spectral", "base"))
  
  all_df_2 <- tidyr::expand(all_df, core_name, y) %>% 
    left_join(all_df %>% select(core_name, ari, y)) %>% 
    left_join(names_map)
  
  return(all_df_2)
}

ARIHeatmap <- function(ari_list, title="", max_lim=0.45) {
  df = cleanARIRes(ari_list)
  p <-
    df %>% 
    mutate(plot_name = factor(plot_name, levels=rev(plot_names))) %>%
    ggplot(aes(y=plot_name, x=y, fill=ari, colour="")) + 
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "red",
                        limit = c(0,max_lim),
                        na.value = "black") + 
    scale_colour_manual(values=NA) +              
    guides(colour=guide_legend("Not Computed", 
                               override.aes=list(colour="dark grey"))) +
    labs(x = "Distance Measure", y = "Sub-Cat", 
         fill = "Adj. Rand Index") + 
    ggtitle(title) + 
    theme_bw()
  
  return(p)
}

ARIBarChart <- function(ari_list, title = "") {
  df = cleanARIRes(ari_list)
  p <-
    df %>% 
    mutate(plot_name = factor(plot_name, levels=rev(plot_names))) %>%
    ggplot(aes(y = plot_name,fill=y, x=ari)) + 
    geom_bar(stat="identity", 
             position = "dodge") + 
    labs(x = "Adj. Rand Index",
         fill = "Sub-Cat", 
         y = "Distance Measure") +
    xlim(0, 1) +
    ggtitle(title) +
    theme_bw()
    
  return(p)
}


##################################################
### Yielding the data:
rel_path <- here("output/tab_muris_sc/hclust_res/")

res_name_format <- "no_low_level_hclust_results_list_"

data_source1 <- "_droplet_and_facs"
data_source2 <- "_facs_only"
data_source3 <- "_droplet_only"

FOLDS=100

files_to_load <- 
  list.files(rel_path, pattern=res_name_format %p% 
               data_source1 %p% "_fold")
files_to_load <- 
  files_to_load[which(grepl(FOLDS%p%".RData",
                            files_to_load))]


res_both <- 
  loadAvgFolds(files_to_load)

### TODO: adjust the file listing here.
# res_facs <- 
#   loadAvgFolds(rel_path, res_name_format, 
#                data_source2, FOLDS = 20)
# 
# res_droplet <- 
#   loadAvgFolds(rel_path, res_name_format, 
#                data_source3, FOLDS = 25)

## Both source
p1 <- ARIHeatmap(res_both$hclust_ari_celltypes,
                 title = "Clustering from Low-Level Cell Types\n" %p%
                   "Both scRNA seq Sources")

p1

p1bar <- ARIBarChart(res_both$hclust_ari_celltypes,
                 title = "Clustering from Low-Level Cell Types\n" %p%
                   "Both scRNA seq Sources")


p2 <- ARIHeatmap(res_both$hclust_ari_celllevel,
                 title = "Clustering from Single Cell Level\n" %p%
                   "Both scRNA seq Sources")

p2

p2bar <- ARIBarChart(res_both$hclust_ari_celllevel,
                 title = "Clustering from Single Cell Level\n" %p%
                   "Both scRNA seq Sources")


## Facs only

p3 <- ARIHeatmap(res_facs$hclust_ari_celltypes,
                 title = "Clustering from Low-Level Cell Types\n" %p%
                   "FACS scRNA seq Source",
                 max_lim=0.5)

p3

p3bar <- ARIBarChart(res_facs$hclust_ari_celltypes,
                 title = "Clustering from Low-Level Cell Types\n" %p%
                   "FACS scRNA seq Source")


p4 <- ARIHeatmap(res_facs$hclust_ari_celllevel,
                 title = "Clustering from Single Cell Level\n" %p%
                   "FACS scRNA seq Source",
                 max_lim=0.5)

p4

p4bar <- ARIBarChart(res_facs$hclust_ari_celllevel,
                 title = "Clustering from Single Cell Level\n" %p%
                   "FACS scRNA seq Source")

## Droplet Only

p5 <- ARIHeatmap(res_droplet$hclust_ari_celltypes,
                 title = "Clustering from Low-Level Cell Types\n" %p%
                   "Droplet scRNA seq Source",
                 max_lim=0.5)

p5

p5bar <- ARIBarChart(res_droplet$hclust_ari_celltypes,
                 title = "Clustering from Low-Level Cell Types\n" %p%
                   "Droplet scRNA seq Source")


p6 <- ARIHeatmap(res_droplet$hclust_ari_celllevel,
                 title = "Clustering from Single Cell Level\n" %p%
                   "Droplet scRNA seq Source",
                 max_lim=0.5)

p6bar <- ARIBarChart(res_droplet$hclust_ari_celllevel,
                 title = "Clustering from Single Cell Level\n" %p%
                   "Droplet scRNA seq Source")

### Hierarchical Clustering Examination Plot:

load(here("output/tab_muris_sc/hclust_res/no_low_level_hclust_results_list" %p% 
            "__droplet_and_facs_fold1of100.RData"),
     verbose=TRUE)

tm_hclust_res_df <- 
  results_list$hclust_res_celltypes$side_ref_g50_dist$hclust_res$clust_res_df

load(here("output/tab_muris_sc/dist_res/full_dist_res__droplet_and_facs_fold1.RData"))
  
umap_res <- full_dist_res$umap_list$pca_25_umap[1:1008,1:2] %>% data.frame()



p <- PlotHClustResNew(tm_hclust_res_df, umap_res)

sum(tm_hclust_res_df$cell_group_high==1)

hclust_res_df <- 
  hclust_res_df %>% 
  filter(cell_group_high <= 17 & cluster_high <= 17)

ClusterPropHeatmap <- function(hclust_res_df) {
  tiss <- length(unique(hclust_res_df$cell_group_high))
  mat = matrix(nrow=tiss, ncol=tiss)
  for(i in seq_len(length(unique(hclust_res_df$cell_group_high)))) {
    for(j in seq_len(length(unique(hclust_res_df$cluster_high)))) {
      mat[i,j] = sum(hclust_res_df$cell_group_high==i & 
                       hclust_res_df$cluster_high==j) /
        sum(hclust_res_df$cell_group_high==i)
    }
  }
  as_dat = data.frame(mat)
  names(as_dat) = tissues
  as_dat$tissues = factor(as.character(1:17), 
                          levels=as.character(1:17))
  
  as_dat2 <- gather(as_dat, tissue, overlap, Bladder:Large_Intestine)
  
  as_dat2 %>% 
    ggplot(aes(x=tissues,y=tissue,fill=overlap)) + 
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "red",
                        na.value = "black") + 
    theme_bw() + 
    labs(fill = "% Overlap", x="Cluster",y="Ground Truth Tissue")
}
heatmap(
  as.matrix(mat), Rowv=NA,
  Colv=as.dendrogram(hclust(dist(t(as.matrix(mat))))),
  labRow = tissues,
  keep.dendro = FALSE
)


PlotHClustResNew <- function(hclust_res_df,
                             umap_res,
                             size = 0.5,
                             n_neighbors = 15,
                             min_dist = 0.01) {
  
  ## get HClust results
  
  ## Get Umap Embedding
  
  
  if(length(unique(hclust_res_df$cell_group_low)) == nrow(hclust_res_df)) {
    p <- 
      ggarrange( 
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cell_group_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("Actual Subgroups") +
          theme(title = element_text(size=8)) + 
          theme(legend.position = "none"),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cluster_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("HClust Subgroups") +
          theme(title = element_text(size=8)) + 
          theme(legend.position = "none"),
        ncol = 2
      )
  } else{
    p <- 
      ggarrange( 
        # umap_res %>% 
        #   mutate(celltype = factor(hclust_res_df$cell_group_low)) %>%
        #   ggplot(aes(x = X1, y = X2, col = celltype)) + 
        #   geom_point(size = size) + 
        #   theme_bw() +
        #   labs(x="UMAP1", y="UMAP2") +
        #   theme(legend.position = "none") +
        #   ggtitle("Actual Low-Level Groups") +
        #   theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cell_group_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("Actual High-Level Groups") +
          theme(legend.position = "none") + 
          theme(title = element_text(size=8)),
        umap_res %>% 
          mutate(celltype = factor(hclust_res_df$cluster_high)) %>%
          ggplot(aes(x = X1, y = X2, col = celltype)) + 
          geom_point(size = size) + 
          theme_bw() + 
          labs(x="UMAP1", y="UMAP2") +
          ggtitle("HClust Subgroups") +
          theme(legend.position = "none") +
          theme(title = element_text(size=8)),
        ncol = 2
      )
  }
  
  return(p)
  
}
  

### Visualize the dendrogram of cell groups



list.files(here("output/tab_muris_sc/hclust_res"))


load(here("output/tab_muris_sc/hclust_res/" %p% 
            "no_low_level_hclust_results_list_" %p% 
            "_droplet_only_fold1of500.RData"),
     verbose=TRUE)

load(here("output/tab_muris_sc/hclust_res/" %p% 
            "no_low_level_hclust_results_list_" %p% 
            "_droplet_and_facs_fold1of50.RData"),
     verbose=TRUE)

load(here("output/tab_muris_sc/hclust_res/" %p% 
            "no_low_level_hclust_results_list_" %p% 
            "_facs_only_fold1of20.RData"),
     verbose=TRUE)



library(dendextend)
library(circlize)
library(ggdendro)


## ARI bar chart samples:
ARIBarChart(results_list$hclust_ari_celltypes)
ARIBarChart(results_list$hclust_ari_celllevel)

hclust_res <- 
  results_list$hclust_res_celltypes$side_ref_g150_dist$hclust_res$hclust_res

meta_data <- 
  results_list$meta_data

clust_res_df <-
  results_list$hclust_res_celltypes[[1]]$hclust_res$clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)

dendro_df <- 
  dendro_data(hclust_res)

n_celltypes <- 
  length(unique(
    results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_low %p% 
      results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_high
  ))

n_tissues <- 
  length(unique(
    results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_high))

### Creating dendrogram from cell type to "ground truth" tissue

## TODO: filler distance matrix for now.:
dist_mat <- matrix(0.1, nrow(clust_res_sample), nrow(clust_res_sample))

diag(dist_mat) = 0

tiss_hclust <- 
  HClust_with_start_clusters(
    clusters = as.numeric(factor(clust_res_sample$tissue,
                                 levels = unique(clust_res_sample$tissue))), 
    subg_inds = rep(1, nrow(clust_res_sample)),
    dissim_mat = dist_mat, 
    cell_ids = clust_res_sample$cell_index)


tiss_hclust_res <- tiss_hclust$hclust_res
{
  dend_tiss <- 
    tiss_hclust_res %>% 
    as.dendrogram()
  
  cell_type_names <-
    str_to_title(clust_res_sample$cell_type[labels(dend_tiss)])
  
  tissues <- 
    str_to_title(clust_res_sample$tissue[labels(dend_tiss)])
  
  tissues <- 
    case_when(
      tissues == "Brain_non-Myeloid" ~ "Brain, Non-Myeloid",
      tissues == "Brain_myeloid" ~ "Brain, Myeloid",
      tissues == "Limb_muscle" ~ "Limb Muscle", 
      tissues == "Large_intestine" ~ "Large Intestine",
      tissues == "Mammary_gland" ~ "Mammary Gland",
      TRUE ~ tissues
      
    )
  
  cell_type_names <- 
    gsub("Of", "of", cell_type_names)
  
  cell_type_names <- 
    case_when(
      grepl("luminal epithelial cell of mammary gland",
            cell_type_names, 
            ignore.case = TRUE) ~ "Luminal Epithelial Cell\n" %p%
        "       of Mammary Gland",
      grepl("endothelial cell of hepatic sinusoid",
            cell_type_names, 
            ignore.case = TRUE) ~ "Endothelial Cell \n" %p%
        "       of Hepatic Sinusoid",
      grepl("professional antigen presenting cell",
            cell_type_names, 
            ignore.case = TRUE) ~ "Professional Antigen\n" %p%
        "       Presenting Cell",
      grepl("Slamf1-Negative Multipotent Progenitor Cell",
            cell_type_names, 
            ignore.case = TRUE) ~ "Slamf1-Negative\n" %p%
        "       Multipotent\n" %p% 
        "       Progenitor Cell",
      TRUE ~ cell_type_names
    )
  
  labels(dend_tiss) <- cell_type_names
  dend_tiss <- 
    dend_tiss %>% 
    set("branches_k_color", 
        k = n_tissues) %>%
    set("labels_cex", 0.4) #%>%
  #color_labels()
  
  labels(dend_tiss) <- 
    as.numeric(factor(get_leaves_branches_col(dend_tiss),
                      levels = unique(
                        get_leaves_branches_col(dend_tiss)
                        ))) %p% 
    ". " %p% 
    labels(dend_tiss)
  
  #plot(dend, horiz=TRUE)
  
  
  ggd1 <- as.ggdend(dend_tiss)
  
  ggplot(ggd1, horiz=TRUE) + 
    theme(axis.text = element_text(size=2))
  
  circlize_dendrogram(dend_tiss, 
                      labels_track_height = 0.4, 
                      dend_track_height = .3,
                      labels=TRUE)
  legend("topleft", 
         legend = seq_len(n_tissues) %p% ". " %p% 
           unique(tissues), 
         fill = unique(get_leaves_branches_col(dend_tiss)),
         cex=0.5, #pch=1, 
         pt.cex = 1)
  
}



## dendro plotting function:
prod_circlized_dendrogram <- 
  function(hclust_res, clust_res_df, sample_clust_res_df = NULL) {
    dend <- 
      hclust_res %>% 
      as.dendrogram()
    
    if(!is.null(sample_clust_res_df)) {
      dend <- 
        prune(dend, setdiff(labels(dend), 
                            sample_clust_res_df$cell_index))
    }
    # labels(dend) <- gsub("[^B ][^T ] cell", "", 
    #                      clust_res_df$cell_type[labels(dend)])
    
    ## convert to numeric
    #labels(dend) <- clust_res_df$[labels(dend), ]
    rownames(clust_res_df) <- clust_res_df$cell_index
    clust_res_df[labels(dend),]$cell_type
    
    
    
    cell_type_names <-
      str_to_title(clust_res_df[labels(dend),]$cell_type)
    
    cell_type_names <- 
      gsub("Of", "of", cell_type_names)
    
    cell_type_names <- 
      case_when(
        grepl("luminal epithelial cell of mammary gland",
              cell_type_names, 
              ignore.case = TRUE) ~ "Luminal Epithelial Cell\n" %p%
          "       of Mammary Gland",
        grepl("endothelial cell of hepatic sinusoid",
              cell_type_names, 
              ignore.case = TRUE) ~ "Endothelial Cell \n" %p%
          "       of Hepatic Sinusoid",
        grepl("professional antigen presenting cell",
              cell_type_names, 
              ignore.case = TRUE) ~ "Professional Antigen\n" %p%
          "       Presenting Cell",
        grepl("Slamf1-Negative Multipotent Progenitor Cell",
              cell_type_names, 
              ignore.case = TRUE) ~ "Slamf1-Negative\n" %p%
          "       Multipotent\n" %p% 
          "       Progenitor Cell",
        TRUE ~ cell_type_names
      )
    
    labels(dend) <- cell_type_names
    dend <- 
      dend %>% 
      set("branches_k_color", 
          k = n_tissues) %>%
      set("labels_cex", 0.4) #%>%
    #color_labels()
    
    labels(dend) <- 
      as.numeric(as.factor(get_leaves_branches_col(dend))) %p% 
      ". " %p% 
      labels(dend)
    
    #plot(dend, horiz=TRUE)
    
    circle_dend_p <- 
      circlize_dendrogram(dend, 
                          labels_track_height = 0.4, 
                          dend_track_height = .3,
                          labels=TRUE)
    
    return(list(dend=dend,
                circle_dend_p = circle_dend_p))
    
}

side_ref_g150_circlized <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_celltypes$
      side_ref_g150_dist$hclust_res$hclust_res, 
    clust_res_df, sample_clust_res_df = clust_res_sample)

pca10_circlized <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_celltypes$
      pca_10_dist$hclust_res$hclust_res, 
    clust_res_df, sample_clust_res_df = clust_res_sample)


tm_tanglegram <- function(dend1, dend2, 
                          dend1name = "", 
                          dend2name = "") {
  ## remove numbers from labels
  labels(dend1) <- 
    gsub("^\\d{1,}\\. ", "", labels(dend1))
  labels(dend1) <- 
    gsub(" {2,}", "", labels(dend1))
  
  labels(dend2) <- 
    gsub("^\\d{1,}\\. ", "", labels(dend2))
  labels(dend2) <- 
    gsub(" {2,}", "", labels(dend2))
  labels(dend2)
  
  dl <- dendlist(dend1, 
                 dend2)
  tanglegram(dl, sort = TRUE,
             common_subtrees_color_lines = FALSE, 
             highlight_distinct_edges  = TRUE, 
             main_left = dend1name, 
             main_right = dend2name,
             highlight_branches_lwd = FALSE,
             common_subtrees_color_branches  = TRUE,
             main = "Entanglement: " %p% round(entanglement(dend1, dend2), 2))
  
  
  
}


tm_tanglegram(dend1 = side_ref_g150_circlized$dend,
              dend2 = pca10_circlized$dend,
              dend1name = "SIDEREF (150)",
              dend2name = "PCA (10)")



###############################################################################
### Code for Splatter Results
###############################################################################
list.files(here("output/splatter_sim/hclust_res"))

SPLAT_FILE_NAME <- "splatter_sim_skin_eq_grps_lowvar_rep1_hclustres.RData"

load(here("output/splatter_sim/hclust_res/" %p% SPLAT_FILE_NAME),
     verbose = TRUE)


ARIBarChart(results_list$hclust_ari_celltypes)
ARIBarChart(results_list$hclust_ari_celllevel)


### TODO: function to compute within cell subgroup ARI.

splatter_subgroup_map <-
  data.frame(group_num = levels(results_list$meta_data$SubGroup)) %>% 
  mutate(cell_group_high = row_number()) %>% 
  mutate(subgroup_label = c(
    "[1] High DE Prob",
    "[2] Low DE Prob",
    "[3,4,5] High Ind. DE, Low Shared DE",
    "[6,7,8] High Ind. DE, High Shared DE",
    "[9,10,11] Low Ind DE, High Shared DE",
    "[12, 13, 14] Low Ind DE, Low Shared DE",
    "[15, 16, 17] Low Ind DE, Low Shared DE, High DE Variance",
    "[18, 19, 20] High Ind DE, Low Shared DE, High DE Variance"
  ))

splatter_new_group_map <- 
  data.frame(
    cell_group_low = 1:20,
    new_cell_group_low = c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                           2, 
                           12, 13, 14, 15, 16, 17, 18, 19, 20)
  )


splatterDendrogram <- 
  function(hclust_res, 
           clust_res_df,
           splatter_subgroup_map,
           y_adj = 10) {
    
    clust_res_df <- 
      clust_res_df %>% 
      left_join(splatter_subgroup_map) %>% 
      left_join(splatter_new_group_map)
    
    sample_clust_res_df <-
      clust_res_df %>% 
      mutate(cell_num = seq_len(nrow(clust_res_df))) %>%  
      group_by(cell_group_low) %>% 
      sample_n(1)
    
    dend <- 
      hclust_res %>% 
      as.dendrogram()
    
    dend <- 
      dend %>% 
      color_branches(
        col=colorspace::rainbow_hcl(n=nrow(splatter_subgroup_map),c=90,l=50)[
          unique(clust_res_df$cell_group_high[labels(dend)])],
        clusters=clust_res_df$cell_group_high[labels(dend)]) %>%
      set("labels_cex", 0.4) 
    
    if(!is.null(sample_clust_res_df)) {
      dend <- 
        prune(dend, setdiff(seq_len(nrow(clust_res_df)), 
                            sample_clust_res_df$cell_num))
    }
    
    plot_subgroups <- 
      clust_res_df$subgroup_label[labels(dend)]
    
    cell_type_names <- 
      "Group " %p% 
      clust_res_df$new_cell_group_low[labels(dend)]
    
    
    
    labels(dend) <- cell_type_names
    plot(dend, 
         xlim=c(0,length(cell_type_names)+3),
         ylim=c(0,get_nodes_attr(dend, "height")[1]+y_adj))
    
    dend %>% rect.dendrogram(k=nrow(splatter_subgroup_map),border="grey")
    #abline(h = heights_per_k.dendrogram(dend)[as.character(nrow(splatter_subgroup_map))],
    #       lty=2)
    
    legend(x="topright", 
           legend = splatter_subgroup_map$subgroup_label, 
           fill = colorspace::rainbow_hcl(n=nrow(splatter_subgroup_map),c=90,l=50),
           cex=0.65, #pch=1, 
           pt.cex = 0.7)
  }


### 2000 Variable Genes
load(here("output/splatter_sim/hclust_res/" %p% SPLAT_FILE_NAME),
     verbose = TRUE)

splatterDendrogram(
  results_list$hclust_res_celltypes$pca_wtd10_dist$
    hclust_res$hclust_res, 
  results_list$hclust_res_celltypes$pca_wtd10_dist$
    hclust_res$clust_res_df,
  splatter_subgroup_map)


splatterDendrogram(
  results_list$hclust_res_celltypes$pca_wtd10_dist_spectral_wtd_d10_dist$
    hclust_res$hclust_res, 
  results_list$hclust_res_celltypes$pca_wtd10_dist_spectral_wtd_d10_dist$
    hclust_res$clust_res_df,
  splatter_subgroup_map)


### SIDEREF
splatterDendrogram(
  results_list$hclust_res_celltypes$side_ref_g150_dist$
    hclust_res$hclust_res, 
  results_list$hclust_res_celltypes$side_ref_g150_dist$
    hclust_res$clust_res_df,
  splatter_subgroup_map, y_adj = 150)


splatterDendrogram(
  results_list$hclust_res_celltypes$side_ref_g300_dist_spectral_wtd_d10_dist$
    hclust_res$hclust_res, 
  results_list$hclust_res_celltypes$side_ref_g300_dist_spectral_wtd_d10_dist$
    hclust_res$clust_res_df,
  splatter_subgroup_map)







### Load average results 
splat_hclust_out_path <- 
  here("output/splatter_sim/hclust_res/")

## Get filenames:
splat_hclust_out_files <- 
  list.files(splat_hclust_out_path)
    
splat_hclust_lowvar_2000g_files <- 
  splat_hclust_out_files[
    str_detect(splat_hclust_out_files, "g\\.RData", negate = TRUE) 
  ]

splat_hclust_lowvar_1000g_files <- 
  splat_hclust_out_files[
    str_detect(splat_hclust_out_files, "1000g") 
  ]

splat_hclust_lowvar_allg_files <- 
  splat_hclust_out_files[
    str_detect(splat_hclust_out_files, "allg") 
  ]


splat_hclust_highvar_2000g_files <- 
  splat_hclust_out_files[
    str_detect(splat_hclust_out_files, "highvar.*2000g") 
  ]

### Load average folds:
avg_splat_res_low_var_2000g <- 
  loadAvgFolds(splat_hclust_out_path %p% 
                 splat_hclust_lowvar_2000g_files) 


avg_splat_res_low_var_1000g <- 
  loadAvgFolds(splat_hclust_out_path %p% 
                 splat_hclust_lowvar_1000g_files)

avg_splat_res_low_var_allg <- 
  loadAvgFolds(splat_hclust_out_path %p% 
                 splat_hclust_lowvar_allg_files)


avg_splat_res_high_var_2000g <- 
  loadAvgFolds(splat_hclust_out_path %p% 
                 splat_hclust_highvar_2000g_files)


### Plot results (they're all the same.):
ARIBarChart(avg_splat_res_low_var_2000g$hclust_ari_celltypes,
            title = "From Cell Types, Standard Noise")

ARIBarChart(avg_splat_res_low_var_2000g$hclust_ari_celllevel,
            title = "From Cell Level, Standard Noise")


### Umap plot of splatter -----------------------------------------------------

load(splat_hclust_out_path %p% splat_hclust_lowvar_2000g_files[1],
     verbose=TRUE)

results_list$dist_res$plot_list$pca_25_plot_list

results_list$dist_res$plot_list$pca_wtd10_plot_list

pear_wtd_umap_df <- results_list$dist_res$umap_list$pear_dist_spectral_wtd_d10_dist %>% 
  data.frame() 

results_list$dist_res$plot_list$euclid_plot_list

pear_wtd_umap_df %>% 
  mutate(celltype = results_list$dist_res$plot_list[[1]]$data$celltype) %>% 
  ggplot(aes(x=UMAP1,y=UMAP2,col=celltype))+geom_point() + 
  theme_bw() + ggtitle("Pearson Spectral (10 Dims)") + 
  labs(col="Cell Type")

spear_wtd_umap_df <- results_list$dist_res$umap_list$spearman_dist_spectral_wtd_d10_dist %>% 
  data.frame()

names(spear_wtd_umap_df) <- c("UMAP1", "UMAP2")
spear_wtd_umap_df %>%
  mutate(celltype = results_list$dist_res$plot_list[[1]]$data$celltype) %>% 
  ggplot(aes(x=UMAP1,y=UMAP2,col=celltype))+geom_point() + 
  theme_bw() + ggtitle("Spearman Spectral (10 Dims)") + 
  labs(col="Cell Type")

pca_spectral_umap_df <- results_list$dist_res$umap_list$pca_wtd25_dist_spectral_wtd_d10_dist %>% 
  data.frame()

names(pca_spectral_umap_df) <- c("UMAP1", "UMAP2")
pca_spectral_umap_df %>%
  mutate(celltype = results_list$dist_res$plot_list[[1]]$data$celltype) %>% 
  ggplot(aes(x=UMAP1,y=UMAP2,col=celltype))+geom_point() + 
  theme_bw() + ggtitle("Euclidean Applied to 25Wtd. PCs\nSpectral (10 Dims)") + 
  labs(col="Cell Type")



### Umap plot of splatter -----------------------------------------------------

# ARIBarChart(avg_splat_res_low_var_1000g$hclust_ari_celltypes,
#             title = "1000 Genes, From Cell Types")
# 
# ARIBarChart(avg_splat_res_low_var_1000g$hclust_ari_celllevel,
#             title = "1000 Genes, From Cell Level")


ARIBarChart(avg_splat_res_high_var_2000g$hclust_ari_celltypes,
            title = "From Cell Types, High Noise")

ARIBarChart(avg_splat_res_high_var_2000g$hclust_ari_celllevel,
            title = "From Cell Level, High Noise")



### TM:  ----------------------------------------------------------------------


tm_hclust_out_path <- 
  here("output/tab_muris_sc/hclust_res/")

## Get filenames:
tm_hclust_out_files <- 
  list.files(tm_hclust_out_path)

tm_hclust_droplet_files <- 
  tm_hclust_out_files[
    str_detect(tm_hclust_out_files, "droplet") 
  ]

tm_hclust_facs_files <- 
  tm_hclust_out_files[
    str_detect(tm_hclust_out_files, "facs") 
  ]

avg_tm_droplet_res <- 
  loadAvgFolds(tm_hclust_out_path %p% 
                 tm_hclust_droplet_files)

avg_tm_facs_res <- 
  loadAvgFolds(tm_hclust_out_path %p% 
                 tm_hclust_facs_files)



ARIBarChart(avg_tm_droplet_res$hclust_ari_celltypes,
            title = "From Cell Types, Droplet Datasets")

ARIBarChart(avg_tm_droplet_res$hclust_ari_celllevel,
            title = "From Cell Level, Droplet Datasets")

ARIBarChart(avg_tm_facs_res$hclust_ari_celltypes,
            title = "From Cell Types, FACS Datasets")

ARIBarChart(avg_tm_facs_res$hclust_ari_celllevel,
            title = "From Cell Level, FACS Datasets")



### TABULA MURIS DENDROGRAMS --------------------------------------------------

load(tm_hclust_out_path %p% tm_hclust_facs_files[1], 
     verbose=TRUE)

meta_data <- results_list$meta_data

n_celltypes <- 
  length(unique(
    results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_low %p% 
      results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_high
  ))

n_tissues <- 
  length(unique(
    results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_high))

clust_res_df <-
  results_list$hclust_res_celltypes$pear_dist$
  hclust_res$clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)


pear_circlized_facs <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_celltypes$pear_dist$
      hclust_res$hclust_res, 
    clust_res_df, 
    sample_clust_res_df = clust_res_sample)



############ Pear spectral #######################


clust_res_df <-
  results_list$hclust_res_celltypes$pear_dist_spectral_wtd_d10_dist$
  hclust_res$clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)


pear_spect_circlized_facs <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_celltypes$pear_dist_spectral_wtd_d10_dist$
      hclust_res$hclust_res, 
    clust_res_df, 
    sample_clust_res_df = clust_res_sample)


############ Tissue results #######################

clust_res_df <-
  results_list$hclust_res_tissue$pear_dist_spectral_wtd_d10_dist$
  clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)


tissue_spect_circlized_facs <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_tissue$pear_dist_spectral_wtd_d10_dist$
      hclust_res, 
    clust_res_df, 
    sample_clust_res_df = clust_res_sample)



# tm_tanglegram(dend1 = pear_circlized_facs$dend,
#               dend2 = tissue_spect_circlized_facs$dend,
#               dend1name = "Pearson (FACS)",
#               dend2name = "Tissue Labels (FACS)")
# 
# 
# tm_tanglegram(dend1 = pear_circlized_facs$dend,
#               dend2 = tissue_spect_circlized_facs$dend,
#               dend1name = "Pearson (FACS)",
#               dend2name = "Tissue Labels (FACS)")


### TABULA MURIS DENDROGRAMS (DROPLET) ---------------------------------------

load(tm_hclust_out_path %p% tm_hclust_droplet_files[1], 
     verbose=TRUE)

meta_data <- results_list$meta_data

n_celltypes <- 
  length(unique(
    results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_low %p% 
      results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_high
  ))

n_tissues <- 
  length(unique(
    results_list$hclust_res_celltypes[[1]]$
      hclust_res$
      clust_res_df$
      cell_group_high))

clust_res_df <-
  results_list$hclust_res_celltypes$pear_dist$
  hclust_res$clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)


pear_circlized_droplet <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_celltypes$pear_dist$
      hclust_res$hclust_res, 
    clust_res_df, 
    sample_clust_res_df = clust_res_sample)



############ Pear spectral #######################


clust_res_df <-
  results_list$hclust_res_celltypes$pear_dist_spectral_wtd_d10_dist$
  hclust_res$clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)


pear_spect_circlized_droplet <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_celltypes$pear_dist_spectral_wtd_d10_dist$
      hclust_res$hclust_res, 
    clust_res_df, 
    sample_clust_res_df = clust_res_sample)


############ Tissue results #######################

clust_res_df <-
  results_list$hclust_res_tissue$pear_dist_spectral_wtd_d10_dist$
  clust_res_df %>% 
  left_join(meta_data %>% 
              select(cell_index = cell_id, tissue, cell_type))

clust_res_sample <- 
  clust_res_df %>% 
  mutate(cell_num = seq_len(nrow(clust_res_df))) %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::sample_n(1)


tissue_spect_circlized_droplet <-
  prod_circlized_dendrogram(
    hclust_res = results_list$hclust_res_tissue$pear_dist_spectral_wtd_d10_dist$
      hclust_res, 
    clust_res_df, 
    sample_clust_res_df = clust_res_sample)




tm_tanglegram(dend1 = pear_circlized_droplet$dend,
              dend2 = tissue_spect_circlized_facs$dend,
              dend1name = "Pearson (Droplet)",
              dend2name = "Tissue Labels (Droplet)")


tm_tanglegram(dend1 = pear_circlized_droplet$dend,
              dend2 = tissue_spect_circlized_facs$dend,
              dend1name = "Pearson (FACS)",
              dend2name = "Tissue Labels (Droplet)")


################################################################################
### Group-wise Heatmaps
################################################################################



### SPLATTER GROUPWISEDISTANCEHEATMAP

source(here("R/subgroup_composition_gene_contrib.R"))

df_group_list <- 
  data.frame(cell_group_low = 
               as.numeric(gsub("Group", "", 
                               as.character(results_list$meta_data$Group)))) %>% 
  left_join(splatter_new_group_map)
  
  

groupwiseDistanceHeatmap(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                         results_list$dist_res$dist_list$side_ref_g150_dist,
                         title = "SIDEREF (150)",
                         hclust = FALSE,
                         preset_levels = as.character(1:20))



groupwiseDistanceHeatmap(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                         results_list$dist_res$dist_list$side_ref_g150_dist_spectral_wtd_d10_dist,
                         title = "SIDEREF (150), Spectral (10 Dim)",
                         hclust = FALSE,
                         preset_levels = as.character(1:20))



groupwiseDistanceHeatmap(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                         results_list$dist_res$dist_list$pca_wtd25_dist,
                         title = "PCA (25)",
                         hclust = TRUE,
                         preset_levels = as.character(1:20))



groupwiseDistanceHeatmap(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                         results_list$dist_res$dist_list$pca_wtd10_dist_spectral_wtd_d10_dist,
                         title = "PCA (25), Spectral (10 Dim)",
                         hclust = TRUE,
                         preset_levels = as.character(1:20))



### Legend for reference of subgroups
plot.new()
legend("topright", 
       legend = splatter_subgroup_map$subgroup_label,
       cex=0.5, #pch=1, 
       pt.cex = 1)



###############################################################################
### Dendrogram from the relative group distances?
###############################################################################


relativeDistDendro <- 
  function(group_labels, 
           dist_mat,
           group_subgroup_map,
           hclust_method = "ward.D2",
           y_adj = 2) {
    
    group_dist_mat <- getDistDf(group_labels, dist_mat, gather = FALSE) %>% 
      dplyr::select(-group)
    
    ## symmetrize
    group_dist_mat = 0.5*(group_dist_mat+t(group_dist_mat))
    diag(group_dist_mat) = 0
    
    hclust_res <- stats::hclust(d = as.dist(group_dist_mat), 
                                method = hclust_method)
    
    ### create dendrogram plot
    dend <- 
      hclust_res %>% 
      as.dendrogram()
    
    labels(dend) <- as.numeric(labels(dend))
    
    group_subgroup_map <- group_subgroup_map %>% 
      dplyr::arrange(new_cell_group_low)
    
    dend <- 
      dend %>% 
      color_branches(
        col = colorspace::rainbow_hcl(n=nrow(splatter_subgroup_map),c=90,l=50)[
          unique(group_subgroup_map$cell_group_high[labels(dend)])],
        clusters = group_subgroup_map$cell_group_high[labels(dend)]) %>%
      set("labels_cex", 0.4) 
    
    
    plot_subgroups <- 
      group_subgroup_map$subgroup_label[labels(dend)]
    
    cell_type_names <- 
      "Group " %p% 
      group_subgroup_map$new_cell_group_low[labels(dend)]
    
    
    
    labels(dend) <- cell_type_names
    plot(dend, 
         xlim=c(0,length(cell_type_names)+3),
         ylim=c(0,get_nodes_attr(dend, "height")[1]+y_adj))
    
    dend %>% 
      rect.dendrogram(
        k = length(unique(group_subgroup_map$cell_group_high)), 
        border="grey")
    
    legend(x="topright", 
           legend = unique(group_subgroup_map$subgroup_label), 
           fill = colorspace::rainbow_hcl(n=nrow(splatter_subgroup_map),c=90,l=50),
           cex=0.65, #pch=1, 
           pt.cex = 0.7)
}



group_subgroup_map <- clust_res_df %>% 
  dplyr::distinct(cell_group_low, cell_group_high) %>% 
  left_join(splatter_new_group_map) %>% 
  left_join(subgroup_label_map)


relativeDistDendro(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                   results_list$dist_res$dist_list$pca_wtd10_dist_spectral_wtd_d10_dist,
                   group_subgroup_map,
                   hclust_method = "ward.D2",
                   y_adj = 1)


relativeDistDendro(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                   results_list$dist_res$dist_list$pca_10_dist,
                   group_subgroup_map,
                   hclust_method = "ward.D2",
                   y_adj = 1)

relativeDistDendro(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                   results_list$dist_res$dist_list$pca_wtd25_dist,
                   group_subgroup_map,
                   hclust_method = "ward.D2",
                   y_adj = 1)


relativeDistDendro(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                   results_list$dist_res$dist_list$side_ref_g150_dist,
                   group_subgroup_map,
                   hclust_method = "ward.D2",
                   y_adj = 1)


relativeDistDendro(gsub("Group", "", as.character(df_group_list$new_cell_group_low)), 
                   results_list$dist_res$dist_list$side_ref_g150_dist_spectral_d10_dist,
                   group_subgroup_map,
                   hclust_method = "ward.D2",
                   y_adj = 1)



###############################################################################
### Real Data Groupwise Heatmaps
###############################################################################


### load splatter res


# load(here("output/tab_muris_sc/hclust_res/no_low_level_hclust_results_list_" %p% 
#             "_droplet_only_fold1of20.RData"))


load(here("output/tab_muris_sc/dist_res/full_dist_res__droplet_only_fold1.RData"),
     verbose = TRUE)


### SIDEREF 150 (HClust Ordered)
cell_types <- unique(full_dist_res$plot_list$side_ref_g150_plot_list$
                       data$celltype) %>% 
  sort()

keep_groups <-
  as.logical(table(full_dist_res$plot_list$side_ref_g150_plot_list$
                     data$celltype) > 5)

keep_cell_types <- cell_types[keep_groups]

keep_cells <- which(full_dist_res$plot_list$side_ref_g150_plot_list$
                      data$celltype %in% keep_cell_types)

group_labels <- 
  full_dist_res$plot_list$side_ref_g150_plot_list$data$celltype[
    keep_cells
  ]

dist_mat <- full_dist_res$dist_list$side_ref_g150_dist[keep_cells,
                                                       keep_cells]


groupwiseDistanceHeatmap(group_labels, 
                         dist_mat,
                         title = "SIDEREF (150)",
                         hclust = TRUE,
                         numeric_x=TRUE) #+ 
  #theme(axis.text.x=element_text(angle=-30))
  
group_labels <-
  data.frame(cell_type = group_labels) %>% 
  left_join(meta_data_filt)

group_labels = group_labels$cell_type %p% ", " %p% 
  group_labels$tissue


groupwiseDistanceHeatmap(group_labels, 
                         dist_mat,
                         title = "SIDEREF (150)",
                         hclust = FALSE,
                         preset_levels = 
                           meta_data_filt$cell_type %p% 
                           ", " %p% meta_data_filt$tissue,
                         numeric_x=TRUE)


groupwiseDistanceHeatmap(group_labels[more_proximal_to_keep], 
                         dist_mat[more_proximal_to_keep, 
                                  more_proximal_to_keep],
                         title = "SIDEREF (150)",
                         hclust = FALSE,
                         preset_levels = 
                           meta_data_filt$cell_type %p% 
                           ", " %p% meta_data_filt$tissue,
                         numeric_x=TRUE)


### PCA 25 (HClust Ordered)
cell_types <- unique(full_dist_res$plot_list$pca_25_plot_list$
                       data$celltype) %>% 
  sort()

keep_groups <-
  as.logical(table(full_dist_res$plot_list$pca_25_plot_list$
                     data$celltype) > 5)

keep_cell_types <- cell_types[keep_groups]

keep_cells <- which(full_dist_res$plot_list$pca_25_plot_list$
                      data$celltype %in% keep_cell_types)

group_labels <- 
  full_dist_res$plot_list$pca_25_plot_list$data$celltype[
    keep_cells
  ]

dist_mat <- full_dist_res$dist_list$pca_25_dist[keep_cells,
                                                     keep_cells]


groupwiseDistanceHeatmap(group_labels, 
                         dist_mat,
                         title = "PCA (25)",
                         hclust = TRUE,
                         numeric_x=TRUE)


### ORDERED BY TISSUES INSTEAD
meta_data_filt <- results_list$meta_data %>% 
  dplyr::filter(cell_type %in% keep_cell_types) %>% 
  group_by(cell_type) %>% 
  sample_n(1) %>% 
  dplyr::distinct(tissue, cell_type) %>% 
  dplyr::arrange(tissue)


group_labels <-
  data.frame(cell_type = group_labels) %>% 
  left_join(meta_data_filt)

group_labels = group_labels$cell_type %p% ", " %p% 
  group_labels$tissue


more_proximal_to_keep <-
  which(!group_labels %in% 
          c("hepatocyte, Liver",
            "kidney proximal straight tubule epithelial cell, Kidney",
            "proerythroblast, Marrow"))

groupwiseDistanceHeatmap(group_labels, 
                         dist_mat,
                         title = "PCA (25)",
                         hclust = FALSE,
                         preset_levels = 
                           meta_data_filt$cell_type %p% 
                           ", " %p% meta_data_filt$tissue,
                         numeric_x=TRUE)



groupwiseDistanceHeatmap(group_labels[more_proximal_to_keep], 
                         dist_mat[more_proximal_to_keep, 
                                  more_proximal_to_keep],
                         title = "PCA (25)",
                         hclust = FALSE,
                         preset_levels = 
                           meta_data_filt$cell_type %p% 
                           ", " %p% meta_data_filt$tissue,
                         numeric_x=TRUE)





