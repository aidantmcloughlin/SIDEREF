################################################################################
### Title:        Gene Ontology cleaning, usage
### Author:       Aidan McLoughlin
################################################################################

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
## Statistics
library(Hmisc)
library(diffusr)


`%p%` <- function(x, y) paste0(x, y)

### load gene ontology data 

## TERM: central nervous system development
## ID: GO:0007417
## http://www.informatics.jax.org/go/term/GO:0007417
brain_go_term_summary <- 
  read_xlsx(path = here("data/gene_ontology/" %p% 
                          "GO_term_summary_20210218_185754.xlsx")) %>%
  clean_names()


## Proteoform: Different forms of a protein produced from variety of 
##    sequence variations, splice isoforms, post-translational modifications
table(brain_go_term_summary$proteoform, exclude = NULL)
## all NA
table(brain_go_term_summary$qualifier, exclude = NULL)
brain_go_term_summary <- 
  brain_go_term_summary %>% 
  dplyr::select(-qualifier)

## EVIDENCE ACRONYMS ----------------------------------------------------------

table(brain_go_term_summary$evidence, exclude = NULL)

## HOMOLOGY RELATED:
##    IBA: Inferred from biological aspect of ancestor (#4)
##    ISO: Inferred from sequence orthology (#3)
##    ISS: Inferred from sequence or structural similarity



## EXPERIMENTAL RELATED:
##    IDA: Inferred from direct assay (#5)
##    IGI: Inferred from genetic interaction (#2)
##    IMP: Inferred from mutant phenotype (#1)



## AUTOMATED:
##    IEA: Inferred from electronic annotation (automated)

## OTHER:
##    IC:  Inferred by curator
##    NAS: Non-Traceable author statement
##    TAS: Traceable author statement

brain_go_term_summary <- 
  brain_go_term_summary %>% 
  mutate(evidence_group = 
           case_when(
             evidence %in% c("EXP", "HMP", "HGI", "HDA", "HEP",
                             "IDA", "IEP", "IGI", "IMP", "IPI") ~ "Experimental",
             evidence %in% c("IAS", "IBA", "IBD", "IKR", "IMR", 
                             "IRD", "ISA", "ISM", "ISO", "ISS") ~ "Homology", 
             evidence %in% c("IEA", "RCA") ~ "Automated",
             evidence %in% c("IC", "NAS", "ND", "TAS") ~ "Other"
           ))



brain_go_term_summary %>% 
  group_by(annotated_term) %>% 
  summarize(count=n()) %>% 
  arrange(-count) %>% head(n = 10)


## For example extract cerebellum development genes

## many annotations for a gene?
many_annotations <- 
  brain_go_term_summary %>% 
  group_by(mgi_gene_marker_id) %>% 
  mutate(count = n()) %>% 
  dplyr::filter(count >= 25)



## number of gene marker IDs vs symbol vs name?
length(unique(brain_go_term_summary$mgi_gene_marker_id))
length(unique(brain_go_term_summary$symbol))
length(unique(brain_go_term_summary$name))
## Nice. 931 unique gene-id-symbol-name-pairings..



## Connect Gene names to annotation from multiome data ------------------------

load(here("data/share_seq/brain/brain.rda"))

p1 <- DimPlot(brain, 
              reduction = "umap.rna", 
              group.by = "seurat_clusters", label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("RNA")

p2 <- DimPlot(brain, 
              reduction = "umap.atac", 
              group.by = "seurat_clusters", label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("ATAC")

p3 <- DimPlot(brain, 
              reduction = "wnn.umap", 
              group.by = "seurat_clusters", 
              label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN")


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
  dplyr::select(barcode = rna.bc, celltype)

barcodes.brain <- 
  barcodes.brain %>% 
  left_join(celltype_df)



brain@meta.data <-cbind(brain@meta.data, celltype = barcodes.brain$celltype)

p3 <- DimPlot(brain, 
              reduction = "wnn.umap", 
              group.by = "seurat_clusters", 
              label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN-clusters")

p4 <- DimPlot(brain, 
              reduction = "wnn.umap", 
              group.by = "celltype", 
              label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN-celltypes")

p1 + p2 + p3 & NoLegend() & 
  theme(plot.title = element_text(hjust = 0.5))


p3 + p4 & NoLegend()


## Annotating cell types
table(barcodes.brain$celltype)




## Filter to Gene set of interest ---------------------------------------------

brain_go_term_summary %>% 
  dplyr::filter(grepl("cerebellum", annotated_term)) %$%
  table(annotated_term)


## morphogenesis: "the biological process that causes a cell, tissue or 
##                  organism to develop its shape"

geneset_filtered <- 
  brain_go_term_summary %>% 
  dplyr::filter(grepl("cerebellum dev", annotated_term)) %>% 
  dplyr::filter(evidence_group %in% c("Experimental", "Homology"))


expression_matrix <- 
  GetAssayData(brain, slot = "counts")

## Access normalized count data
umi_counts_norm <- 
  as.matrix(brain@assays$SCT@data) %>% 
  as.data.frame()


## Access cluster labelings
cluster_metadata <- 
  brain@meta.data


table(cluster_metadata$seurat_clusters)


## Example cell cluster examination
umi_counts_subset <- 
  umi_counts_norm %>% 
  dplyr::select(
    cluster_metadata %>% 
      dplyr::filter(seurat_clusters == 10) %>% 
      row.names()
  )



## remove bottom x percent of genes by expression

GENE_PCTILE_CUTOFF <- 0.98

umi_counts_subset$mean_express <- 
  apply(umi_counts_subset, 1, mean)

## Histogram of expression rankings 
umi_counts_subset %>% 
  mutate(in_subset = row.names(umi_counts_subset) %in%
           geneset_filtered$symbol,
         mean_express_rank = rank(-mean_express)) %>%
  dplyr::filter(in_subset==TRUE) %>%
  ggplot(aes(x = mean_express_rank)) + 
  geom_histogram(color="black",fill="white") + 
  ggtitle("Mean expression rank of GO Cerebellum Development\n" %p% 
            "genes in same cell cluster from WNN cell clustering\n" %p% 
            "16,392 Genes") + 
  theme_bw()


umi_cutoff <- 
  quantile(umi_counts_subset$mean_express, 
           GENE_PCTILE_CUTOFF)

umi_counts_random <- 
  umi_counts_subset %>% 
  dplyr::filter(!row.names(umi_counts_subset) %in% 
                  geneset_filtered$symbol & 
                  mean_express > 0) %>% 
  sample_n(100)

umi_counts_subset <- 
  umi_counts_subset %>% 
  dplyr::filter(row.names(umi_counts_subset) %in% 
                  geneset_filtered$symbol) %>% 
  ## no genes with 0 expression
  dplyr::filter(mean_express > 0) %>%
  ## join randomly selected genes
  bind_rows(umi_counts_random) %>%
  dplyr::select(-mean_express)
  

## Gene-Gene coexpression

corr_matrix <-
  rcorr(t(as.matrix(umi_counts_subset)),
        type = "spearman")$r

corr_matrix_df <- 
  corr_matrix %>% 
  as.data.frame() %>% 
  mutate(gene1 = row.names(corr_matrix))

corr_matrix_df <- 
  corr_matrix_df %>%
  gather(gene2, corr, 
         names(corr_matrix_df)[1:(length(names(corr_matrix_df))-1)])

## Too long to render and noisy  
corr_matrix_df %>%
  ggplot(aes(x=gene1, y=gene2, fill=corr)) + geom_tile()+
  scale_fill_gradient(low = "blue", high = "red")


geneset_filtered$symbol %in%
  row.names(umi_counts_norm)

## full genelist: 
full_genes <- brain@assays$RNA@counts@Dimnames[[1]]

## Extract the lacking genes
missing_genes <- 
  geneset_filtered %>% 
  dplyr::filter(!symbol %in% full_genes)


search_genes_in_subset <- 
  geneset_filtered$symbol[geneset_filtered$symbol %in%
                            row.names(umi_counts_subset)]




### Normalize the Laplacian

laplace_matrix <- 
  normalize.laplacian(abs(corr_matrix))

gene_pca_results <-
  prcomp(x = laplace_matrix, 
       ## already centered and normalized
       center = FALSE, scale = FALSE)

gene_pc_df <- 
  gene_pca_results$x %>% 
  as.data.frame() %>% 
  mutate(gene_name = row.names(corr_matrix),
         in_subset = gene_name %in% geneset_filtered$symbol)

gene_pc_df %>% 
  ggplot(aes(x=PC1, y=PC2,
             color = in_subset)) + 
  geom_point() + 
  ggtitle("100 Randomly Expressed Genes, keep all 30 from GO cerebellum dev.\n" %p% 
            "Pick cell cluster of size 106 from joint WNN cell clusters") + 
  theme_bw()
  


## Relevant Pieces of the Pipeline constructed as function --------------------
## ----------------------------------------------------------------------------

SpectralClusterSCSubset <- function(seurat_cluster,
                                    norm_expression_matrix,
                                    cell_cluster_df,
                                    cell_group = "seurat",
                                    gene_sample_method = "random",
                                    pctile = 0.98) {
  ## Example cell cluster examination
  umi_counts_subset <- 
    norm_expression_matrix %>% 
    dplyr::select(
      cell_cluster_df %>% 
        dplyr::filter(seurat_clusters == seurat_cluster) %>% 
        row.names()
    )
  
  
  
  ## remove bottom x percent of genes by expression
  
  GENE_PCTILE_CUTOFF <- 0.98
  
  umi_counts_subset$mean_express <- 
    apply(umi_counts_subset, 1, mean)
  
  
  umi_cutoff <- 
    quantile(umi_counts_subset$mean_express, 
             GENE_PCTILE_CUTOFF)
  
  
  if(gene_sample_method == "random") {
    umi_counts_random <- 
      umi_counts_subset %>% 
      dplyr::filter(!row.names(umi_counts_subset) %in% 
                      geneset_filtered$symbol & 
                      mean_express > 0) %>% 
      sample_n(100)
    
    umi_counts_subset <- 
      umi_counts_subset %>% 
      dplyr::filter(row.names(umi_counts_subset) %in% 
                      geneset_filtered$symbol) %>% 
      ## no genes with 0 expression
      dplyr::filter(mean_express > 0) %>%
      ## join randomly selected genes
      bind_rows(umi_counts_random) %>%
      dplyr::select(-mean_express) 
  }
  
  
  ## Gene-Gene coexpression
  
  corr_matrix <-
    rcorr(t(as.matrix(umi_counts_subset)),
          type = "spearman")$r
  
  
  ### Normalize the Laplacian
  
  laplace_matrix <- 
    normalize.laplacian(abs(corr_matrix))
  
  gene_pca_results <-
    prcomp(x = laplace_matrix, 
           ## already centered and normalized
           center = FALSE, scale = FALSE)
  
  gene_pc_df <- 
    gene_pca_results$x %>% 
    as.data.frame() %>% 
    mutate(gene_name = row.names(corr_matrix),
           in_subset = gene_name %in% geneset_filtered$symbol)
  
  p <-
    gene_pc_df %>% 
    ggplot(aes(x=PC1, y=PC2,
               color = in_subset)) + 
    geom_point() + 
    #ggtitle("100 Randomly Expressed Genes, keep all 30 from GO cerebellum dev.\n" %p% 
    #          "Pick cell cluster of size 106 from joint WNN cell clusters") + 
    theme_bw()
  
  return(list(gene_pc_df = gene_pc_df, pc_plot = p))
  
}




SpectralClusterSCType <- function(celltype,
                                    norm_expression_matrix,
                                    cell_cluster_df,
                                    cell_group = "seurat",
                                    gene_sample_method = "random",
                                    pctile = 0.98) {
  ## Example cell cluster examination
  umi_counts_subset <- 
    norm_expression_matrix %>% 
    dplyr::select(
      cell_cluster_df %>% 
        dplyr::filter(celltype == celltype) %>% 
        row.names()
    )
  
  
  
  ## remove bottom x percent of genes by expression
  
  GENE_PCTILE_CUTOFF <- 0.98
  
  umi_counts_subset$mean_express <- 
    apply(umi_counts_subset, 1, mean)
  
  
  umi_cutoff <- 
    quantile(umi_counts_subset$mean_express, 
             GENE_PCTILE_CUTOFF)
  
  
  if(gene_sample_method == "random") {
    umi_counts_random <- 
      umi_counts_subset %>% 
      dplyr::filter(!row.names(umi_counts_subset) %in% 
                      geneset_filtered$symbol & 
                      mean_express > 0) %>% 
      sample_n(100)
    
    umi_counts_subset <- 
      umi_counts_subset %>% 
      dplyr::filter(row.names(umi_counts_subset) %in% 
                      geneset_filtered$symbol) %>% 
      ## no genes with 0 expression
      dplyr::filter(mean_express > 0) %>%
      ## join randomly selected genes
      bind_rows(umi_counts_random) %>%
      dplyr::select(-mean_express) 
  }
  
  
  ## Gene-Gene coexpression
  
  corr_matrix <-
    rcorr(t(as.matrix(umi_counts_subset)),
          type = "spearman")$r
  
  
  ### Normalize the Laplacian
  
  laplace_matrix <- 
    normalize.laplacian(abs(corr_matrix))
  
  gene_pca_results <-
    prcomp(x = laplace_matrix, 
           ## already centered and normalized
           center = FALSE, scale = FALSE)
  
  gene_pc_df <- 
    gene_pca_results$x %>% 
    as.data.frame() %>% 
    mutate(gene_name = row.names(corr_matrix),
           in_subset = gene_name %in% geneset_filtered$symbol)
  
  p <-
    gene_pc_df %>% 
    ggplot(aes(x=PC1, y=PC2,
               color = in_subset)) + 
    geom_point() + 
    #ggtitle("100 Randomly Expressed Genes, keep all 30 from GO cerebellum dev.\n" %p% 
    #          "Pick cell cluster of size 106 from joint WNN cell clusters") + 
    theme_bw() + 
    ggtitle(celltype)
  
  return(list(gene_pc_df = gene_pc_df, pc_plot = p))
  
}



## Run the operation for each cell cluster:

## Seurat-Generated Clusters
cell_clust_results <- vector(mode = "list",
                             length = length(unique(
                               cluster_metadata$seurat_clusters)))

index <- 1

clust_names <- 
  as.numeric(as.character(sort(unique(cluster_metadata$seurat_clusters))))

for(clust in clust_names) {
  cell_clust_results[[index]] <-
    SpectralClusterSCSubset(seurat_cluster = clust,
                            norm_expression_matrix = umi_counts_norm,
                            cell_cluster_df = cluster_metadata)$pc_plot
  
  index <- index + 1
}
names(cell_clust_results) <- clust_names

ggpubr::ggarrange(plotlist = cell_clust_results,
                  ncol = 5, nrow = 5)


## TODO: what is going on with this one?
cell_clust_results$`1`

table(cluster_metadata$seurat_clusters)

## Share-Seq Annotated Cell Types
length(table(cluster_metadata$celltype))

cell_type_results <- vector(mode = "list",
                             length = length(unique(
                               cluster_metadata$celltype)))

index <- 1

cell_type_names <- 
  as.character(sort(unique(cluster_metadata$celltype)))


for(type in cell_type_names) {
  cell_type_results[[index]] <-
    SpectralClusterSCType(celltype = type,
                            norm_expression_matrix = umi_counts_norm,
                            cell_cluster_df = cluster_metadata)$pc_plot
  
  index <- index + 1
  print(index)
}
names(cell_type_results) <- cell_type_names

ggpubr::ggarrange(plotlist = cell_type_results,
                  ncol = 5, nrow = 4)







