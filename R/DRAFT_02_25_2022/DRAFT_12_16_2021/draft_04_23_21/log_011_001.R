# log_011_001.r
# e18

rm(list=ls())
setwd("/Users/yurikoharigaya/Documents/uncch/jiang_lab/project_alignment/log/log_011")

# https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html

library(Seurat)
library(Signac)
BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)

## create a seurat object based on the gene expression data
## and then add in the atac-seq data as a second assay.

file <- "e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5"
path <- paste(dir.in, file, sep="")
inputdata <- Read10X_h5(path)
# Genome matrix has multiple modalities, returning a list of matrices for this genome

list.files(dir.out)
file <- "inputdata.rda"
path <- paste(dir.out, file, sep="")
save(inputdata, file=path)


list.files(dir.out)
file <- "inputdata.rda"
path <- paste(dir.out, file, sep="")
load(path)
rna_counts <- inputdata$`Gene Expression`
atac_counts <- inputdata$Peaks

e18 <- CreateSeuratObject(counts = rna_counts)
e18[["percent.mt"]] <- PercentageFeatureSet(e18, pattern = "^MT-")
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
dir <- "/Users/yurikoharigaya/Dropbox/download_database/10xgenomics/e18_mouse_brain_fresh/"
file <- "e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz"
path <- paste(dir, file, sep="")
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = path,
  min.cells = 10,
  annotation = annotations
)
e18[["ATAC"]] <- chrom_assay

## perform basic qc based on the number of detected molecules for each modality
## as well as mitochondrial percentage
VlnPlot(e18, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
e18 <- subset(
  x = e18,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)
e18
# An object of class Seurat 
# 170928 features across 3867 samples within 2 assays 
# Active assay: RNA (32285 features, 0 variable features)
# 1 other assay present: ATAC

## perform pre-processing and dimensional reduction on both assays independently
## using standard approaches for rna and atac-seq data.
DefaultAssay(e18) <- "RNA"
e18 <- SCTransform(e18, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(e18) <- "ATAC"
e18 <- FindTopFeatures(e18, min.cutoff = 5)
e18 <- RunTFIDF(e18)
e18 <- RunSVD(e18)
e18 <- RunUMAP(e18, reduction = 'lsi', dims = 1:50,
                reduction.name = "umap.atac",
                reduction.key = "atacUMAP_")
## calculate a wnn graph, representing a weighted combination of rna and atac-seq modalities

load(here("data/EmbryonicMouseBrain_10x_multiomics/e18.rda"))



e18 <- 
  FindMultiModalNeighbors(e18, 
                          reduction.list = list("pca", "lsi"), 
                          dims.list = list(1:50, 2:40))


e18 <- RunUMAP(e18, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")

### ATM: they are using the SLM algorithm here.
e18 <- FindClusters(e18, 
                    graph.name = "wsnn", 
                    algorithm = 3, verbose = FALSE)





p1 <- DimPlot(e18, 
              reduction = "umap.rna", 
              group.by = "seurat_clusters", label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("RNA")

p2 <- DimPlot(e18, 
              reduction = "umap.atac", 
              group.by = "seurat_clusters", label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("ATAC")

## joint UMAP with autocluster labels
p3 <- DimPlot(e18, 
              reduction = "wnn.umap", 
              group.by = "seurat_clusters", 
              label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN Auto-Cluster")

p4 <- DimPlot(e18, 
              reduction = "wnn.umap",
              group.by = "celltype",
              label = TRUE, 
              label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN Cell Names")

p1 + p2 + p3 & NoLegend() & 
  theme(plot.title = element_text(hjust = 0.5))


p3 & NoLegend()
p4 & NoLegend()



## Cell types column

e18_mouse_brain_fresh_5k_atac_peak_annotation.tsv



