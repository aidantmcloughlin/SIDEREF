## Preprocessing our genomics data sets of interest.



## 10X PBMC Data --------------------------------------------------------------

counts <- Read10X_h5(here("data/PBMC_10X_multiomics/pbmc_granulocyte/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))
fragpath <- here("data/PBMC_10X_multiomics/pbmc_granulocyte/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# create a Seurat object containing the RNA data
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)


DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)


VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)



# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

## TODO: my stopping point
#save(pbmc, counts, annotation, file = here("output/pbmc_intermed.RData"))

# call peaks using MACS2
#peaks <- CallPeaks(pbmc, macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2") ## path to macs2 package installation


## TODO: will need to place this under unified analysis pipeline
## If Narrow Peaks file is generated through terminal:
if(!"peaks" %in% ls()) {
  # read in narrowpeak file
  df <- read.table(
    file = here("data/PBMC_10X_multiomics/pbmc_granulocyte/macs2_peaks.narrowPeak"),
    col.names = c("chr", "start", "end", "name",
                  "score", "strand", "fold_change",
                  "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                  "relative_summit_position")
  )
  peaks <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
}



# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)



### Normalize Gene Expression values ------------------------------------------
DefaultAssay(pbmc) <- "RNA"
pbmc <- FindVariableFeatures(pbmc, nfeatures = 5000)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 30)
#pbmc <- SCTransform(pbmc)
#pbmc <- RunPCA(pbmc)


### Normalize ATAC Peak counts values -----------------------------------------
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)



### Gene Activity Matrix ------------------------------------------------------
gene.activities <- GeneActivity(pbmc)

pbmc[['ATAC_gene']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'ATAC_gene',
  normalization.method = 'LogNormalize'
)


### Cell Annotations ----------------------------------------------------------
data("pbmc.rna")
pbmc.rna


## Get Barcodes and save

labels <- pbmc.rna@meta.data$seurat_annotations
cellnames <- colnames(pbmc.rna@assays$RNA@counts)

labels_xwalk <- data.frame(cellnames, labels)


#cell_names_data <- colnames(pbmc@assays$RNA@counts)
cells_data <- data.frame(cellnames = row.names(pbmc@meta.data))
cells_data <- cells_data %>% left_join(labels_xwalk)

## append back to Seurat object
pbmc@meta.data <- cbind(pbmc@meta.data, cells_data$labels)

save(pbmc, annotation, file = here("output/pbmc_intermed.RData"))







## SNAREseq -------------------------------------------------------------------

set.seed(1234)

# load processed data matrices for each assay
rna <- Read10X(here("data/SNARE_seq/GSE126074_AdBrainCortex_rna/"),
               gene.column = 1)
atac <- Read10X(here("data/SNARE_seq/GSE126074_AdBrainCortex_atac/"),
                gene.column = 1)

fragments <- here("data/SNARE_seq/fragments.sort.bed.gz")

# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)
snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(snare[["ATAC"]]) <- annotations



### Quality control ---------------------------
DefaultAssay(snare) <- "ATAC"
snare <- TSSEnrichment(snare)
snare <- NucleosomeSignal(snare)
snare$blacklist_fraction <- FractionCountsInRegion(
  object = snare,
  assay = 'ATAC',
  regions = blacklist_mm10
)



Idents(snare) <- "all"  # group all cells together, rather than by replicate
VlnPlot(
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)



snare <- subset(
  x = snare,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)
snare

VlnPlot(
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)



### Gene Expression pre-processing ----------------------------------



DefaultAssay(snare) <- "RNA"

snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 30)



### Chromatin Accessibility pre-processing --------------------------

DefaultAssay(snare) <- 'ATAC'

snare <- FindTopFeatures(snare, min.cutoff = 10)
snare <- RunTFIDF(snare)
snare <- RunSVD(snare)

snare_atac_variable <-
  RunTFIDF(snare@assays$ATAC@counts[snare@assays$ATAC@var.features, ])

atac_peaks_variable <- FindVariableFeatures(snare_atac_variable, 
                                            selection.method = "disp")

atac_peaks_variable <- 
  row.names(atac_peaks_variable[order(atac_peaks_variable$mvp.dispersion, 
                                      decreasing = TRUE), 
                                , 
                                drop = FALSE][1:3000, ])

snare_atac_variable <- snare_atac_variable[atac_peaks_variable, ]

## save results
#save(snare_atac_variable, file = here("output/snareseq_atac_peaks_variable.RData"))
#save(snare, file = here("output/snareseq_intermed.RData"))
#load(here("output/snareseq_intermed.RData"))


### Gene Activity Matrix ------------------------------------------------------
gene.activities <- GeneActivity(snare)

snare[['ATAC_gene']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(snare) <- "ATAC_gene"
## TODO: create constant:
snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)

save(snare, file = here("output/snareseq_intermed.RData"))





### Cell Type Labels ------------------------------------------------

# label transfer from Allen brain
allen <- readRDS(here("data/SNARE_seq/allen_brain.RDS"))

# use the RNA assay in the SNARE-seq data for integration with scRNA-seq
DefaultAssay(snare) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = snare,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = snare[['pca']],
  dims = 1:30
)

snare <- AddMetaData(object = snare, metadata = predicted.labels)


# label clusters based on predicted ID
new.cluster.ids <- c(
  "L2/3 IT",
  "L4",
  "L6 IT",
  "L5 CT",
  "L4",
  "L5 PT",
  "Pvalb",
  "Sst",
  "Astro",
  "Oligo",
  "Vip/Lamp5",
  "L6 IT.2",
  "L6b",
  "NP"
)
names(x = new.cluster.ids) <- levels(x = snare)
snare <- RenameIdents(object = snare, new.cluster.ids)
snare$celltype <- Idents(snare)


## save results
#save(snare, file = here("output/snareseq_intermed.RData"))



###############################################################################
###############################################################################
### Training, Testing Splits
###############################################################################
###############################################################################

set.seed(876143987)

load(here("output/snareseq_intermed.RData"))


snare_test_cells <- sort(sample(ncol(snare), floor(ncol(snare) * 0.10)))

snare_train_cells <- setdiff(1:ncol(snare), snare_test_cells)

snare_train <- snare[, snare_train_cells]
snare_test <- snare[, snare_test_cells]


SaveH5Seurat(snare_train, filename = here("output/snareseq_train.h5Seurat"))
SaveH5Seurat(snare_test, filename = here("output/snareseq_test.h5Seurat"))


rm(snare_train, snare_test, snare)
gc()

load(file = here("output/pbmc_intermed.RData"))


pbmc_test_cells <- sort(sample(ncol(pbmc), floor(ncol(pbmc) * 0.10)))

pbmc_train_cells <- setdiff(1:ncol(pbmc), pbmc_test_cells)

pbmc_test <- pbmc[, pbmc_test_cells]
pbmc_train <- pbmc[, pbmc_train_cells]
# couldn't allocate the space for train
SaveH5Seurat(pbmc, filename = here("output/pbmc_train.h5Seurat"))
SaveH5Seurat(pbmc_test, filename = here("output/pbmc_test.h5Seurat"))

