## Preprocessing our genomics data sets of interest.


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
save(snare_atac_variable, file = here("output/snareseq_atac_peaks_variable.RData"))


### Gene Activity Matrix ------------------------------------------------------
gene.activities <- GeneActivity(snare)

snare[['ATAC_gene']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(snare) <- "ATAC_gene"

snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)




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
save(snare, file = here("output/snareseq_final.RData"))



set.seed(209827)
snare_sample <- subset(snare, cells = sample(Cells(snare), 1500))

save(snare_sample,
     file = here("output/snare_sample_1500.RData"))

snare_test <-
  snare[, setdiff(colnames(snare), colnames(snare_sample))] 

save(snare_test,
     file = here("output/snare_test.RData"))
