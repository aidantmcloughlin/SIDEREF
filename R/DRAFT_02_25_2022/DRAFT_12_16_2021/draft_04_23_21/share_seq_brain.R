
ct=read.table(here('data/share_seq/brain/celltype_brain.txt'), header = T)
# first column: atac barcode
# second column: rna barcode
# third column: annotated cell type from the paper

rna.count=read.table(here('data/share_seq/brain/GSM4156610_brain.rna.counts.txt.gz'), 
                     header = T, row.names = 1)
rna.count[1:5, 1:5]

# match the rna.count with the cell type labels
any(is.na(match(ct$rna.bc, colnames(rna.count))))
rna.count=rna.count[,match(ct$rna.bc, colnames(rna.count))]

# The barcode file is the same as the annotation ct; omitted
peak.barcodes <- scan(here('data/share_seq/brain/GSM4156599_brain.barcodes.txt.gz'), what="")
all(ct$rna.bc==peak.barcodes)
rm(peak.barcodes)

# Read in the peak matrix: stored as Market matrix
# https://atlas.gs.washington.edu/mouse-atac/docs/
library(Matrix)
library(GenomicRanges)
peak.count= readMM(here("data/share_seq/brain/GSM4156599_brain.counts.txt.gz"))
dim(peak.count)
peak.bed=read.delim(here('data/share_seq/brain/GSM4156599_brain.peaks.bed.gz'), header=FALSE)
peak.granges=GRanges(seqnames=peak.bed$V1, ranges=IRanges(st=peak.bed$V2, end=peak.bed$V3))
rm(peak.bed)
rownames(peak.count)=paste(peak.granges)
colnames(peak.count)=ct$rna.bc # we will use the rna barcode from here and on

dim(ct) # barcode and cell type
dim(peak.count) # peak count
length(peak.granges) # GRanges for peaks
dim(rna.count) # rna count


library(Seurat)
library(Signac)
# library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79) # For mouse brain
library(dplyr)
library(ggplot2)

# # None is returned with the command below
# rownames(rna.count)[grep('^mt-', rownames(rna.count), ignore.case = TRUE)]
# # Try another approach https://www.michaelchimenti.com/2019/03/calculate-mitochondrial-for-mouse-scrna-seq/
# library(tidyverse)
# mouse_mito = as.tibble(read.csv("Mouse.MitoCarta2.0.txt", header = TRUE, sep='\t'))
# mouse_mito = mouse_mito %>% select(c(Symbol, MCARTA2.0_score)) %>% slice(1:100)
# mito.genes = as.character(mouse_mito$Symbol)
# mito.genes = mito.genes[mito.genes %in% rownames(rna.count)]
# rownames(rna.count)[match(mito.genes, rownames(rna.count))]=paste('MT-', mito.genes, sep='')
# # The results do not look good. The list may not be as accurate.

# Create Seurat object
brain <- CreateSeuratObject(counts = rna.count)
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")


# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.use <- seqnames(peak.granges) %in% standardChromosomes(peak.granges)
peak.count <- peak.count[as.vector(grange.use), ]
peak.granges <- peak.granges[as.vector(grange.use)]
rm(grange.use)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Yuriko mentioned that there is issue loading the fragment file
# Leave the fragment file out for now.
# frag.file <- "GSM4156599_brain.atac.fragments.bed.gz"
chrom_assay <- CreateChromatinAssay(
  counts = peak.count,
  sep = c(":", "-"),
  genome = 'mm10',
  # fragments = frag.file,
  min.cells = 10,
  annotation = annotations)

brain[["ATAC"]] <- chrom_assay
brain@assays$RNA
brain@assays$ATAC


## perform basic qc based on the number of detected molecules for each modality
## as well as mitochondrial percentage
VlnPlot(brain, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
brain <- subset(
  x = brain,
  subset = nCount_ATAC < 3000 &
    nCount_ATAC > 100 &
    nCount_RNA < 25000 &
    nCount_RNA > 500 &
    percent.mt < 20
)
VlnPlot(brain, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()



## perform pre-processing and dimensional reduction on both assays independently
## using standard approaches for rna and atac-seq data
DefaultAssay(brain) <- "RNA"
brain <- SCTransform(brain, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(brain) <- "ATAC"
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(brain)
brain <- RunUMAP(brain, reduction = 'lsi', dims = 2:50,
                reduction.name = "umap.atac",
                reduction.key = "atacUMAP_")
## calculate a wnn graph, representing a weighted combination of rna and atac-seq modalities
## use this graph for umap visualization and clustering
brain <- FindMultiModalNeighbors(brain, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
brain <- RunUMAP(brain, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
brain <- FindClusters(brain, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
table(brain$seurat_clusters)

p1 <- DimPlot(brain, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(brain, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(brain, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
save(brain, file='brain.rda')



# We will now get the gene activity matrix
# https://satijalab.org/signac/articles/brain_vignette.html
DefaultAssay(brain) <- "ATAC"
# gene.activities <- GeneActivity(brain) # This would not run because it requires fragment file to be loaded.

# Below is a short-cut but is not very precise: we only aggregate reads in the peak regions
# But there are reads off peak regions but within gene bodies that are defined
annotation <- Annotation(object = brain)
transcripts <- Signac:::CollapseToLongestTranscript(ranges = annotation)
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"]
transcripts <- Extend(x = transcripts, upstream = 2000, 
                      downstream = 0)
gene.activities=matrix(nrow=length(transcripts), ncol=ncol(brain))
rownames(gene.activities)=transcripts$gene_name
colnames(gene.activities)=colnames(brain)
peak.count=brain@assays$ATAC@counts
peak.granges=brain@assays$ATAC@ranges
for(i in 1:nrow(gene.activities)){
  if(i %%200 ==0) cat(i,' ')
  peak.index=which(countOverlaps(peak.granges, transcripts[i])>0)
  gene.activities[i,]=apply(peak.count[peak.index,, drop=FALSE],2, sum)
}
gene.activities=gene.activities[apply(gene.activities,1,sum)>10,]
gene.activities=gene.activities[!duplicated(rownames(gene.activities)),] # remove duplicate genes
gene.activities=gene.activities[-which(rownames(gene.activities)==''),]

brain[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
## NOTE: WE WILL STILL NEED TO LOOK INTO THE FRAGMENT FILE;
## This calculation is not accurate since we only focus on the peak regions

# Perform normalization/scaling of the gene activity matrix
DefaultAssay(brain) <- "ACTIVITY"
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes)
brain@assays$ACTIVITY@scale.data[1:5,1:5]
save(brain, file='brain.rda')



