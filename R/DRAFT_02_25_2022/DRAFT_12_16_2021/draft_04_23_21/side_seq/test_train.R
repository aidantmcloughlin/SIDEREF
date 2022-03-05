################ Build test train splits


set.seed(876143987)

load("snareseq_final.RData")


snare_test_cells <- sort(sample(ncol(snare), floor(ncol(snare) * 0.10)))

snare_train_cells <- setdiff(1:ncol(snare), snare_test_cells)

snare_train <- snare[, snare_train_cells]
snare_test <- snare[, snare_test_cells]


SaveH5Seurat(snare_train, filename = "output/snareseq_train.h5Seurat", overwrite = TRUE)
SaveH5Seurat(snare_test, filename = "output/snareseq_test.h5Seurat", overwrite = TRUE)


rm(snare_train, snare_test, snare)
gc()

load(file = ("pbmc_intermed.RData"))


pbmc_test_cells <- sort(sample(ncol(pbmc), floor(ncol(pbmc) * 0.10)))

pbmc_train_cells <- setdiff(1:ncol(pbmc), pbmc_test_cells)

pbmc_test <- pbmc[, pbmc_test_cells]
pbmc_train <- pbmc[, pbmc_train_cells]
# couldn't allocate the space for train
SaveH5Seurat(pbmc_train, filename = "output/pbmc_train.h5Seurat", overwrite = TRUE)
SaveH5Seurat(pbmc_test, filename = "output/pbmc_test.h5Seurat", overwrite = TRUE)

