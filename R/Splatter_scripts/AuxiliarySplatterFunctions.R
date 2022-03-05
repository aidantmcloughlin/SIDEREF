### Pasting internal functions

splatSimLibSizes <- function(sim, params) {
  
  nCells <- getParam(params, "nCells")
  lib.loc <- getParam(params, "lib.loc")
  lib.scale <- getParam(params, "lib.scale")
  lib.norm <- getParam(params, "lib.norm")
  if (lib.norm) {
    exp.lib.sizes <- rnorm(nCells, lib.loc, lib.scale)
    min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
    exp.lib.sizes[exp.lib.sizes < 0] <- min.lib/2
  }
  else {
    exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
  }
  colData(sim)$ExpLibSize <- exp.lib.sizes
  return(sim)
}




splatSimBatchCellMeans <- function(sim, params) {
  
  nBatches <- getParam(params, "nBatches")
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  gene.means <- rowData(sim)$GeneMean
  if (nBatches > 1) {
    batches <- colData(sim)$Batch
    batch.names <- unique(batches)
    batch.facs.gene <- rowData(sim)[, paste0("BatchFac", 
                                             batch.names)]
    batch.facs.cell <- as.matrix(batch.facs.gene[, as.numeric(factor(batches))])
  }
  else {
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    batch.facs.cell <- matrix(1, ncol = nCells, nrow = nGenes)
  }
  batch.means.cell <- batch.facs.cell * gene.means
  colnames(batch.means.cell) <- cell.names
  rownames(batch.means.cell) <- gene.names
  assays(sim)$BatchCellMeans <- batch.means.cell
  return(sim)

}



splatSimGeneMeans <- function(sim, params, use.genes.means) {
  
  nGenes <- getParam(params, "nGenes")
  mean.shape <- getParam(params, "mean.shape")
  mean.rate <- getParam(params, "mean.rate")
  out.prob <- getParam(params, "out.prob")
  out.facLoc <- getParam(params, "out.facLoc")
  out.facScale <- getParam(params, "out.facScale")
  base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
  outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc, 
                                  out.facScale)
  median.means.gene <- median(base.means.gene)
  outlier.means <- median.means.gene * outlier.facs
  is.outlier <- outlier.facs != 1
  means.gene <- base.means.gene
  means.gene[is.outlier] <- outlier.means[is.outlier]
  rowData(sim)$BaseGeneMean <- base.means.gene
  rowData(sim)$OutlierFactor <- outlier.facs
  rowData(sim)$GeneMean <- means.gene
  return(sim)
}




splatSimTrueCounts <- function(sim, params) {
  
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  nGenes <- getParam(params, "nGenes")
  nCells <- getParam(params, "nCells")
  cell.means <- assays(sim)$CellMeans
  true.counts <- matrix(rpois(as.numeric(nGenes) * as.numeric(nCells), 
                              lambda = cell.means), nrow = nGenes, ncol = nCells)
  colnames(true.counts) <- cell.names
  rownames(true.counts) <- gene.names
  assays(sim)$TrueCounts <- true.counts
  return(sim)
}


### GROUP FUNCTIONS -----------------------------------------------------------
splatSimGroupDE <- function(sim, params) {

  nGenes <- getParam(params, "nGenes")
  nGroups <- getParam(params, "nGroups")
  de.prob <- getParam(params, "de.prob")
  de.downProb <- getParam(params, "de.downProb")
  de.facLoc <- getParam(params, "de.facLoc")
  de.facScale <- getParam(params, "de.facScale")
  means.gene <- rowData(sim)$GeneMean
  
  
  
  for (idx in seq_len(nGroups)) {
    de.facs <- getLNormFactors(nGenes, de.prob[idx], de.downProb[idx], 
                               de.facLoc[idx], de.facScale[idx])
    group.means.gene <- means.gene * de.facs
    rowData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
  }
  return(sim)
}


splatSimGroupDEAlt <- function(sim, params,
                               group_share_list,
                               group_share_probs,
                               group_share_facLocs,
                               group_share_facScales,
                               group_share_downProb) {
  
  nGenes <- getParam(params, "nGenes")
  nGroups <- getParam(params, "nGroups")
  de.prob <- getParam(params, "de.prob")
  de.downProb <- getParam(params, "de.downProb")
  de.facLoc <- getParam(params, "de.facLoc")
  de.facScale <- getParam(params, "de.facScale")
  means.gene <- rowData(sim)$GeneMean

  ### assign standard groups MFs
  for (idx in seq_len(nGroups)) {
    de.facs <- getLNormFactors(nGenes, de.prob[idx], de.downProb[idx], 
                               de.facLoc[idx], de.facScale[idx])
    group.means.gene <- means.gene * de.facs
    rowData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
  }
  
  ### assign shared DE genes btw global groups
  for(i in seq_len(length(group_share_list))) {
    shared_groups <- group_share_list[[i]]
    de_cols <- which(names(rowData(sim)) %in% paste0("DEFacGroup", shared_groups))
    genes_non_de <- which(apply(rowData(sim)[, de_cols], 1, 
                                function(x) return(sum(abs(1 - x) > 1e-4) == 0)))
    
    sel.prob  <- group_share_probs[i]
    fac.loc   <- group_share_facLocs[i]
    fac.scale <- group_share_facScales[i]
    neg.prob  <- group_share_downProb[i]
    
    ## select from pool and compute multiplicative factors 
    is.selected <- as.logical(rbinom(seq_len(length(genes_non_de)), 
                                     1, sel.prob))
    n.selected <- sum(is.selected)
    dir.selected <- (-1)^rbinom(n.selected, 1, neg.prob)
    facs.selected <- matrix(rlnorm(n.selected * length(shared_groups), 
                                   fac.loc, fac.scale), ncol = length(shared_groups))
    ## reverse direction if needed.
    facs.selected <- apply(facs.selected, 2, function(x) ifelse(x < 1, 1/x, x))
    
    final_fac_mat <- apply(facs.selected, 2, function(x)x^dir.selected)
    rowData(sim)[genes_non_de[is.selected], de_cols] <- final_fac_mat
    
    ## track these genes
    rowData(sim)[["shared_genes" %p% paste(shared_groups, collapse = "_")]] <- 
      seq_len(nGenes) %in% genes_non_de[is.selected]
    

  }
  
  return(sim)
}



getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale) {
  
  is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
  n.selected <- sum(is.selected)
  dir.selected <- (-1)^rbinom(n.selected, 1, neg.prob)
  facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
  dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
  factors <- rep(1, n.facs)
  factors[is.selected] <- facs.selected^dir.selected
  return(factors)
}




splatSimGroupCellMeans <- function(sim, params) {
  nGroups <- getParam(params, "nGroups")
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  groups <- colData(sim)$Group
  group.names <- levels(groups)
  exp.lib.sizes <- colData(sim)$ExpLibSize
  batch.means.cell <- assays(sim)$BatchCellMeans
  group.facs.gene <- rowData(sim)[, paste0("DEFac", group.names)]
  cell.facs.gene <- as.matrix(group.facs.gene[, paste0("DEFac", 
                                                       groups)])
  cell.means.gene <- batch.means.cell * cell.facs.gene
  cell.props.gene <- t(t(cell.means.gene)/colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  assays(sim)$BaseCellMeans <- base.means.cell
  return(sim)
}




############  -----------------------------------------------------------------
### Grabbing the actual UMI counts --------------------------------------------
############  -----------------------------------------------------------------

splatSimBCVMeans <- function(sim, params) {
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  nGenes <- getParam(params, "nGenes")
  nCells <- getParam(params, "nCells")
  bcv.common <- getParam(params, "bcv.common")
  bcv.df <- getParam(params, "bcv.df")
  base.means.cell <- assays(sim)$BaseCellMeans
  if (is.finite(bcv.df)) {
    bcv <- (bcv.common + (1/sqrt(base.means.cell))) * sqrt(bcv.df/rchisq(nGenes, 
                                                                         df = bcv.df))
  }
  else {
    warning("'bcv.df' is infinite. This parameter will be ignored.")
    bcv <- (bcv.common + (1/sqrt(base.means.cell)))
  }
  means.cell <- matrix(rgamma(as.numeric(nGenes) * as.numeric(nCells), 
                              shape = 1/(bcv^2), scale = base.means.cell * (bcv^2)), 
                       nrow = nGenes, ncol = nCells)
  colnames(means.cell) <- cell.names
  rownames(means.cell) <- gene.names
  assays(sim)$BCV <- bcv
  assays(sim)$CellMeans <- means.cell
  return(sim)
}


splatSimTrueCounts <- function(sim, params) {
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  nGenes <- getParam(params, "nGenes")
  nCells <- getParam(params, "nCells")
  cell.means <- assays(sim)$CellMeans
  true.counts <- matrix(rpois(as.numeric(nGenes) * as.numeric(nCells), 
                              lambda = cell.means), nrow = nGenes, ncol = nCells)
  colnames(true.counts) <- cell.names
  rownames(true.counts) <- gene.names
  assays(sim)$TrueCounts <- true.counts
  return(sim)
}





############  -----------------------------------------------------------------
### DROPOUT FUNCTION ----------------------------------------------------------
############  -----------------------------------------------------------------

splatSimDropout <- function(sim, params) {
  dropout.type <- getParam(params, "dropout.type")
  true.counts <- assays(sim)$TrueCounts
  dropout.mid <- getParam(params, "dropout.mid")
  dropout.shape <- getParam(params, "dropout.shape")
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nBatches <- getParam(params, "nBatches")
  nGroups <- getParam(params, "nGroups")
  cell.means <- assays(sim)$CellMeans
  switch(dropout.type, experiment = {
    if ((length(dropout.mid) != 1) || length(dropout.shape) != 
        1) {
      stop("dropout.type is set to 'experiment' but dropout.mid ", 
           "and dropout.shape aren't length 1")
    }
    dropout.mid <- rep(dropout.mid, nCells)
    dropout.shape <- rep(dropout.shape, nCells)
  }, batch = {
    if ((length(dropout.mid) != nBatches) || length(dropout.shape) != 
        nBatches) {
      stop("dropout.type is set to 'batch' but dropout.mid ", 
           "and dropout.shape aren't length equal to nBatches ", 
           "(", nBatches, ")")
    }
    batches <- as.numeric(factor(colData(sim)$Batch))
    dropout.mid <- dropout.mid[batches]
    dropout.shape <- dropout.shape[batches]
  }, group = {
    if ((length(dropout.mid) != nGroups) || length(dropout.shape) != 
        nGroups) {
      stop("dropout.type is set to 'group' but dropout.mid ", 
           "and dropout.shape aren't length equal to nGroups ", 
           "(", nGroups, ")")
    }
    if ("Group" %in% colnames(colData(sim))) {
      groups <- as.numeric(colData(sim)$Group)
    } else {
      stop("dropout.type is set to 'group' but groups have not ", 
           "been simulated")
    }
    dropout.mid <- dropout.mid[groups]
    dropout.shape <- dropout.shape[groups]
  }, cell = {
    if ((length(dropout.mid) != nCells) || length(dropout.shape) != 
        nCells) {
      stop("dropout.type is set to 'cell' but dropout.mid ", 
           "and dropout.shape aren't length equal to nCells ", 
           "(", nCells, ")")
    }
  })
  if (dropout.type != "none") {
    drop.prob <- vapply(seq_len(nCells), function(idx) {
      eta <- log(cell.means[, idx])
      return(logistic(eta, x0 = dropout.mid[idx], k = dropout.shape[idx]))
    }, as.numeric(seq_len(nGenes)))
    keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob), 
                   nrow = nGenes, ncol = nCells)
    counts <- true.counts * keep
    colnames(drop.prob) <- cell.names
    rownames(drop.prob) <- gene.names
    colnames(keep) <- cell.names
    rownames(keep) <- gene.names
    assays(sim)$DropProb <- drop.prob
    assays(sim)$Dropout <- !keep
  }
  else {
    counts <- true.counts
  }
  BiocGenerics::counts(sim) <- counts
  return(sim)
}


sparsifyMatrices <- function(matrix.list, auto = TRUE, threshold = 0.95, verbose = TRUE) {
  if (!auto) {
    if (verbose) {
      message("Converting all matrices to sparse format")
    }
    matrix.list <- lapply(matrix.list, as, Class = "dgCMatrix")
    return(matrix.list)
  }
  if (verbose) {
    message("Automatically converting to sparse matrices, ", 
            "threshold = ", threshold)
  }
  for (mat.name in names(matrix.list)) {
    mat <- matrix.list[[mat.name]]
    prop.zero <- sum(mat == 0)/length(mat)
    if (is.integer(mat)) {
      size.factor <- 3 - 3 * prop.zero
    }
    else if (is.numeric(mat)) {
      size.factor <- 1.5 - 1.5 * prop.zero
    }
    else if (is(mat, "dgCMatrix")) {
      if (verbose) {
        message("Skipping '", mat.name, "' as it is already a dgCMatrix")
      }
      next
    }
    else {
      warning("matrix '", mat.name, "' is class '", 
              class(mat), "', unable to estimate size reduction factor")
      size.factor <- NA
    }
    if (is.na(size.factor) | size.factor < threshold) {
      if (verbose) {
        message("Converting '", mat.name, "' to sparse matrix: ", 
                "estimated sparse size ", round(size.factor, 
                                                2), " * dense matrix")
      }
      matrix.list[[mat.name]] <- as(mat, "dgCMatrix")
    }
    else {
      if (verbose) {
        message("Skipping '", mat.name, "': estimated sparse size ", 
                round(size.factor, 2), " * dense matrix")
      }
    }
  }
  return(matrix.list)
}




###############################################################################
### Dropout
###############################################################################


logistic <-function (x, x0, k) {
  1/(1 + exp(-k * (x - x0)))
}
