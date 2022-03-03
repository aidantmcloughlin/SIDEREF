splatSimulateAlt <- function(params = newSplatParams(),
                             method = "groups",
                             sparsify = TRUE, verbose = TRUE,
                             group_share_list,
                             group_share_probs,
                             group_share_facLocs,
                             group_share_facScales,
                             group_share_downProb) {
  
  checkmate::assertClass(params, "SplatParams")
  
  #method <- match.arg(method)
  
  if (verbose) {message("Getting parameters...")}
  params <- setParams(params)
  
  
  #params <- splatter::expandParams(params)
  validObject(params)
  
  # Set random seed
  seed <- getParam(params, "seed")
  set.seed(seed)
  
  # Get the parameters we are going to use
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nBatches <- getParam(params, "nBatches")
  batch.cells <- getParam(params, "batchCells")
  nGroups <- getParam(params, "nGroups")
  group.prob <- getParam(params, "group.prob")
  
  
  # Set up name vectors
  if (verbose) {message("Creating simulation object...")}
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  batch.names <- paste0("Batch", seq_len(nBatches))
  group.names <- paste0("Group", seq_len(nGroups))
  
  
  # Create SingleCellExperiment to store simulation
  cells <-  data.frame(Cell = cell.names)
  rownames(cells) <- cell.names
  features <- data.frame(Gene = gene.names)
  rownames(features) <- gene.names
  
  
  sim <- SingleCellExperiment(rowData = features, colData = cells,
                              metadata = list(Params = params))
  
  # Make batches vector which is the index of param$batchCells repeated
  # params$batchCells[index] times
  batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                    b = batch.cells)
  batches <- unlist(batches)
  colData(sim)$Batch <- batch.names[batches]
  
  ## creating Group names for the Cells
  if (method != "single") {
    groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                     replace = TRUE)
    colData(sim)$Group <- factor(group.names[groups], levels = group.names)
  }
  
  
  sim <- splatSimLibSizes(sim, params)
  
  
  if (verbose) {message("Simulating gene means...")}
  sim <- splatSimGeneMeans(sim, params, use.genes.means=NULL)
  # if (nBatches > 1) {
  #   if (verbose) {message("Simulating batch effects...")}
  #   sim <- splatSimBatchEffects(sim, params)
  # }
  
  
  sim <- splatSimBatchCellMeans(sim, params)
  if (method == "single") {
    sim <- splatSimSingleCellMeans(sim, params)
  } else if (method == "groups") {
    
    ## GROUP SIMULATIONS ------------------------------------------------------
    if (verbose) {message("Simulating group DE...")}
    ## DE GENES ------------------------
    sim <- splatSimGroupDEAlt(sim, params,
                              group_share_list,
                              group_share_probs,
                              group_share_facLocs,
                              group_share_facScales,
                              group_share_downProb)
    
    ## CELL MEANS ----------------------
    if (verbose) {message("Simulating cell means...")}
    sim <- splatSimGroupCellMeans(sim, params)
  }
  
  if (verbose) {message("Simulating BCV...")}
  sim <- splatSimBCVMeans(sim, params)
  if (verbose) {message("Simulating counts...")}
  sim <- splatSimTrueCounts(sim, params)
  if (verbose) {message("Simulating dropout (if needed)...")}
  sim <- splatSimDropout(sim, params)
  
  if (sparsify) {
    if (verbose) {message("Sparsifying assays...")}
    assays(sim) <- sparsifyMatrices(assays(sim), auto = TRUE,
                                    verbose = verbose)
  }
  
  if (verbose) {message("Done!")}
  return(sim)
  
}



