file_based_idx <- lapply(unique(meta_data_full$file_source),
                         function(x) which(meta_data_full$file_source == x))


file_idx_shuf <-
  lapply(file_based_idx, function(x){
    sample(x, replace=FALSE)})


cell_idx_splits_nested <- lapply(file_idx_shuf, function(x){
  split(x,
        cut(seq_along(x), FOLDS,
            labels = FALSE))
  
})

cell_idx_splits <- vector(mode="list", length=FOLDS)

for(f in seq_len(FOLDS)) {
  for(file in seq_len(length(cell_idx_splits_nested))) {
    cell_idx_splits[[f]] = c(cell_idx_splits[[f]], 
                             cell_idx_splits_nested[[file]][[f]])
  }
}

## Take first "REPS" values.
cell_idx_splits <- cell_idx_splits[1:REPS]

