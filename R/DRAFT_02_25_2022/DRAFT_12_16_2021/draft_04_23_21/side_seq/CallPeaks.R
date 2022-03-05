## Call Peaks Usage

## args for CallPeaks() function
macs2.path = "MACS2PATH "
outdir = tempdir()
broad = FALSE
effective.genome.size = 2.7e9
extsize = 200
shift = -extsize/2
additional.args = NULL
name = "macs2"
cleanup = TRUE
verbose = TRUE
  
frags <- Fragments(object = pbmc)
# get all fragment file paths
allfragpaths <- sapply(X = frags, FUN = GetFragmentData, slot = "path")
allfragpaths <- Reduce(f = paste, x = allfragpaths)

object = allfragpaths


# if list of paths given, collapse to a single space-separated string
if (length(x = object) > 1) {
  object <- Reduce(f = paste, x = object)
}

broadstring <- ifelse(test = broad, yes = " --broad ", no = "")

cmd <- paste0(
  macs2.path,
  " callpeak -t ",
  object,
  " -g ",
  as.character(x = effective.genome.size),
  broadstring,
  " -f BED --nomodel --extsize ",
  as.character(x = extsize),
  " --shift ",
  as.character(x = shift),
  " -n ",
  as.character(x = name),
  " --outdir ",
  outdir,
  " ",
  additional.args
)

# call macs2
system(
  command = cmd,
  wait = TRUE,
  ignore.stderr = !verbose,
  ignore.stdout = !verbose
)

if (broad) {
  # read in broadpeak
  df <- read.table(
    file = paste0(outdir, .Platform$file.sep, name, "_peaks.broadPeak"),
    col.names = c("chr", "start", "end", "name",
                  "score", "strand", "fold_change",
                  "neg_log10pvalue_summit", "neg_log10qvalue_summit")
  )
  files.to.remove <- paste0(
    name,
    c("_peaks.broadPeak", "_peaks.xls", "_peaks.gappedPeak")
  )
} else {
  # read in narrowpeak file
  df <- read.table(
    file = paste0(outdir, .Platform$file.sep, name, "_peaks.narrowPeak"),
    col.names = c("chr", "start", "end", "name",
                  "score", "strand", "fold_change",
                  "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                  "relative_summit_position")
  )
  files.to.remove <- paste0(
    name,
    c("_peaks.narrowPeak", "_peaks.xls", "_summits.bed")
  )
}

gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
if (cleanup) {
  files.to.remove <- paste0(outdir, .Platform$file.sep, files.to.remove)
  for (i in files.to.remove) {
    if (file.exists(i)) {
      file.remove(i)
    }
  }
}
