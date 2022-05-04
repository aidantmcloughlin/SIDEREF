### This script downloads the necessary Tabula Muris files to run SIDEREF analysis.
library(here)

SAVE_LOC <- here("data/tab_muris_sc/seurat/")

## If your computer/internet is too slow, you may have to increase this value.
options(timeout=10000)

### URLs of file locations:
DATA_URLS <- c(
  "https://figshare.com/ndownloader/files/13088531",
  "https://figshare.com/ndownloader/files/13088642",
  "https://figshare.com/ndownloader/files/13088849",
  "https://figshare.com/ndownloader/files/13089194",
  "https://figshare.com/ndownloader/files/13089197",
  "https://figshare.com/ndownloader/files/13089323",
  "https://figshare.com/ndownloader/files/13089629",
  "https://figshare.com/ndownloader/files/13089821",
  "https://figshare.com/ndownloader/files/13090478",
  "https://figshare.com/ndownloader/files/13090580",
  "https://figshare.com/ndownloader/files/13090877",
  "https://figshare.com/ndownloader/files/13091135"
  )

### Corresponding file names:
FILE_NAMES <- c(
  "droplet_Bladder_seurat_tiss.Robj",
  "droplet_Heart_and_Aorta_seurat_tiss.Robj",
  "droplet_Kidney_seurat_tiss.Robj",
  "droplet_Limb_Muscle_seurat_tiss.Robj",
  "droplet_Liver_seurat_tiss.Robj",
  "droplet_Lung_seurat_tiss.Robj",
  "droplet_Mammary_Gland_seurat_tiss.Robj",
  "droplet_Marrow_seurat_tiss.Robj",
  "droplet_Spleen_seurat_tiss.Robj",
  "droplet_Thymus_seurat_tiss.Robj",
  "droplet_Tongue_seurat_tiss.Robj",
  "droplet_Trachea_seurat_tiss.Robj"
)

for(d in 7:length(DATA_URLS)) { 
  download.file(DATA_URLS[d], 
                paste0(SAVE_LOC, FILE_NAMES[d]))
  }



