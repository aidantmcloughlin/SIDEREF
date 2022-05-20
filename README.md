# SIDEREF
SIDEREF dissimilarity matrix for examining scRNA seq global relationships.

*Note: Package development for key functions in this repository is in-progress.*

To conserve space on the repository, the computationally expensive 
outputs run on our HPC are available as a compressed file [here](https://figshare.com/projects/SIDEREF_results_objects/139609). 
The placeholder directory/ies for each downloaded output folder is ``output``.



You may reproduce output figures and tables that use these data by running:
  ``library(here)``
  
  ``source(here("R/post_hoc_analysis.R"))``
  

In order to reproduce all analysis, you will first need to download the following 
files from *Tabula Muris*  [here](https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263). The placeholder directory for these raw data files is ``data/tab_muris_sc/seurat``.

A brief description of this data is located  [here](https://tabula-muris.ds.czbiohub.org).

We provide an R script to automatically download and place the data files:
  ``library(here)``
  
  ``source(here("R/download_data.R"))``

The files we download for our analysis are as follows:

* droplet_Bladder_seurat_tiss.Robj
* droplet_Heart_and_Aorta_seurat_tiss.Robj
* droplet_Kidney_seurat_tiss.Robj
* droplet_Limb_Muscle_seurat_tiss.Robj
* droplet_Liver_seurat_tiss.Robj
* droplet_Lung_seurat_tiss.Robj
* droplet_Mammary_Gland_seurat_tiss.Robj
* droplet_Marrow_seurat_tiss.Robj
* droplet_Spleen_seurat_tiss.Robj
* droplet_Thymus_seurat_tiss.Robj
* droplet_Tongue_seurat_tiss.Robj
* droplet_Trachea_seurat_tiss.Robj

You may also be interested in exploring the FACS scRNA seq files from this link.

Once the data is downloaded and placed in the specified folder, the entire analysis may be run using the following two commands, after opening ``sc_multiome.Rproj``:

  ``library(here)``
  
  ``source(here("R/main.R"))``
  
Package dependencies are located in ``R/libraries.R``.
