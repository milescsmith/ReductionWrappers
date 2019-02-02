[![Build Status](https://travis-ci.com/milescsmith/dim.reduction.wrappers.svg?branch=master)](https://travis-ci.com/milescsmith/dim.reduction.wrappers)
# dim.reduction.wrappers
R wrappers around dimensionality reduction methods found in Python modules.  Uses the [reticulate](https://github.com/rstudio/reticulate) package to expose functionality.  Additionally provides bridging functions that let these work as drop-in replacements when working with Seurat (verions 2 or 3) and SingleCellExperiment objects.  Currently wraps:
  * [opt-SNE](https://github.com/omiq-ai/Multicore-opt-SNE)
  * [openTSNE](https://github.com/pavlin-policar/openTSNE)
  * [UMAP (Uniform Manifold Approximation and Projection)](https://github.com/lmcinnes/umap)
  * [PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)](https://www.biorxiv.org/content/early/2017/03/24/120378)
  * [PhenoGraph](https://github.com/jacoblevine/PhenoGraph)
