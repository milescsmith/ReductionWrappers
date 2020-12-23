<!-- badges: start -->
  [![R build status](https://github.com/milescsmith/ReductionWrappers/workflows/pkgdown/badge.svg)](https://github.com/milescsmith/ReductionWrappers/actions)
[![Build Status](https://travis-ci.com/milescsmith/ReductionWrappers.svg?branch=master)](https://travis-ci.com/milescsmith/ReductionWrappers)
<!-- badges: end -->

# ReductionWrappers

R wrappers around dimensionality reduction methods found in Python modules.  Uses the [reticulate](https://github.com/rstudio/reticulate) package to expose functionality.  Additionally provides bridging functions that let these work as drop-in replacements when working with Seurat (verions 3) and SingleCellExperiment objects.  Currently wraps:
  * [ForceAtlas2](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679), using the Python implementation found [here](https://github.com/bhargavchippada/forceatlas2).  Currently works for Seurat objects only (and will remain so until I can find a fast algorithm for computing simultaneous nearest neighbors for SingleCellExperiment objects).
  * [opt-SNE](https://github.com/omiq-ai/Multicore-opt-SNE)
  * [openTSNE](https://github.com/pavlin-policar/openTSNE)
  * [UMAP (Uniform Manifold Approximation and Projection)](https://github.com/lmcinnes/umap)
  * [PaCMAP (Pairwise Controlled Manifold Approximation)](https://github.com/YingfanWang/PaCMAP)
  * [PAGA (Partition-Based Graph Abstraction)](https://github.com/theislab/paga) through Scanpy.
  * [PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)](https://www.biorxiv.org/content/early/2017/03/24/120378)
  * [PhenoGraph](https://github.com/jacoblevine/PhenoGraph)
