# dim.reduction.wrappers
R wrappers around dimensionality reduction methods found in Python modules.  Uses the [reticulate](https://github.com/rstudio/reticulate) package to expose functionality.  Additionally provides bridging functions that let these work as drop-in replacements when working with Seurat objects.  Currently wraps:
  * [Multicore-TSNE](https://github.com/DmitryUlyanov/Multicore-TSNE)
  * [UMAP (Uniform Manifold Approximation and Projection)](https://github.com/lmcinnes/umap)
  * [PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)](https://www.biorxiv.org/content/early/2017/03/24/120378)
  * [PhenoGraph](https://github.com/jacoblevine/PhenoGraph)
    * [Deep count autoencoder](https://github.com/theislab/dca)
  
To be added shortly:
  * [UNCURL](https://github.com/yjzhang/uncurl_python)
