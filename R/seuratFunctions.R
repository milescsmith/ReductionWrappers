#' dim.from.python
#'
#' Helper function to assist entering dimensional reduction data from Python
#' reduction methods
#'
#' @param seuratObj
#' @param python.dataframe Dataframe returned by a Python function
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...)
#' @param assay.used Assay from which the data that is dimensionally reduced comes
#'
#' @importFrom Seurat CreateDimReducObject
#' @importFrom glue glue
#'
#' @return A Seurat object with the dataframe stored in the
#'   seuratObj@dr$reduction.use slot
#' @export
#'
#' @examples
dim.from.python <- function(seuratObj, python.dataframe, reduction.save, assay.used){
  dim.xfer <- as.matrix(python.dataframe)
  rownames(dim.xfer) <- colnames(seuratObj)
  reduction.data <- CreateDimReducObject(embeddings = dim.xfer,
                                         assay = assay.used,
                                         key = glue("{reduction.save}_"))
  seuratObj[[reduction.save]] <- reduction.data
  return(seuratObj)
}

#' python.dim.reduction.bridge
#'
#' Generalized helper function that pulls the data from a Seurat object, passes
#' the dataframe to a Python function and places the resulting dataframe in the
#' appropriate slot
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...)
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...)
#' @param function.use Dimensional reduction function to call.
#' @param ... Extra parameters to pass to the dimensional reduction function.
#'
#' @importFrom Seurat GetDimReduction
#'
#' @return
#' @export
#'
#' @examples
python.dim.reduction.bridge <- function(seuratObj,
                                        reduction.use,
                                        reduction.save,
                                        function.use,
                                        ...){

  if (reduction.use %in% names(seuratObj)) {
    cell.embeddings <- Embeddings(object[[reduction]])
    assay <- DefaultAssay(object = object[[reduction]])
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }

  dim.df <- function.use(cell.embeddings, ...)
  seuratObj <- dim.from.python(seuratObj = seuratObj,
                               python.dataframe = dim.df,
                               reduction.save = reduction.save,
                               assay.used = assay)
  return(seuratObj)
}

#' GoMCtSNE
#'
#' Run tSNE projection on a Seurat object using the MulticoreTSNE function
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). (default: pca)
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). (default: tsne)
#' @param ... Extra parameters to pass to the multicoreTSNE function.
#'
#' @return Seurat object with the tSNE projection stored in the
#'   seuratObj@dr$tsne slot (unless otherwise specified)
#' @export
#'
#' @examples
GoMCtSNE <- function(seuratObj,
                        reduction.use = 'pca',
                        reduction.save = 'tsne',
                        ...){
  seuratObj <- python.dim.reduction.bridge(seuratObj,
                                           reduction.use = reduction.use,
                                           reduction.save = reduction.save,
                                           function.use = multicoreTSNE,
                                           ...)
  return(seuratObj)
}

#' GofitSNE
#'
#' Run tSNE projection on a Seurat object using the FIt-SNE function
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). (default: pca)
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). (default: fitsne)
#' @param ... Extra parameters to pass to the FIt-SNE function.
#'
#' @return Seurat object with the tSNE projection stored in the
#'   seuratObj@dr$tsne slot (unless otherwise specified)
#' @export
#'
#' @examples
GofitSNE <- function(seuratObj,
                     reduction.use = 'pca',
                     reduction.save = 'fitsne',
                     ...){
  seuratObj <- python.dim.reduction.bridge(seuratObj,
                                           reduction.use = reduction.use,
                                           reduction.save = reduction.save,
                                           function.use = fitsne,
                                           ...)
  return(seuratObj)
}

#' GoUMAP
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). (default: pca)
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). (default: umap)
#' @param ... Extra parameters to pass to the umap function.
#'
#' @return Seurat object with the UMAP projection stored in the
#'   seuratObj@dr$umap slot (unless otherwise specified)
#' @export
#'
#' @examples
GoUMAP <- function(seuratObj,
                      reduction.use = 'pca',
                      reduction.save = 'umap',
                      ...){
  seuratObj <- python.dim.reduction.bridge(seuratObj,
                                           reduction.use = reduction.use,
                                           reduction.save = reduction.save,
                                           function.use = umap,
                                           ...)
  return(seuratObj)
}


#' GoPHATE
#'
#' Project trajectory-based dimensional reduction using PHATE
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). (default: pca)
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). (default: phate)
#' @param ... Extra parameters to pass to the phate function.
#'
#' @return Seurat object with the PHATE projection stored in the
#'   seuratObj@dr$phate slot (unless otherwise specified)
#' @export
#'
#' @examples
GoPHATE <- function(seuratObj,
                       reduction.use = 'pca',
                       reduction.save = 'phate',
                       ...){
  seuratObj <- python.dim.reduction.bridge(seuratObj,
                                           reduction.use = reduction.use,
                                           reduction.save = reduction.save,
                                           function.use = phate,
                                           ...)
  return(seuratObj)
}

#' GoPhenoGraph
#'
#' Perform community clustering using Phenograph
#'
#' @param seuratObj
#' @param reduction.use Dimensional reduction to use for clustering calculations
#'   (i.e. pca, ica, cca, etc...)
#' @param k Number of nearest neighbors to use in first step of graph
#'   construction (default = 30) If a list of integers is passed, Phenograph
#'   will be run with each value and the last will be used to set
#'   seuratObj@ident.
#' @param prefix String prefix to used as in the column name entered in the
#'   meta.data slot
#' @param ... Extra parameters to pass to the phenograph function.
#'
#' @importFrom Seurat GetDimReduction SetAllIdent
#'
#' @return A Seurat object with community information stored in
#'   seuratObj@meta.data$prefix# columns, where # = k
#' @export
#'
#' @examples
GoPhenoGraph <- function(seuratObj,
                         reduction.use = 'pca',
                         k = 30,
                         prefix = "community",
                         ...){
  ce <- GetDimReduction(seuratObj,
                        reduction.type = reduction.use,
                        slot = "cell.embeddings")
  for (value in k){
    communities <- phenograph(ce, k = value, ...)
    seuratObj@meta.data[,paste0(prefix,value)] <- communities
  }
  seuratObj <- SetAllIdent(seuratObj, paste0(prefix,value))
  return(seuratObj)
}

#' GoDCA
#'
#' Use deep count autoencoder to impute data
#'
#' @param seuratObj
#' @param ... Extra parameters to pass to the dca function.
#'
#' @importFrom Seurat SetAssayData
#' @importFrom magrittr %>%
#'
#' @return A Seurat object with imputed data stored in seuratObj@dr$dca slot
#' @export
#'
#' @examples
GoDCA <- function(seuratObj, ...){
  exprs <- seuratObj@raw.data %>% as.matrix() %>% t()
  dca_exprs <- dca(exprs)
  seuratObj <- SetAssayData(object = seuratObj,
                            assay.type = "dca",
                            slot = "data",
                            new.data = dca_exprs)
  return(seuratObj)
}
