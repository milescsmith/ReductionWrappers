#' @title PushData
#'
#' @description Helper function to assist entering dimensional reduction data from Python
#' reduction methods
#'
#' @param seuratObj
#' @param python.dataframe Dataframe returned by a Python function
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...)
#' @param assay.used Assay from which the data that is dimensionally reduced comes
#'
#' @importFrom glue glue
#'
#' @return A Seurat object with the dataframe stored in the
#'   seuratObj@dr$reduction.use slot
#' @export
#'
#' @examples
PushData <- function(object, ...) {
  UseMethod("PushData")
}

#' @rdname PushData
#' @method PushData Seurat
#' @importFrom Seurat CreateDimReducObject
#' @return
#' @export
PushData.Seurat <- function(seuratObj,
                            python.dataframe,
                            reduction.save,
                            assay.used) {
  dim.xfer <- as.matrix(python.dataframe)
  rownames(dim.xfer) <- colnames(seuratObj)
  reduction.data <- CreateDimReducObject(
    embeddings = dim.xfer,
    assay = assay.used,
    key = glue("{reduction.save}_")
  )
  seuratObj[[reduction.save]] <- reduction.data
  return(seuratObj)
}

#' @rdname PushData
#' @method PushData seurat
#' @importFrom Seurat SetDimReduction
#' @return
#' @export
PushData.seurat <- function(seuratObj,
                            python.dataframe,
                            reduction.save) {
  dim.xfer <- as.matrix(python.dataframe)
  rownames(dim.xfer) <- rownames(seuratObj@meta.data)
  seuratObj <- SetDimReduction(
    object = seuratObj,
    reduction.type = reduction.save,
    slot = "cell.embeddings",
    new.data = dim.xfer
  )
  seuratObj <- SetDimReduction(
    object = seuratObj,
    reduction.type = reduction.save,
    slot = "key",
    new.data = glue("{reduction.save}_") %>% as.character()
  )
  return(seuratObj)
}

#' @title ReductionBridge
#'
#' @description Generalized helper function that pulls the data from a Seurat object, passes
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
#' @importFrom Seurat Embeddings DefaultAssay
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
ReductionBridge <- function(object, ...) {
  UseMethod("ReductionBridge")
}

#' @rdname ReductionBridge
#' @method ReductionBridge seurat
#' @importFrom Seurat GetDimReduction
#' @return
#' @export
ReductionBridge.seurat <- function(seuratObj,
                                   reduction.use,
                                   reduction.save,
                                   function.use,
                                   ...) {
  cell.embeddings <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "cell.embeddings"
  )
  dim.df <- function.use(cell.embeddings, ...)
  seuratObj <- PushData(
    seuratObj = seuratObj,
    python.dataframe = dim.df,
    reduction.save = reduction.save
  )
  return(seuratObj)
}

#' @rdname ReductionBridge
#' @method ReductionBridge Seurat
#' @importFrom Seurat Embeddings DefaultAssay
#' @return
#' @export
ReductionBridge.Seurat <- function(seuratObj,
                                   reduction.use,
                                   reduction.save,
                                   function.use,
                                   ...) {
  if (reduction.use %in% names(seuratObj)) {
    cell.embeddings <- Embeddings(seuratObj[[reduction.use]])
    assay <- DefaultAssay(object = seuratObj[[reduction.use]])
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }

  dim.df <- function.use(cell.embeddings, ...)
  seuratObj <- PushData(
    seuratObj = seuratObj,
    python.dataframe = dim.df,
    reduction.save = reduction.save,
    assay.used = assay
  )
  return(seuratObj)
}

#' @title DooptSNE
#'
#' @description Perform tSNE projection on a Seurat object using the
#' Multicore-opt-SNE function
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: tsne
#' @param ... Extra parameters to pass to the multicoreTSNE function.
#'
#' @return Seurat object with the tSNE projection stored in the
#'   seuratObj@dr$tsne slot (unless otherwise specified)
#' @export
#'
#' @examples
DooptSNE <- function(seuratObj,
                     reduction.use = "pca",
                     reduction.save = "optsne",
                     ...) {
  seuratObj <- ReductionBridge(seuratObj,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = optSNE,
    ...
  )
  return(seuratObj)
}

#' @title DoopenTSNE
#'
#' @description Perform tSNE projection on a Seurat object using the openTSNE
#' library, with FIt-SNE selected by default
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: openTSNE
#' @param ... Extra parameters to pass to the openTSNE function.
#'
#' @return Seurat object with the tSNE projection stored in the
#'   seuratObj@dr$tsne slot (unless otherwise specified)
#'
#' @export
#'
#' @examples
DoopenTSNE <- function(seuratObj,
                       reduction.use = "pca",
                       reduction.save = "openTSNE",
                       ...) {
  seuratObj <- ReductionBridge(seuratObj,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = openTSNE,
    ...
  )
  return(seuratObj)
}

#' @title DoUMAP
#'
#' @description Perform UMAP dimentional reduction
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: umap
#' @param ... Extra parameters to pass to the umap function.
#'
#' @return Seurat object with the UMAP projection stored in the
#'   seuratObj@dr$umap slot (unless otherwise specified)
#' @export
#'
#' @examples
DoUMAP <- function(seuratObj,
                   reduction.use = "pca",
                   reduction.save = "umap",
                   ...) {
  seuratObj <- ReductionBridge(seuratObj,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = umap,
    ...
  )
  return(seuratObj)
}


#' @title DoPHATE
#'
#' @description Project trajectory-based dimensional reduction using PHATE
#'
#' @param seuratObj
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: phate
#' @param ... Extra parameters to pass to the phate function.
#'
#' @return Seurat object with the PHATE projection stored in the
#'   seuratObj@dr$phate slot (unless otherwise specified)
#' @export
#'
#' @examples
DoPHATE <- function(seuratObj,
                    reduction.use = "pca",
                    reduction.save = "phate",
                    ...) {
  seuratObj <- ReductionBridge(seuratObj,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = phate,
    ...
  )
  return(seuratObj)
}

#' @title DoPhenoGraph
#'
#' @description Perform community clustering using Phenograph
#'
#' @param seuratObj
#' @param reduction.use Dimensional reduction to use for clustering calculations
#'   (i.e. pca, ica, cca, etc...)
#' @param k Number of nearest neighbors to use in first step of graph
#'   construction.  If a list of integers is passed, Phenograph
#'   will be Do with each value and the last will be used to set
#'   seuratObj@ident.  Default = 30
#' @param prefix String prefix to used as in the column name entered in the
#'   meta.data slot
#' @param ... Extra parameters to pass to the phenograph function.
#'
#' @importFrom Seurat AddMetaData Embeddings Idents
#' @importFrom glue glue
#'
#' @return A Seurat object with community information stored in
#'   seuratObj@meta.data$prefix# columns, where # = k
#' @export
#'
#' @examples
DoPhenoGraph <- function(seuratObj,
                         reduction.use = "pca",
                         k = 30,
                         prefix = "community",
                         ...) {
  if (reduction.use %in% names(seuratObj)) {
    cell.embeddings <- Embeddings(seuratObj[[reduction.use]])
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }


  for (value in k) {
    cluster_name <- glue("{prefix}{value}")
    communities <- phenograph(cell.embeddings, k = value, ...)
    seuratObj <- AddMetaData(
      object = seuratObj,
      metadata = communities,
      col.name = cluster_name
    )
    Idents(seuratObj) <- cluster_name
  }

  return(seuratObj)
}
