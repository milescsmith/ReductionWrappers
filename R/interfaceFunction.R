#' @title PushData
#'
#' @description Helper function to assist entering dimensional reduction data from Python
#' reduction methods
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param python_df Dataframe returned by a Python function
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...)
#' @param assay.used Assay from which the data that is dimensionally reduced comes
#' @param ... Arguments passed to specific downstream methods
#'
#' @importFrom glue glue
#' @importFrom magrittr %<>%
#'
#' @return
#' @examples
PushData <- function(object, ...) {
  UseMethod("PushData")
}

#' @rdname PushData
#' @method PushData Seurat
#' @import Seurat
#' @return
PushData.Seurat <- function(object,
                            python_df,
                            reduction.save,
                            assay.used) {
  python_df %<>% as.matrix()
  rownames(python_df) <- colnames(object)
  reduction.data <- CreateDimReducObject(
    embeddings = python_df,
    assay = assay.used,
    key = as.character(glue("{reduction.save}_"))
  )
  object[[reduction.save]] <- reduction.data
  return(object)
}

#' @rdname PushData
#' @method PushData seurat
#' @import Seurat
#' @return
PushData.seurat <- function(object,
                            python_df,
                            reduction.save) {
  python_df %<>% as.matrix()
  rownames(python_df) <- rownames(object@meta.data)
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.save,
    slot = "cell.embeddings",
    new.data = python_df
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.save,
    slot = "key",
    new.data = as.character(glue("{reduction.save}_"))
  )
  return(object)
}

#' @rdname PushData
#' @method PushData SingleCellExperiment
#' @import SingleCellExperiment
#' @return
PushData.SingleCellExperiment <- function(object,
                                          python_df,
                                          reduction.save) {
  python_df %<>% as.matrix()
  rownames(python_df) <- colnames(object)
  reducedDim(x = object, type = toupper(reduction.save)) <- python_df
  return(object)
}

#' @title ReductionBridge
#'
#' @description Generalized helper function that pulls the data from an object, passes
#' the dataframe to a Python function, and places the resulting dataframe in the
#' appropriate slot
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...)
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...)
#' @param function.use Dimensional reduction function to call.
#' @param ... Extra parameters to pass to the dimensional reduction function.
#'
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
#' @import Seurat
#' @return
ReductionBridge.seurat <- function(object,
                                   reduction.use = "pca",
                                   reduction.save,
                                   function.use,
                                   ...) {
  cell.embeddings <- GetCellEmbeddings(
    object = object,
    reduction.type = reduction.use,
  )
  python_df <- function.use(cell.embeddings, ...)
  object <- PushData(
    object = object,
    python_df = python_df,
    reduction.save = reduction.save
  )
  return(object)
}

#' @rdname ReductionBridge
#' @method ReductionBridge Seurat
#' @import Seurat
#' @return
ReductionBridge.Seurat <- function(object,
                                   reduction.use = "pca",
                                   reduction.save,
                                   function.use,
                                   ...) {
  if (reduction.use %in% names(object)) {
    cell.embeddings <- Embeddings(object[[reduction.use]])
    assay <- DefaultAssay(object = object[[reduction.use]])
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }

  python_df <- function.use(cell.embeddings, ...)
  object <- PushData(
    object = object,
    python_df = python_df,
    reduction.save = reduction.save,
    assay.used = assay
  )
  return(object)
}

#' @rdname ReductionBridge
#' @method ReductionBridge SingleCellExperiment
#' @import SingleCellExperiment
#' @return
ReductionBridge.SingleCellExperiment <- function(object,
                                                 reduction.use = "PCA",
                                                 reduction.save,
                                                 function.use,
                                                 ...) {
  if (toupper(reduction.use) %in% reducedDimNames(object)) {
    cell.embeddings <- reducedDim(x = object, type = toupper(reduction.use))
  }
  else {
    stop(glue("{reduction.use} has not yet been performed"))
  }

  python_df <- function.use(cell.embeddings, ...)
  object <- PushData(
    object = object,
    python_df = python_df,
    reduction.save = toupper(reduction.save)
  )
  return(object)
}

#' @title DooptSNE
#'
#' @description Perform tSNE projection on a Seurat object using the
#' Multicore-opt-SNE function
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: tsne
#' @param ... Extra parameters to pass to the multicoreTSNE function.
#'
#' @return
#' @export
#'
#' @examples
DooptSNE <- function(object,
                     reduction.use = "pca",
                     reduction.save = "optsne",
                     ...) {
  object <- ReductionBridge(object,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = optSNE,
    ...
  )
  return(object)
}

#' @title DoopenTSNE
#'
#' @description Perform tSNE projection on a Seurat object using the openTSNE
#' library, with FIt-SNE selected by default
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: openTSNE
#' @param ... Extra parameters to pass to the openTSNE function.
#'
#' @return
#' @export
#'
#' @examples
DoopenTSNE <- function(object,
                       reduction.use = "pca",
                       reduction.save = "openTSNE",
                       ...) {
  object <- ReductionBridge(object,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = openTSNE,
    ...
  )
  return(object)
}

#' @title DoUMAP
#'
#' @description Perform UMAP dimentional reduction
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: umap
#' @param ... Extra parameters to pass to the umap function.
#'
#' @return
#' @export
#'
#' @examples
DoUMAP <- function(object,
                   reduction.use = "pca",
                   reduction.save = "umap",
                   ...) {
  object <- ReductionBridge(object,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = umap,
    ...
  )
  return(object)
}


#' @title DoPHATE
#'
#' @description Project trajectory-based dimensional reduction using PHATE
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction.use Prior dimensional reduction to use for calculations
#'   (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction.save Name to use for the reduction (i. e. tsne, umap,
#'   etc...). Default: phate
#' @param ... Extra parameters to pass to the phate function.
#'
#' @return
#' @export
#'
#' @examples
DoPHATE <- function(object,
                    reduction.use = "pca",
                    reduction.save = "phate",
                    ...) {
  object <- ReductionBridge(object,
    reduction.use = reduction.use,
    reduction.save = reduction.save,
    function.use = phate,
    ...
  )
  return(object)
}

#' @title DoPhenoGraph
#'
#' @description Perform community clustering using Phenograph
#'
#' @param object A Seurat or SingleCellExperiment object with data to be clustered.
#' @param reduction.use Dimensional reduction to use for clustering calculations
#'   (i.e. pca, ica, cca, etc...)
#' @param k Number of nearest neighbors to use in first step of graph
#'   construction.  If a list of integers is passed, Phenograph
#'   will be Do with each value and the last will be used to set
#'   object@ident.  Default = 30
#' @param prefix String prefix to used as in the column name entered in the
#'   meta.data slot
#' @param ... Extra parameters to pass to the phenograph function.
#'
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
DoPhenoGraph <- function(object, ...) {
  UseMethod("DoPhenoGraph")
}

#' @rdname DoPhenoGraph
#' @method DoPhenoGraph Seurat
#' @import Seurat
#' @return
#' @export
DoPhenoGraph.Seurat <- function(object,
                                reduction.use = "pca",
                                k = 30,
                                prefix = "community",
                                ...) {
  if (reduction.use %in% names(object)) {
    cell.embeddings <- Embeddings(object[[reduction.use]])
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }

  for (value in k) {
    cluster_name <- glue("{prefix}{value}")
    communities <- phenograph(cell.embeddings, k = value, ...)
    object <- AddMetaData(
      object = object,
      metadata = communities,
      col.name = cluster_name
    )
    Idents(object) <- cluster_name
  }

  return(object)
}

#' @rdname DoPhenoGraph
#' @method DoPhenoGraph seurat
#' @import Seurat
#' @return
#' @export
DoPhenoGraph.seurat <- function(object,
                                reduction.use = "pca",
                                k = 30,
                                prefix = "community",
                                ...) {
  if (reduction.use %in% names(object@dr)) {
    cell.embeddings <- GetCellEmbeddings(
      object = object,
      reduction.type = reduction.use,
    )
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }

  for (value in k) {
    cluster_name <- glue("{prefix}{value}")
    communities <- phenograph(cell.embeddings, k = value, ...)
    object <- AddMetaData(
      object = object,
      metadata = communities,
      col.name = cluster_name
    )
    object <- SetAllIdent(object, cluster_name)
  }

  return(object)
}

#' @rdname DoPhenoGraph
#' @method DoPhenoGraph SingleCellExperiment
#' @import SingleCellExperiment
#' @return
#' @export
DoPhenoGraph.SingleCellExperiment <- function(object,
                                              reduction.use = "pca",
                                              k = 30,
                                              prefix = "community",
                                              ...) {
  if (reduction.use %in% reducedDimNames(object)) {
    cell.embeddings <- reducedDim(x = object, type = reduction.use)
  }
  else {
    message(glue("{reduction.use} has not yet been performed"))
    stop()
  }

  for (value in k) {
    cluster_name <- glue("{prefix}{value}")
    communities <- phenograph(cell.embeddings, k = value, ...)
    colData(object)[[cluster_name]] <- communities
  }

  return(object)
}
