#' @title PushData
#'
#' @description Helper function to assist entering dimensional reduction data from Python
#' reduction methods
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param python_df Dataframe returned by a Python function
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...)
#' @param assay_used Assay from which the data that is dimensionally reduced comes
#' @param ... Arguments passed to specific downstream methods
#'
#' @importFrom glue glue
#' @importFrom magrittr %<>%
#'
PushData <- function(object, ...) {
  UseMethod("PushData")
}

#' @rdname PushData
#' @method PushData Seurat
#' @importFrom Seurat CreateDimReducObject
#' @return Seurat object
PushData.Seurat <- function(object,
                            python_df,
                            reduction_save,
                            assay_used) {
  python_df %<>% as.matrix()
  rownames(python_df) <- colnames(object)
  reduction_data <- CreateDimReducObject(
    embeddings = python_df,
    assay = assay_used,
    key = as.character(glue("{reduction_save}_"))
  )
  object[[reduction_save]] <- reduction_data
  return(object)
}

#' @rdname PushData
#' @method PushData SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDim<-
#' @return SingleCellExperiment object
PushData.SingleCellExperiment <- function(object,
                                          python_df,
                                          reduction_save) {
  python_df %<>% as.matrix()
  rownames(python_df) <- colnames(object)
  reducedDim(x = object, type = toupper(reduction_save)) <- python_df
  return(object)
}


#' @title ReductionBridge
#'
#' @description Generalized helper function that pulls the data from an object, passes
#' the dataframe to a Python function, and places the resulting dataframe in the
#' appropriate slot
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction_use Prior dimensional reduction to use for calculations
#' (i.e. pca, ica, cca, etc...)
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...)
#' @param function_use Dimensional reduction function to call.
#' @param dims_use Dimensions of `reduction_use` to pass to `function_use`
#' @param ... Extra parameters to pass to the dimensional reduction function.
#'
#' @importFrom glue glue
#' @importFrom rlang %||%
#'
ReductionBridge <- function(object, ...) {
  UseMethod("ReductionBridge")
}


#' @rdname ReductionBridge
#' @method ReductionBridge Seurat
#' @importFrom Seurat Embeddings DefaultAssay
#' @return Seurat object
ReductionBridge.Seurat <- function(object,
                                   reduction_use = "pca",
                                   reduction_save,
                                   function_use,
                                   dims_use = NULL,
                                   ...) {
  if (reduction_use %in% names(object)) {
    cell_embeddings <- Embeddings(object[[reduction_use]])
    adjacencies <- object@graphs
    assay <- DefaultAssay(object = object[[reduction_use]])
    snn <- as.matrix(object@graphs[[glue("{assay}_snn")]])
    snn <- snn[rownames(cell_embeddings),rownames(cell_embeddings)]
  }
  else {
    message(glue("{reduction_use} has not yet been performed"))
    stop()
  }

  dims_use = dims_use %||% 1:ncol(cell_embeddings)

  if (!all(dims_use %in% 1:ncol(cell_embeddings))) {
    stop(glue("You have selected dimensions that are outside the bounds of {reduction_use}"))
  }

  if (match.call()[5] == "fa2()") {
    python_df <- function_use(cell_embeddings[,dims_use], snn, ...)  
  } else {
    python_df <- function_use(cell_embeddings[,dims_use], ...)
  }
  
  object <- PushData(
    object = object,
    python_df = python_df,
    reduction_save = reduction_save,
    assay_used = assay
  )
  return(object)
}

#' @rdname ReductionBridge
#' @method ReductionBridge SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @return SingleCellExperiment object
ReductionBridge.SingleCellExperiment <- function(object,
                                                 reduction_use = "PCA",
                                                 reduction_save,
                                                 function_use,
                                                 dims_use = NULL,
                                                 ...) {
  if (toupper(reduction_use) %in% reducedDimNames(object)) {
    cell_embeddings <- reducedDim(x = object, type = toupper(reduction_use))
  }
  else {
    stop(glue("{reduction_use} has not yet been performed"))
  }

  if (match.call()[5] == "fa2()") {
    stop("Sorry, ForceAtlas2 projection is not currently implemented for SingleCellAssay objects and won't be until
         I can figure out how to get a simultaneous nearest network calculated for it.")
  }

  dims_use = dims_use %||% 1:ncol(cell_embeddings)

  if (!all(dims_use %in% 1:ncol(cell_embeddings))) {
    stop(glue("You have selected dimensions that are outside the bounds of {reduction_use}"))
  }

  python_df <- function_use(cell_embeddings[,dims_use], ...)
  object <- PushData(
    object = object,
    python_df = python_df,
    reduction_save = toupper(reduction_save)
  )
  return(object)
}


#' @title DooptSNE
#'
#' @description Perform tSNE projection on a Seurat object using the
#' Multicore-opt-SNE function
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction_use Prior dimensional reduction to use for calculations
#' (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...). Default: tsne
#' @param dims_use Dimensions from `reduction_use` to pass to PHATE
#' @param ... Extra parameters to pass to the multicoreTSNE function.
#'
#' @export
#'
DooptSNE <- function(object,
                     reduction_use = "pca",
                     reduction_save = "optsne",
                     dims_use = NULL,
                     ...) {
  object <- ReductionBridge(object,
                            reduction_use = reduction_use,
                            reduction_save = reduction_save,
                            function_use = optSNE,
                            dims_use = dims_use,
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
#' @param reduction_use Prior dimensional reduction to use for calculations
#' (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...). Default: openTSNE
#' @param dims_use Dimensions from `reduction_use` to pass to PHATE
#' @param ... Extra parameters to pass to the openTSNE function.
#'
#' @export
#'
DoopenTSNE <- function(object,
                       reduction_use = "pca",
                       reduction_save = "openTSNE",
                       dims_use = NULL,
                       ...) {
  object <- ReductionBridge(object,
                            reduction_use = reduction_use,
                            reduction_save = reduction_save,
                            function_use = openTSNE,
                            dims_use = dims_use,
                            ...
  )
  return(object)
}

#' @title DoUMAP
#'
#' @description Perform UMAP dimentional reduction
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction_use Prior dimensional reduction to use for calculations
#' (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...). Default: umap
#' @param dims_use Dimensions from `reduction_use` to pass to PHATE
#' @param ... Extra parameters to pass to the umap function.
#'
#' @export
#'
DoUMAP <- function(object,
                   reduction_use = "pca",
                   reduction_save = "umap",
                   dims_use = NULL,
                   ...) {
  object <- ReductionBridge(object,
                            reduction_use = reduction_use,
                            reduction_save = reduction_save,
                            function_use = umap,
                            dims_use = dims_use,
                            ...
  )
  return(object)
}


#' @title DoForceAtlas2
#'
#' @description Perform ForceAtlas2 dimentional reduction
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction_use Prior dimensional reduction to use for calculations
#' (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...). Default: umap
#' @param dims_use Dimensions from `reduction_use` to pass to PHATE
#' @param ... Extra parameters to pass to the umap function.
#'
#' @export
#'
DoForceAtlas2 <- function(object,
                          reduction_use = "pca",
                          reduction_save = "fa2",
                          dims_use = NULL,
                          ...) {
  object <- ReductionBridge(object,
                            reduction_use = reduction_use,
                            reduction_save = reduction_save,
                            function_use = fa2,
                            dims_use = dims_use,
                            ...
  )
  return(object)
}


#' @title DoPHATE
#'
#' @description Project trajectory-based dimensional reduction using PHATE
#'
#' @param object A Seurat or SingleCellExperiment object to be transformed.
#' @param reduction_use Prior dimensional reduction to use for calculations
#' (i.e. pca, ica, cca, etc...). Default: pca
#' @param reduction_save Name to use for the reduction (i. e. tsne, umap,
#' etc...). Default: phate
#' @param dims_use Dimensions from `reduction_use` to pass to PHATE
#' @param ... Extra parameters to pass to the phate function.
#'
#' @export
#'
DoPHATE <- function(object,
                    reduction_use = "pca",
                    reduction_save = "phate",
                    dims_use = NULL,
                    ...) {
  object <- ReductionBridge(object,
                            reduction_use = reduction_use,
                            reduction_save = reduction_save,
                            function_use = phate,
                            dims_use = dims_use,
                            ...
  )
  return(object)
}

#' @title DoPhenoGraph
#'
#' @description Perform community clustering using Phenograph
#'
#' @param object A Seurat or SingleCellExperiment object with data to be clustered.
#' @param reduction_use Dimensional reduction to use for clustering calculations
#' (i.e. pca, ica, cca, etc...)
#' @param k Number of nearest neighbors to use in first step of graph
#' construction.  If a list of integers is passed, Phenograph
#' will be Do with each value and the last will be used to set
#' object@ident.  Default = 30
#' @param prefix String prefix to used as in the column name entered in the
#'  meta.data slot
#' @param ... Extra parameters to pass to the phenograph function.
#'
#' @importFrom glue glue
#' @export
#'
DoPhenoGraph <- function(object, ...) {
  UseMethod("DoPhenoGraph")
}

#' @rdname DoPhenoGraph
#' @method DoPhenoGraph Seurat
#' @importFrom Seurat Embeddings AddMetaData Idents<-
#' @return Seurat object
#' @export
DoPhenoGraph.Seurat <- function(object,
                                reduction_use = "pca",
                                k = 30,
                                prefix = "community",
                                ...) {
  if (reduction_use %in% names(object)) {
    cell_embeddings <- Embeddings(object[[reduction_use]])
  }
  else {
    message(glue("{reduction_use} has not yet been performed"))
    stop()
  }

  for (value in k) {
    cluster_name <- glue("{prefix}{value}")
    communities <- phenograph(cell_embeddings, k = value, ...)
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
#' @method DoPhenoGraph SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom SummarizedExperiment colData<- colData
#' @return SingleCellExperiment object
#' @export
DoPhenoGraph.SingleCellExperiment <- function(object,
                                              reduction_use = "pca",
                                              k = 30,
                                              prefix = "community",
                                              ...) {
  if (reduction_use %in% reducedDimNames(object)) {
    cell_embeddings <- reducedDim(x = object, type = reduction_use)
  }
  else {
    message(glue("{reduction_use} has not yet been performed"))
    stop()
  }

  for (value in k) {
    cluster_name <- glue("{prefix}{value}")
    communities <- phenograph(cell_embeddings, k = value, ...)
    colData(object)[[cluster_name]] <- communities
  }

  return(object)
}
