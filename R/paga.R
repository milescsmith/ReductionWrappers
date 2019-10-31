#' @title PAGA
#'
#' @description Use scanpy.tl.paga() to produce a partition-based graph abstraction for a Seurat
#' object and use that to initialize a UMAP.  Additionally, runs cluster determination via the
#' 'leiden' or 'louvain' algorithms.
#'
#' If dimensional reduction has already been performed (PCA, ICA, or harmony), that is used to
#' find neighbors, otherwise PCA is run.
#'
#' Parameters are prefixed by the step to which they correspond (i.e. "neighbors_" are
#' passed to scanpy.pp.neighbors())
#'
#' Heavily based on the fantastic walk through found at https://romanhaa.github.io/blog/paga_to_r/
#'
#' @param seurat_obj
#' @param assay Seurat object assay to use when converting to Scanpy object
#' @param slim Temporarily discard all unnecessary data from the Seurat object (i.e. keep only the normalized data for the assay and reduction used for neighborhood calculation).  May help when performing PAGA on large objects. (Default: FALSE)
#' @param seurat_grouping Force PAGA to use this metadata grouping variable. (Default: NULL)
#' @param set_ident Set the cluster identity for each cell when returning the object? (Default: TRUE)
#' @param clustering_algorithm Whether to use the "louvain" or "leiden" algorithms (Default: "leiden")
#' @param clustering_resolution Resolution to pass to the clustering algorith (Default: 1.0)
#' @param reduction_name dimensional reduction name, `umap` by default
#' @param reduction_key dimensional reduction key, specifies the string before the number for the dimension names. `umap` by default
#' @param edge_filter_weight Set edges with a weight below this threshold to NA (Default: 0.1)
#' @param neighbors_n_neighbors
#' @param neighbors_n_pcs
#' @param neighbors_use_rep
#' @param neighbors_knn
#' @param neighbors_random_state
#' @param neighbors_method
#' @param neighbors_metric
#' @param clustering_restrict_to
#' @param clustering_random_state
#' @param clustering_key_added
#' @param clustering_adjacency
#' @param clustering_directed
#' @param clustering_use_weights
#' @param clustering_n_iterations
#' @param clustering_partition_type
#' @param paga_show
#' @param paga_threshold
#' @param paga_layout
#' @param paga_init_pos
#' @param paga_root
#' @param paga_single_component
#' @param paga_random_state
#' @param umap_min_dist
#' @param umap_spread
#' @param umap_n_components
#' @param umap_alpha
#' @param umap_gamma
#' @param umap_negative_sample_rate
#' @param umap_init_pos
#'
#' @return
#' @export
#'
#' @importFrom s2a convert_to_anndata
#' @importFrom glue glue
#' @importFrom reticulate import
#' @importFrom Seurat DietSeurat Idents<-
#'
#' @examples
PAGA <- function(seurat_obj,
                 assay = "RNA",
                 slim = FALSE,
                 seurat_grouping = NULL,
                 set_ident = TRUE,
                 clustering_algorithm = "leiden",
                 reduction_name = "umap",
                 reduction_key = "umap_",
                 edge_filter_weight = 0.10,

                 neighbors_n_neighbors = 15,
                 neighbors_n_pcs = NULL,
                 neighbors_use_rep = "pca",
                 neighbors_knn = TRUE,
                 neighbors_random_state = 0,
                 neighbors_method = 'umap',
                 neighbors_metric = 'euclidean',

                 clustering_resolution = 1.0,
                 clustering_restrict_to=NULL,
                 clustering_random_state=0,
                 clustering_key_added=NULL,
                 clustering_adjacency=NULL,
                 clustering_directed=TRUE,
                 clustering_use_weights=TRUE,
                 clustering_n_iterations=-1,
                 clustering_partition_type=NULL,

                 paga_show = FALSE,
                 paga_plot = FALSE,
                 paga_add_pos = TRUE,
                 paga_threshold=0.01,
                 paga_layout=NULL,
                 paga_init_pos=NULL,
                 paga_root=0.0,
                 paga_single_component=NULL,
                 paga_random_state=0.0,

                 umap_min_dist=0.5,
                 umap_spread=1.0,
                 umap_n_components=2,
                 umap_alpha=1.0,
                 umap_gamma=1.0,
                 umap_negative_sample_rate=5,
                 umap_init_pos='spectral'){

  if (isTRUE(slim)){
    slimmed_obj <- DietSeurat(object = seurat_obj,
                              assay = assay,
                              dimreducs = neighbors_use_rep,
                              counts = FALSE)
    alpha <- convert_to_anndata(object = slimmed_obj,
                                assay = assay)
  } else {
    alpha <- convert_to_anndata(object = seurat_obj,
                                assay = assay)
  }

  sc <- import("scanpy",
               delay_load = TRUE)

  if (glue("X_{neighbors_use_rep}") %in% alpha$obsm_keys()){
    sc$pp$neighbors(adata = alpha,
                    n_neighbors = as.integer(neighbors_n_neighbors),
                    n_pcs = neighbors_n_pcs,
                    use_rep = glue("X_{neighbors_use_rep}"),
                    knn = neighbors_knn,
                    random_state = as.integer(neighbors_random_state),
                    method = neighbors_method,
                    metric = neighbors_metric)
  } else {
    if (length(alpha$obsm_keys()) > 0) {
      message(glue("{neighbors_use_rep} was not found.  Performing PCA..."))
    } else {
      message("No reduced dimensional reductions found.  Performing PCA...")
    }
    sc$tl$pca(alpha)
  }

  clustering_key_added <- clustering_key_added %||% clustering_algorithm

  if (!clustering_algorithm %in% c("leiden", "louvain")){
    stop("Unknown clustering algorithm specified.")
  }

  if (is.null(seurat_grouping)){
    grouping <- clustering_algorithm
    sc$tl[[grouping]](adata = alpha,
                                  resolution = as.numeric(clustering_resolution),
                                  restrict_to = clustering_restrict_to,
                                  random_state = as.integer(clustering_random_state),
                                  key_added = clustering_key_added,
                                  adjacency = clustering_adjacency,
                                  directed = clustering_directed,
                                  use_weights = clustering_use_weights,
                                  n_iterations = as.integer(clustering_n_iterations),
                                  partition_type = clustering_partition_type)
    sc$tl$paga(adata = alpha,
               groups = clustering_key_added)
    seurat_obj@meta.data[[grouping]] <- alpha$obs[[grouping]]
  } else {
    grouping <- seurat_grouping
    sc$tl$paga(adata = alpha,
               groups = grouping)
  }

  sc$pl$paga(adata = alpha,
             show = paga_show,
             threshold=as.numeric(paga_threshold),
             layout=paga_layout,
             init_pos=paga_init_pos,
             root=paga_root,
             single_component=paga_single_component,
             random_state=as.integer(paga_random_state,
             plot=FALSE,
             add_pos=TRUE))

  sc$tl$umap(adata = alpha,
             init_pos = "paga",
             min_dist=as.numeric(umap_min_dist),
             spread=as.numeric(umap_spread),
             n_components=as.integer(umap_n_components),
             alpha=as.numeric(umap_alpha),
             gamma=as.numeric(umap_gamma),
             negative_sample_rate=as.integer(umap_negative_sample_rate))

  paga <- list(
    connectivities = alpha$uns$paga$connectivities$todense() %>%
      `rownames<-`(levels(alpha$obs[[alpha$uns$paga$groups]])) %>%
      `colnames<-`(levels(alpha$obs[[alpha$uns$paga$groups]])),
    connectivities_tree = alpha$uns$paga$connectivities$todense(),
    group_name = alpha$uns$paga$groups,
    groups = levels(alpha$obs[[alpha$uns$paga$groups]]),
    group_colors = setNames(alpha$uns[[glue("{alpha$uns$paga$groups}_colors")]],
                            0:(nrow(alpha$uns$paga$pos)-1)),
    position = as_tibble(
      cbind(
        levels(alpha$obs[[alpha$uns$paga$groups]]),
        alpha$uns$paga$pos),
      .name_repair = ~make.names(c("group","x", "y"))
    ),
    umap = as_tibble(alpha$obsm['X_umap'],
                     .name_repair = ~make.names(names = paste0("UMAP_",
                                                               1:ncol(alpha$obsm['X_umap'])),
                                                unique = TRUE))
  )

  paga$edges <- tibble(
    group1 = paga$groups[row(paga$connectivities)[upper.tri(paga$connectivities)]],
    group2 = paga$groups[col(paga$connectivities)[upper.tri(paga$connectivities)]],
    weight = paga$connectivities[upper.tri(paga$connectivities)]
  ) %>%
    mutate(
      x1 = paga$position$x[match(.$group1, rownames(paga$position))],
      y1 = paga$position$y[match(.$group1, rownames(paga$position))],
      x2 = paga$position$x[match(.$group2, rownames(paga$position))],
      y2 = paga$position$y[match(.$group2, rownames(paga$position))]
    ) %>%
    filter(weight >= edge_filter_weight)

  paga_umap <- CreateDimReducObject(embeddings = alpha$obsm[['X_umap']] %>%
                                      `rownames<-`(colnames(seurat_obj[[assay]])) %>%
                                      `colnames<-`(paste0("UMAP_",
                                                          1:ncol(alpha$obsm['X_umap']))),
                                    assay = "RNA",
                                    key = reduction_key)

  seurat_obj[[reduction_name]] <- paga_umap

  seurat_obj@misc$paga <- paga

  if (isTRUE(set_ident)){
    Idents(seurat_obj) <- seurat_obj@meta.data[[grouping]]
  }

  return(seurat_obj)
}


#' @title PAGAplot
#'
#' @description Plot the results from PAGA
#'
#' @param seurat_obj
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_color_manual geom_text labs
#'
#' @return
#' @export
#'
#' @examples
PAGAplot <- function(seurat_obj){
  seurat_obj@misc$paga$position %>%
    ggplot(aes(x, y)) +
    geom_segment(
      data = seurat_obj@misc$paga$edges,
      aes(x = x1,
          y = y1,
          xend = x2,
          yend = y2,
          size = weight*3),
      colour = "black",
      show.legend = FALSE
    ) +
    geom_point(
      aes(color = group),
      size = 7,
      alpha = 1,
      show.legend = FALSE) +
    scale_color_manual(values = seurat_obj@misc$paga$group_colors) +
    geom_text(aes(label = group),
              color = "black",
              fontface = "bold") +
    labs(x = "UMAP_1",
         y = "UMAP_2") +
    theme_cowplot()

}
