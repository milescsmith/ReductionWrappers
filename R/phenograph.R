#' Phenograph
#'
#' Used to cluster high dimensional data. An R wrapper around the Python
#' Phenograph module found at https://github.com/jacoblevine/PhenoGraph
#'
#' @param r.data.frame data to cluster, or sparse matrix of k-nearest neighbor
#'   graph If ndarray, n-by-d array of n cells in d dimensions If sparse matrix,
#'   n-by-n adjacency matrix
#' @param k Number of nearest neighbors to use in first step of graph
#'   construction (default = 30)
#' @param directed  Whether to use a symmetric (default) or asymmetric
#'   ("directed") graph. The graph construction process produces a directed
#'   graph, which is symmetrized by one of two methods (see below)
#' @param prune Whether to symmetrize by taking the average (prune=FALSE) or
#'   product (prune=TRUE) between the graph and its transpose
#' @param min.cluster.size Cells that end up in a cluster smaller than
#'   min_cluster_size are considered outliers and are assigned to -1 in the
#'   cluster labels
#' @param jaccard If TRUE, use Jaccard metric between k-neighborhoods to build
#'   graph. If FALSE, use a Gaussian kernel.
#' @param primary.metric Distance metric to define nearest neighbors. Options
#'   include: {'euclidean', 'manhattan', 'correlation', 'cosine'} Note that
#'   performance will be slower for correlation and cosine.
#' @param n.jobs Nearest Neighbors and Jaccard coefficients will be computed in
#'   parallel using n_jobs. If n_jobs=NULL, the number of jobs is determined
#'   automatically
#' @param q.tol Tolerance (i.e., precision) for monitoring modularity
#'   optimization
#' @param louvain.time.limit Maximum number of seconds to run modularity
#'   optimization. If exceeded the best result so far is returned
#' @param nn_method Whether to use brute force or kdtree for nearest neighbor
#'   search. For very large high-dimensional data sets, brute force (with
#'   parallel computation) performs faster than kdtree.
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#'
#' @return
#' @export
#'
#' @examples
phenograph <- function(r.data.frame,
                       k = 30,
                       directed = FALSE,
                       prune = FALSE,
                       min.cluster.size = 10,
                       jaccard = TRUE,
                       primary.metric = 'euclidean',
                       n.jobs = NULL,
                       q.tol = 0.001,
                       louvain.time.limit = 2000,
                       nn_method = 'kdtree'){
  if(!py_module_available('phenograph')){
    stop("The MulticoreTSNE module is unavailable.  Please activate the appropriate environment or install the module.")
  }
  if (is.null(n.jobs)){
    n.jobs <- detectCores()
  }
  phenograph.module <- import(module = 'phenograph', delay_load = TRUE)
  phenograph.tuple <- phenograph.module$cluster(r.data.frame,
                                                k = as.integer(k),
                                                directed = directed,
                                                prune = prune,
                                                min_cluster_size = as.integer(min.cluster.size),
                                                jaccard = jaccard,
                                                primary_metric = primary.metric,
                                                n_jobs = as.integer(n.jobs),
                                                q_tol = as.numeric(q.tol),
                                                louvain_time_limit = as.integer(louvain.time.limit),
                                                nn_method = 'kdtree')
  communities <- phenograph.tuple[[1]]
  return(communities)
}
