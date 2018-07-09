#' multicoreTSNE
#'
#' An R wrapper for the Multicore-TSNE Python library found at https://github.com/DmitryUlyanov/Multicore-TSNE
#'
#' @param r.data.frame
#' @param n.components
#' @param perplexity
#' @param early_exaggeration
#' @param learning_rate
#' @param n.iter
#' @param n_iter_without_progress
#' @param min_grad_norm
#' @param metric
#' @param init
#' @param verbose
#' @param random_state
#' @param method
#' @param angle
#' @param n.jobs
#' @param cheat_metric
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#'
#' @return dataframe
#' @export
#'
#' @examples
multicoreTSNE <- function(r.data.frame,
                          n.components = 3,
                          perplexity = 30,
                          early.exaggeration = 12,
                          learning.rate = 200,
                          n.iter = 2000,
                          n.iter.without.progress = 30,
                          min.grad.norm = 1e-07,
                          metric = 'euclidean',
                          init = 'random',
                          verbose = 1,
                          random.state = NULL,
                          method = 'barnes_hut',
                          angle = 0.5,
                          n.jobs = NULL,
                          cheat.metric = TRUE
                          ){
  mctsne.module <- import(module = 'MulticoreTSNE', delay_load = TRUE)
  if (is.null(n.jobs)){
    n.jobs <- detectCores()
  }
  mctsne <- mctsne.module$MulticoreTSNE(n_components = as.integer(n.components),
                                        perplexity = as.double(perplexity),
                                        early_exaggeration = as.integer(early.exaggeration),
                                        learning_rate = as.integer(learning.rate),
                                        n_iter = as.integer(n.iter),
                                        n_iter_without_progress = as.integer(n.iter.without.progress),
                                        min_grad_norm = as.integer(min.grad.norm),
                                        metric = metrix,
                                        init = init,
                                        verbose = as.integer(verbose),
                                        random_state = rrandom.state,
                                        method = method,
                                        angle = as.double(angle),
                                        n_jobs = n.jobs,
                                        cheat_metric = cheat.metric)
  mctsne.df <- mctsne$fit_transform(r.data.frame)
  return(mctsne.df)
}

#' umap
#'
#' An R wrapper around the UMAP Python library found at https://github.com/lmcinnes/umap
#'
#' Uniform Manifold Approximation and Projection
#'
#'  Finds a low dimensional embedding of the data that approximates
#'  an underlying manifold.
#'
#' @param r.data.frame
#' @param n.neighbors: double (optional, default 30)
#'      The size of local neighborhood (in terms of number of neighboring
#'      sample points) used for manifold approximation. Larger values
#'      result in more global views of the manifold, while smaller
#'      values result in more local data being preserved. In general
#'      values should be in the range 2 to 100.
#' @param n.components: int (optional, default 3)
#'      The dimension of the space to embed into. This defaults to 3 to
#'      provide easy visualization, but can reasonably be set to any
#'      integer value in the range 2 to 100.
#' @param metric: string or function (optional, default 'euclidean')
#'      The metric to use to compute distances in high dimensional space.
#'      If a string is passed it must match a valid predefined metric. If
#'      a general metric is required a function that takes two 1d arrays and
#'      returns a float can be provided. For performance purposes it is
#'      required that this be a numba jit'd function. Valid string metrics
#'      include:
#'          * euclidean
#'          * manhattan
#'          * chebyshev
#'          * minkowski
#'          * canberra
#'          * braycurtis
#'          * mahalanobis
#'          * wminkowski
#'          * seuclidean
#'          * cosine
#'          * correlation
#'          * haversine
#'          * hamming
#'          * jaccard
#'          * dice
#'          * russelrao
#'          * kulsinski
#'          * rogerstanimoto
#'          * sokalmichener
#'          * sokalsneath
#'          * yule
#'      Metrics that take arguments (such as minkowski, mahalanobis etc.)
#'      can have arguments passed via the metric_kwds dictionary. At this
#'      time care must be taken and dictionary elements must be ordered
#'      appropriately; this will hopefully be fixed in the future.
#' @param negative.sample.rate: int (optional, default 5)
#'      The number of negative edge/1-simplex samples to use per positive
#'      edge/1-simplex sample in optimizing the low dimensional embedding.
#' @param alpha: float (optional, default 1.0)
#'      The initial learning rate for the embedding optimization.
#' @param init: string (optional, default 'spectral')
#'      How to initialize the low dimensional embedding. Options are:
#'          * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
#'          * 'random': assign initial embedding positions at random.
#'          * A numpy array of initial embedding positions.
#' @param min.dist: float (optional, default 0.1)
#'      The effective minimum distance between embedded points. Smaller values
#'      will result in a more clustered/clumped embedding where nearby points
#'      on the manifold are drawn closer together, while larger values will
#'      result on a more even dispersal of points. The value should be set
#'      elative to the ``spread`` value, which determines the scale at which
#'      embedded points will be spread out.
#' @param spread: float (optional, default 1.0)
#'      The effective scale of embedded points. In combination with ``min_dist``
#'      this determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio: float (optional, default 1.0)
#'      Interpolate between (fuzzy) union and intersection as the set operation
#'      used to combine local fuzzy simplicial sets to obtain a global fuzzy
#'      simplicial sets. Both fuzzy set operations use the product t-norm.
#'      The value of this parameter should be between 0.0 and 1.0; a value of
#'      1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy
#'      intersection.
#' @param local.connectivity: int (optional, default 1)
#'      The local connectivity required -- i.e. the number of nearest
#'      neighbors that should be assumed to be connected at a local level.
#'      The higher this value the more connected the manifold becomes
#'      locally. In practice this should be not more than the local intrinsic
#'      dimension of the manifold.
#' @param gamma: float (optional, default 1.0)
#'      Weighting applied to negative samples in low dimensional embedding
#'      optimization. Values higher than one will result in greater weight
#'      being given to negative samples.
#' @param bandwidth: float (optional, default 1.0)
#'      The effective bandwidth of the kernel if we view the algorithm as
#'      similar to Laplacian eigenmaps. Larger values induce more
#'      connectivity and a more global view of the data, smaller values
#'      concentrate more locally.
#' @param random.state: int, RandomState instance or None, optional (default: None)
#'      If int, random_state is the seed used by the random number generator;
#'      If RandomState instance, random_state is the random number generator;
#'      If None, the random number generator is the RandomState instance used
#'      by `np.random`.
#' @param angular.rp.forest: bool (optional, default False)
#'      Whether to use an angular random projection forest to initialise
#'      the approximate nearest neighbor search. This can be faster, but is
#'      mostly on useful for metric that use an angular style distance such
#'      as cosine, correlation etc. In the case of those metrics angular forests
#'      will be chosen automatically.
#' @param verbose: bool (optional, default False)
#'      Controls verbosity of logging.
#'
#' @importFrom reticulate import
#'
#' @return
#' @export
#'
#' @examples
umap <- function(r.data.frame,
                 n.neighbors = 30,
                 min.dist = 0.1,
                 n.components = 3,
                 metric = 'euclidean',
                 init = 'spectral',
                 alpha = 1.0,
                 spread = 1.0,
                 bandwidth = 1.0,
                 random.state = NULL,
                 angular.rp.forest = FALSE,
                 set.op.mix.ratio = 1.0,
                 gamma = 1.0,
                 negative.sample.rate = 5,
                 verbose = TRUE){
  umap.module <- import(module = 'umap', delay_load = TRUE)
  umap.embed <- umap.module$UMAP(n_neighbors = as.integer(n.neighbors),
                                 min_dist = as.double(min.dist),
                                 n_components = as.integer(n.components),
                                 metric = metric,
                                 init = init,
                                 alpha = as.double(alpha),
                                 spread = as.double(spread),
                                 bandwidth = as.double(bandwidth),
                                 random_state = random.state,
                                 angular_rp_forest = angular.rp.forest,
                                 verbose = verbose,
                                 set_op_mix_ratio = as.double(set.op.mix.ratio),
                                 gamma = as.double(gamma),
                                 negative_sample_rate = negative.sample.rate)
  umap.df <- umap.embed$fit_transform(r.data.frame)
  return(umap.df)
}

#' Phenograph
#'
#' Used to cluster high dimensional data.
#' An R wrapper around the Python Phenograph library found at https://github.com/jacoblevine/PhenoGraph
#'
#' @param r.data.frame data to cluster, or sparse matrix of k-nearest neighbor graph
#' If ndarray, n-by-d array of n cells in d dimensions
#' If sparse matrix, n-by-n adjacency matrix
#' @param k Number of nearest neighbors to use in first step of graph construction
#' @param directed  Whether to use a symmetric (default) or asymmetric ("directed") graph.
#' The graph construction process produces a directed graph, which is symmetrized by one of two methods (see below)
#' @param prune Whether to symmetrize by taking the average (prune=False) or product (prune=True) between the graph
# and its transpose
#' @param min.cluster.size Cells that end up in a cluster smaller than min_cluster_size are considered outliers and are assigned to -1 in the cluster labels
#' @param jaccard If True, use Jaccard metric between k-neighborhoods to build graph. If False, use a Gaussian kernel.
#' @param primary.metric Distance metric to define nearest neighbors.
#' Options include: {'euclidean', 'manhattan', 'correlation', 'cosine'}
#' Note that performance will be slower for correlation and cosine.
#' @param n.jobs Nearest Neighbors and Jaccard coefficients will be computed in parallel using n_jobs. If n_jobs=NULL, the number of jobs is determined automatically
#' @param q.tol Tolerance (i.e., precision) for monitoring modularity optimization
#' @param louvain.time.limit Maximum number of seconds to run modularity optimization. If exceeded the best result so far is returned
#' @param nn_method Whether to use brute force or kdtree for nearest neighbor search. For very large high-dimensional data sets, brute force (with parallel computation) performs faster than kdtree.
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
                       jaccard = True,
                       primary.metric = 'euclidean',
                       n.jobs = NULL,
                       q.tol = 0.001,
                       louvain.time.limit = 2000,
                       nn_method = 'kdtree'){
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
                                                q_tol = as.double(q.tol),
                                                louvain_time_limit = as.integer(louvain.time.limit),
                                                nn_method = 'kdtree')
  communities <- phenograph.tuple[[3]]
  return(communities)
}

