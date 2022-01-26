#' @title pacmap
#'
#' @description An R wrapper for the PaCMAP Python module found at
#' https://github.com/YingfanWang/PaCMAP
#'
#' PaCMAP (Pairwise Controlled Manifold Approximation) is a
#' dimensionality reduction method that can be used for
#' visualization, preserving both local and global structure
#' of the data in original space. PaCMAP optimizes the low
#' dimensional embedding using three kinds of pairs of points:
#' neighbor pairs (pair_neighbors), mid-near pair (pair_MN),
#' and further pairs (pair_FP).
#'
#' @param rdf A variable by observation data frame
#' @param n_components integer Dimensions of the embedded space. Default: 3
#' @param perplexity numeric The perplexity is related to the
#' number of nearest neighbors that is used in other manifold learning
#' algorithms. Larger datasets usually require a larger perplexity. Consider
#' selecting a value between 5 and 50. The choice is not extremely critical
#' since t-SNE is quite insensitive to this parameter. Default: 30
#' @param early_exaggeration numeric Controls how tight natural
#' clusters in the original space are in the embedded space and how much space
#' will be between them. For larger values, the space between natural clusters
#' will be larger in the embedded space. Again, the choice of this parameter
#' is not very critical. If the cost function increases during initial
#' optimization, the early exaggeration factor or the learning rate might be
#' too high. Default: 12.0
#' @param learning_rate numeric The learning rate for t-SNE is
#' usually in the range [10.0, 1000.0]. If the learning rate is too high, the
#' data may look like a ‘ball’ with any point approximately equidistant from
#' its nearest neighbours. If the learning rate is too low, most points may
#' look compressed in a dense cloud with few outliers. If the cost function
#' gets stuck in a bad local minimum increasing the learning rate may help. Default: 200.0
#' @param n_iter integer Maximum number of iterations for the
#' optimization. Should be at least 250. Default: 1000
#' @param n_iter_without_progress integer Maximum number of
#' iterations without progress before we abort the optimization, used after
#' 250 initial iterations with early exaggeration. Note that progress is only
#' checked every 50 iterations so this value is rounded to the next multiple
#' of 50. Default: 300
#' @param min_grad_norm numeric If the gradient norm is below
#' this threshold, the optimization will be stopped. Default: 1e-7
#' @param metric character or callable The metric to use when calculating distance
#' between instances in a feature array. If metric is a character, it must be one
#' of the options allowed by scipy.spatial.distance.pdist for its metric
#' parameter, or a metric listed in pairwise.PAIRWISE.DISTANCE.FUNCTIONS. If
#' metric is “precomputed”, X is assumed to be a distance matrix.
#' Alternatively, if metric is a callable function, it is called on each pair
#' of instances (rows) and the resulting value recorded. The callable should
#' take two arrays from X as input and return a value indicating the distance
#' between them. The default is “euclidean” which is interpreted as squared
#' euclidean distance.
#' @param init character or numpy array Initialization of
#' embedding. Possible options are ‘random’, ‘pca’, and a numpy array of shape
#' (n.samples, n.components). PCA initialization cannot be used with
#' precomputed distances and is usually more globally stable than random
#' initialization. Default: “random”
#' @param verbose integer Verbosity level. Default: 0
#' @param random_state int, RandomState instance or NULL If int,
#' random.state is the seed used by the random number generator; If
#' RandomState instance, random.state is the random number generator; If NULL,
#' the random number generator is the RandomState instance used by np.random.
#' Note that different initializations might result in different local minima
#' of the cost function. Default: NULL
#' @param method character By default the gradient
#' calculation algorithm uses Barnes-Hut approximation running in \eqn{O(N log N)}
#' time. method=’exact’ will run on the slower, but exact, algorithm in \eqn{O(N^2)}
#' time. The exact algorithm should be used when nearest-neighbor errors need
#' to be better than 3%. However, the exact method cannot scale to millions of
#' examples. Default: ‘barnes.hut’
#' @param angle numeric Only used if method=’barnes.hut’ This is
#' the trade-off between speed and accuracy for Barnes-Hut T-SNE. ‘angle’ is
#' the angular size (also referred to as theta) of a distant node as
#' measured from a point. If this size is below ‘angle’ then it is used as a
#' summary node of all points contained within it. This method is not very
#' sensitive to changes in this parameter in the range of 0.2 - 0.8. Angle
#' less than 0.2 has quickly increasing computation time and angle greater 0.8
#' has quickly increasing error.#' Default: 0.5
#' @param auto_iter boolean Should optimal parameters be determined?
#' If false, behaves like stock MulticoreTSNE Default: TRUE
#' @param auto_iter_end intNumber of iterations for parameter
#' optimization. Default: 5000
#' @param n_jobs Number of processors to use.  Default: all.
#'
#' @importFrom reticulate import py_module_available
#' @importFrom parallel detectCores
#'
#' @return data.frame with tSNE coordinates
#' @export
#'
pacmap <- function(rdf,
                   n_dims         = 2,
                   n_neighbors    = NULL,
                   MN_ratio       = 0.5,
                   FP_ratio       = 2.0,
                   pair_neighbors = NULL,
                   pair_MN        = NULL,
                   pair_FP        = NULL,
                   distance       = "euclidean",
                   lr             = 1.0,
                   num_iters      = 450,
                   verbose        = FALSE,
                   apply_pca      = TRUE,
                   intermediate   = FALSE){
                     
  if (!reticulate::py_module_available("pacmap")){
    stop("The pacmap module is unavailable.  Please activate the appropriate environment or install the module.")
  }

  pacmap_module <-
    reticulate::import(
      module = "pacmap",
      delay_load = TRUE
    )

  pacmap <-
    pacmap_module$PaCMAP(
      n_dims         = as.integer(n_dims),
      n_neighbors    = n_neighbors,
      MN_ratio       = as.numeric(MN_ratio),
      FP_ratio       = as.numeric(FP_ratio),
      pair_neighbors = pair_neighbors,
      pair_MN        = pair_MN,
      pair_FP        = pair_FP,
      distance       = distance,
      lr             = as.numeric(lr),
      num_iters      = as.integer(num_iters),
      verbose        = verbose,
      apply_pca      = apply_pca,
      intermediate   = intermediate
    )

  pacmap$fit_transform(rdf)
}
