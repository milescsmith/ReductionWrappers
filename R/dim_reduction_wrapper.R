#' multicoreTSNE
#'
#' An R wrapper for the Multicore-TSNE Python module found at
#' https://github.com/DmitryUlyanov/Multicore-TSNE
#'
#' @param r.data.frame
#' @param n.components integer (default: 3) Dimension of the embedded space.
#' @param perplexity numeric (default: 30) The perplexity is related to the
#'   number of nearest neighbors that is used in other manifold learning
#'   algorithms. Larger datasets usually require a larger perplexity. Consider
#'   selecting a value between 5 and 50. The choice is not extremely critical
#'   since t-SNE is quite insensitive to this parameter.
#' @param early.exaggeration numeric (default: 12.0) Controls how tight natural
#'   clusters in the original space are in the embedded space and how much space
#'   will be between them. For larger values, the space between natural clusters
#'   will be larger in the embedded space. Again, the choice of this parameter
#'   is not very critical. If the cost function increases during initial
#'   optimization, the early exaggeration factor or the learning rate might be
#'   too high.
#' @param learning.rate numeric (default: 200.0) The learning rate for t-SNE is
#'   usually in the range [10.0, 1000.0]. If the learning rate is too high, the
#'   data may look like a ‘ball’ with any point approximately equidistant from
#'   its nearest neighbours. If the learning rate is too low, most points may
#'   look compressed in a dense cloud with few outliers. If the cost function
#'   gets stuck in a bad local minimum increasing the learning rate may help.
#' @param n.iter integer (default: 1000) Maximum number of iterations for the
#'   optimization. Should be at least 250.
#' @param n.iter.without.progress integer (default: 300) Maximum number of
#'   iterations without progress before we abort the optimization, used after
#'   250 initial iterations with early exaggeration. Note that progress is only
#'   checked every 50 iterations so this value is rounded to the next multiple
#'   of 50.
#' @param min.grad.norm numeric (default: 1e-7) If the gradient norm is below
#'   this threshold, the optimization will be stopped.
#' @param metric character or callable The metric to use when calculating distance
#'   between instances in a feature array. If metric is a character, it must be one
#'   of the options allowed by scipy.spatial.distance.pdist for its metric
#'   parameter, or a metric listed in pairwise.PAIRWISE.DISTANCE.FUNCTIONS. If
#'   metric is “precomputed”, X is assumed to be a distance matrix.
#'   Alternatively, if metric is a callable function, it is called on each pair
#'   of instances (rows) and the resulting value recorded. The callable should
#'   take two arrays from X as input and return a value indicating the distance
#'   between them. The default is “euclidean” which is interpreted as squared
#'   euclidean distance.
#' @param init character or numpy array (default: “random”) Initialization of
#'   embedding. Possible options are ‘random’, ‘pca’, and a numpy array of shape
#'   (n.samples, n.components). PCA initialization cannot be used with
#'   precomputed distances and is usually more globally stable than random
#'   initialization.
#' @param verbose integer (default: 0) Verbosity level.
#' @param random.state int, RandomState instance or NULL (default: NULL) If int,
#'   random.state is the seed used by the random number generator; If
#'   RandomState instance, random.state is the random number generator; If NULL,
#'   the random number generator is the RandomState instance used by np.random.
#'   Note that different initializations might result in different local minima
#'   of the cost function.
#' @param method character (default: ‘barnes.hut’) By default the gradient
#'   calculation algorithm uses Barnes-Hut approximation running in O(NlogN)
#'   time. method=’exact’ will run on the slower, but exact, algorithm in O(N^2)
#'   time. The exact algorithm should be used when nearest-neighbor errors need
#'   to be better than 3%. However, the exact method cannot scale to millions of
#'   examples.
#' @param angle numeric (default: 0.5) Only used if method=’barnes.hut’ This is
#'   the trade-off between speed and accuracy for Barnes-Hut T-SNE. ‘angle’ is
#'   the angular size (referred to as theta in [3]) of a distant node as
#'   measured from a point. If this size is below ‘angle’ then it is used as a
#'   summary node of all points contained within it. This method is not very
#'   sensitive to changes in this parameter in the range of 0.2 - 0.8. Angle
#'   less than 0.2 has quickly increasing computation time and angle greater 0.8
#'   has quickly increasing error.#' @param cheat_metric
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
                          early.exaggeration = 12.0,
                          learning.rate = 200.0,
                          n.iter = 1000,
                          n.iter.without.progress = 300,
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
  if(!py_module_available('MulticoreTSNE')){
    stop("The MulticoreTSNE module is unavailable.  Please activate the appropriate environment or install the module.")
  }

  mctsne.module <- import(module = 'MulticoreTSNE', delay_load = TRUE)
  if (is.null(n.jobs)){
    n.jobs <- detectCores()
  }
  mctsne <- mctsne.module$MulticoreTSNE(n_components = as.integer(n.components),
                                        perplexity = as.numeric(perplexity),
                                        early_exaggeration = as.numeric(early.exaggeration),
                                        learning_rate = as.numeric(learning.rate),
                                        n_iter = as.integer(n.iter),
                                        n_iter_without_progress = as.integer(n.iter.without.progress),
                                        min_grad_norm = as.numeric(min.grad.norm),
                                        metric = metric,
                                        init = init,
                                        verbose = as.integer(verbose),
                                        random_state = random.state,
                                        method = method,
                                        angle = as.numeric(angle),
                                        n_jobs = n.jobs,
                                        cheat_metric = cheat.metric)
  mctsne.df <- mctsne$fit_transform(r.data.frame)
  return(mctsne.df)
}

#' FIt-SNE
#'
#' An R wrapper for the FIt-SNE Python module found at
#' https://github.com/KlugerLab/FIt-SNE
#'
#' @param X: np.ndarray, shape (samples x dimensions) Input data.
#' @param no_dims: integer (default = 2) Dimensionality of the embedding
#' @param perplexity: numeric (default = 30.0) Perplexity is used to determine
#'   the bandwidth of the Gaussian kernel in the input space. Set to -1 if using
#'   fixed sigma/K (see below)
#' @param sigma: numeric (default = -30) Fixed bandwidth of Gaussian kernel to
#'   use in lieu of the perplexity-based adaptive kernel width typically used in
#'   t-SNE
#' @param K: integer (default = -1) Number of nearest neighbors to get when
#'   using fixed sigma in lieu of perplexity-based adaptive kernel width
#'   typically used in t-SNE
#' @param initialization: np.ndarray, shape (samples x no_dims) (default = NULL)
#'   Initialization of the embedded points to use in lieu of the random
#'   initialization typically used in t-SNE
#' @param theta: numeric (default = 0.5) Set to 0 for exact.  If non-zero, then
#'   will use either Barnes Hut or FIt-SNE based on `fft_not_bh`. If Barnes Hut,
#'   then this determines the accuracy of BH approximation.
#' @param rand_seed: integer (default = -1) Random seed to get deterministic
#'   output
#' @param max_iter: int, (default = 1000) Number of iterations of t-SNE to run.
#' @param stop_lying_iter: int, (default = 200) When to switch off early
#'   exaggeration.
#' @param fft_not_bh: boolean, (default = FALSE) If theta is nonzero, this
#'   determins whether to use FIt-SNE or Barnes Hut approximation.
#' @param ann_not_vptree: boolean, (default=FALSE) This determines whether to
#'   use aproximate (Annoy) or deterministic (vptree) nearest neighbours
#' @param early_exag_coeff: numeric (default = 12.0) When to switch off early
#'   exaggeration. (>1)
#' @param no_momentum_during_exag: boolean (default = FALSE) Set to 0 to use
#'   momentum and other optimization tricks. 1 to do plain, vanilla gradient
#'   descent (useful for testing large exaggeration coefficients)
#' @param start_late_exag_iter: int, (default=-1) When to start late
#'   exaggeration. Set to -1 to not use late exaggeration
#' @param late_exag_coeff: float, (default=-1) Late exaggeration coefficient.
#'   Set to -1 to not use late exaggeration.
#' @param n_trees: integer (default = 50) ANNOY parameter
#' @param search_k: integer (default = -1) ANNOY parameter
#' @param nterms: integer (default = 3) If using FIt-SNE, this is the number of
#'   interpolation points per sub-interval
#' @param intervals_per_integer: numeric (default = 1) See min_num_intervals
#' @param min_num_intervals: integer (default = 50) Let maxloc =
#'   ceil(max(max(X))) and minloc = floor(min(min(X))). i.e. the points are in a
#'   [minloc]^no_dims by [maxloc]^no_dims interval/square.  The number of
#'   intervals in each dimension is either min_num_intervals or ceil((maxloc
#'   -minloc)/opts.intervals_per_integer), whichever is larger.
#'   opts.min_num_intervals must be an integer >0, and
#'   opts.intervals_per_integer must be >0.
#' @param nthreads: unsigned integer (default = 0) Number of threads to be used
#'   in computation of input similarities (both for vptrees and ann). 0 uses the
#'   maximum number of threads supported by the hardware.
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#'
#' @return dataframe
#' @export
#'
#' @examples
fitsne <- function(r.data.frame,
                   no_dims=2,
                   perplexity=30.0,
                   sigma=-30.0,
                   K=-1,
                   initialization=NULL,
                   theta=0.5,
                   rand_seed=-1,
                   max_iter=1000,
                   stop_lying_iter=200,
                   fft_not_bh=TRUE,
                   ann_not_vptree=TRUE,
                   early_exag_coeff=12.0,
                   no_momentum_during_exag=FALSE,
                   start_late_exag_iter=-1,
                   late_exag_coeff=-1,
                   n_trees=50,
                   search_k=-1,
                   nterms=3,
                   intervals_per_integer=1,
                   min_num_intervals=50,
                   nthreads=0
){
  if(!py_module_available('fitsne')){
    stop("The fitsne module is unavailable.
         Please activate the appropriate environment or install the module.")
  }

  fitsne.module <- import(module = 'fitsne', delay_load = TRUE)
  numpy.module <- import(module = 'numpy', delay_load = TRUE)
  if (is.null(nthreads)){
    nthreads <- detectCores()
  }
  fitsne.df <- fitsne.module$FItSNE(X = numpy.module$copy(r.data.frame, order = "C"),
                                    no_dims = as.integer(no_dims),
                                    perplexity = as.numeric(perplexity),
                                    sigma = as.numeric(sigma),
                                    K = as.integer(K),
                                    initialization = initialization,
                                    theta = as.numeric(theta),
                                    rand_seed = as.integer(rand_seed),
                                    max_iter = as.integer(max_iter),
                                    stop_lying_iter = as.integer(stop_lying_iter),
                                    fft_not_bh = fft_not_bh,
                                    ann_not_vptree = ann_not_vptree,
                                    early_exag_coeff = as.numeric(early_exag_coeff),
                                    no_momentum_during_exag = no_momentum_during_exag,
                                    start_late_exag_iter = as.integer(start_late_exag_iter),
                                    late_exag_coeff = as.numeric(late_exag_coeff),
                                    n_trees = as.integer(n_trees),
                                    search_k = as.integer(search_k),
                                    nterms = as.integer(nterms),
                                    intervals_per_integer = as.numeric(intervals_per_integer),
                                    min_num_intervals = as.integer(min_num_intervals),
                                    nthreads = as.integer(nthreads))
  return(fitsne.df)
}

#' umap
#'
#' An R wrapper around the UMAP Python module found at
#' https://github.com/lmcinnes/umap
#'
#' Uniform Manifold Approximation and Projection
#'
#' Finds a low dimensional embedding of the data that approximates an underlying
#' manifold.
#'
#' @param r.data.frame
#' @param n.neighbors numeric (optional, default 30) The size of local
#'   neighborhood (in terms of number of neighboring sample points) used for
#'   manifold approximation. Larger values result in more global views of the
#'   manifold, while smaller values result in more local data being preserved.
#'   In general values should be in the range 2 to 100.
#' @param n.components integer (optional, default 3) The dimension of the space
#'   to embed into. This defaults to 3 to provide easy visualization, but can
#'   reasonably be set to any integer value in the range 2 to 100.
#' @param metric character or function (optional, default 'euclidean') The metric
#'   to use to compute distances in high dimensional space. If a character is
#'   passed it must match a valid predefined metric. If a general metric is
#'   required a function that takes two 1d arrays and returns a numeric can be
#'   provided. For performance purposes it is required that this be a numba
#'   jit'd function. Valid character metrics include: euclidean. manhattan,
#'   chebyshev, minkowski, canberra, braycurtis, mahalanobis, wminkowski,
#'   seuclidean, cosine, correlation, haversine, hamming, jaccard, dice,
#'   russelrao, kulsinski, rogerstanimoto, sokalmichener, sokalsneath, yule/
#'   Metrics that take arguments (such as minkowski, mahalanobis etc.) can have
#'   arguments passed via the metric_kwds dictionary. At this time care must be
#'   taken and dictionary elements must be ordered appropriately; this will
#'   hopefully be fixed in the future.
#' @param negative.sample.rate integer (optional, default 5) The number of
#'   negative edge/1-simplex samples to use per positive edge/1-simplex sample
#'   in optimizing the low dimensional embedding.
#' @param alpha numeric (optional, default 1.0) The initial learning rate for
#'   the embedding optimization.
#' @param init character (optional, default 'spectral') How to initialize the low
#'   dimensional embedding. Options are: * 'spectral': use a spectral embedding
#'   of the fuzzy 1-skeleton * 'random': assign initial embedding positions at
#'   random. * A numpy array of initial embedding positions.
#' @param min.dist numeric (optional, default 0.1) The effective minimum
#'   distance between embedded points. Smaller values will result in a more
#'   clustered/clumped embedding where nearby points on the manifold are drawn
#'   closer together, while larger values will result on a more even dispersal
#'   of points. The value should be set elative to the "spread" value, which
#'   determines the scale at which embedded points will be spread out.
#' @param spread numeric (optional, default 1.0) The effective scale of embedded
#'   points. In combination with "min_dist" this determines how
#'   clustered/clumped the embedded points are.
#' @param set.op.mix.ratio numeric (optional, default 1.0) Interpolate between
#'   (fuzzy) union and intersection as the set operation used to combine local
#'   fuzzy simplicial sets to obtain a global fuzzy simplicial sets. Both fuzzy
#'   set operations use the product t-norm. The value of this parameter should
#'   be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while
#'   0.0 will use a pure fuzzy intersection.
#' @param local.connectivity integer (optional, default 1) The local
#'   connectivity required -- i.e. the number of nearest neighbors that should
#'   be assumed to be connected at a local level. The higher this value the more
#'   connected the manifold becomes locally. In practice this should be not more
#'   than the local intrinsic dimension of the manifold.
#' @param gamma numeric (optional, default 1.0) Weighting applied to negative
#'   samples in low dimensional embedding optimization. Values higher than one
#'   will result in greater weight being given to negative samples.
#' @param bandwidth numeric (optional, default 1.0) The effective bandwidth of
#'   the kernel if we view the algorithm as similar to Laplacian eigenmaps.
#'   Larger values induce more connectivity and a more global view of the data,
#'   smaller values concentrate more locally.
#' @param random.state int, RandomState instance or NULL (default: NULL) If int,
#'   random_state is the seed used by the random number generator; If
#'   RandomState instance, random_state is the random number generator; If NULL,
#'   the random number generator is the RandomState instance used by
#'   `np.random`.
#' @param angular.rp.forest boolean (optional, default FALSE) Whether to use an
#'   angular random projection forest to initialise the approximate nearest
#'   neighbor search. This can be faster, but is mostly on useful for metric
#'   that use an angular style distance such as cosine, correlation etc. In the
#'   case of those metrics angular forests will be chosen automatically.
#' @param verbose: boolean (optional, default FALSE) Controls verbosity of
#'   logging.
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
  if(!py_module_available('umap')){
    stop("The umap module is unavailable.  Please activate the appropriate environment or install the module.")
  }
  umap.module <- import(module = 'umap', delay_load = TRUE)
  umap.embed <- umap.module$UMAP(n_neighbors = as.integer(n.neighbors),
                                 min_dist = as.numeric(min.dist),
                                 n_components = as.integer(n.components),
                                 metric = metric,
                                 init = init,
                                 alpha = as.numeric(alpha),
                                 spread = as.numeric(spread),
                                 bandwidth = as.numeric(bandwidth),
                                 random_state = random.state,
                                 angular_rp_forest = angular.rp.forest,
                                 verbose = verbose,
                                 set_op_mix_ratio = as.numeric(set.op.mix.ratio),
                                 gamma = as.numeric(gamma),
                                 negative_sample_rate = negative.sample.rate)
  umap.df <- umap.embed$fit_transform(r.data.frame)
  return(umap.df)
}

#' phate
#'
#' An R wrapper around the PHATE Python module found at
#' https://github.com/KrishnaswamyLab/PHATE
#'
#' Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)
#' embeds high dimensional single-cell data into two or three dimensions for
#' visualization of biological progressions as described in Moon et al, 2017
#'
#' @param n_components integer (default: 2) number of dimensions in which the
#'   data will be embedded
#' @param k integer (default: 5) number of nearest neighbors on which to build
#'   kernel
#' @param a integer, (default: 15) sets decay rate of kernel tails. If NULL,
#'   alpha decaying kernel is not used
#' @param n_landmark integer, default: 2000 number of landmarks to use in fast
#'   PHATE
#' @param t integer (default: 'auto') power to which the diffusion operator is
#'   powered. This sets the level of diffusion. If 'auto', t is selected
#'   according to the knee point in the Von Neumann Entropy of the diffusion
#'   operator
#' @param gamma numeric (default: 1) Informational distance constant between -1
#'   and 1. "gamma=1" gives the PHATE log potential, "gamma=0` gives a square
#'   root potential.
#' @param n_pca integer (default: 100) Number of principal components to use for
#'   calculating neighborhoods. For extremely large datasets, using n_pca < 20
#'   allows neighborhoods to be calculated in roughly log(n_samples) time.
#' @param knn_dist character (default: 'euclidean') recommended values:
#'   'euclidean', 'cosine', 'precomputed' Any metric from
#'   "scipy.spatial.distance` can be used distance metric for building kNN
#'   graph. If 'precomputed', "data` should be an n_samples x n_samples distance
#'   or affinity matrix. Distance matrices are assumed to have zeros down the
#'   diagonal, while affinity matrices are assumed to have non-zero values down
#'   the diagonal. This is detected automatically using "data[0,0]`. You can
#'   override this detection with "knn_dist='precomputed_distance'` or
#'   "knn_dist='precomputed_affinity'`.
#' @param mds_dist character (default: 'euclidean') recommended values:
#'   'euclidean' and 'cosine' Any metric from "scipy.spatial.distance` can be
#'   used distance metric for MDS
#' @param mds character (default: 'metric') choose from ['classic', 'metric',
#'   'nonmetric']. Selects which MDS algorithm is used for dimensionality
#'   reduction
#' @param n_jobs integer (default: 1) The number of jobs to use for the
#'   computation. If -1 all CPUs are used. If 1 is given, no parallel computing
#'   code is used at all, which is useful for debugging. For n_jobs below -1,
#'   (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are
#'   used
#' @param random_state integer or numpy.RandomState (default: NULL) The
#'   generator used to initialize SMACOF (metric, nonmetric) MDS If an integer
#'   is given, it fixes the seed defaults to the global "numpy` random number
#'   generator
#' @param verbose integer or boolean (default: 1) If "TRUE` or "> 0", print
#'   status messages
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#'
#' @return
#' @export
#'
#' @examples
phate <- function(n.components = 3,
                  k = 5,
                  a = 15,
                  n.landmark = 2000,
                  t = 'auto',
                  gamma = 1,
                  n.pca = 100,
                  knn.dist = 'euclidean',
                  mds.dist = 'euclidean',
                  mds = 'metric',
                  n.jobs = NULL,
                  random.state = NULL,
                  verbose = 1){
  if(!py_module_available('phate')){
    stop("The phate module is unavailable.  Please activate the appropriate environment or install the module.")
  }
  phate.module <- import(module = 'phate', delay_load = TRUE)
  if (is.null(n.jobs)){
    n.jobs <- detectCores()
  }
  phate.embed <- phate.module$PHATE(n_components = as.integer(n.components),
                                    k = as.integer(k),
                                    a = as.integer(a),
                                    n_landmark = as.integer(n.landmark),
                                    t = as.integer(t),
                                    gamma = as.numeric(gamma),
                                    n_pca = as.integer(n.pca),
                                    knn_dist = knn.dist,
                                    mds_dist = mds.dist,
                                    mds = mds,
                                    n_jobs = as.integer(n.jobs),
                                    random_state = random.state,
                                    verbose = as.integer(verbose)
  )
  phate.df <- phate.embed$fit_transform(r.data.frame)
  return(phate.df)
}

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

#' Deep count autoencoder
#'
#' Performs denoising of gene expression using dca
#' (https://github.com/theislab/dca)
#' Some arguments have been removed due to being troublesome in the R-Python
#' interface
#'
#' @param exprDat Unnomalized, unscaled gene expression matrix, with cell names
#'   as rows and gene names as columns
#' @param mode character. "denoise"(default), or "latent".  "denoise"
#'   overwrites "adata.X" with denoised expression values.  In "latent" mode DCA
#'   adds "adata.obsm["X_dca"]" to given adata object. This matrix represent
#'   latent representation of cells via DCA.
#' @param ae_type character. "zinb-conddisp"(default), "zinb",
#'   "nb-conddisp" or "nb".  Type of the autoencoder. Return values and the
#'   architecture is determined by the type e.g. "nb" does not provide dropout
#'   probabilities.
#' @param normalize_per_cell boolean (default: TRUE).  If true,
#'   library size normalization is performed using the
#'   "sc.pp.normalize_per_cell" function in Scanpy and saved into adata object.
#'   Mean layer is re-introduces library size differences by scaling the mean
#'   value of each cell in the output layer. See the manuscript for more
#'   details.
#' @param scale boolean (default: TRUE).  If true, the input of the
#'   autoencoder is centered using "sc.pp.scale" function of Scanpy. Note that
#'   the output is kept as raw counts as loss functions are designed for the
#'   count data.
#' @param log1p boolean (default: TRUE).  If true, the input of the
#'   autoencoder is log transformed with a pseudocount of one using
#'   "sc.pp.log1p" function of Scanpy.
#' @param batchnorm boolean (default: TRUE).  If true, batch
#'   normalization is performed.
#' @param activation str (default: "relu").  Activation function of
#'   hidden layers.
#' @param init str (default: "glorot_uniform").  Initialization
#'   method used to initialize weights.
#' @param epochs integer (default: 300).  Number of total epochs in
#'   training.
#' @param reduce_lr integer (default: 10).  Reduces learning rate if
#'   validation loss does not improve in given number of epochs.
#' @param early_stop integer (default: 15).  Stops training if
#'   validation loss does not improve in given number of epochs.
#' @param batch_size integer (default: 32).  Number of samples in the
#'   batch used for SGD.
#' @param optimizer str (default: "rmsprop").  Type of optimization
#'   method used for training.
#' @param random_state integer (default: 0).  Seed for python, numpy
#'   and tensorflow.
#' @param threads integer or NULL (default: NULL).  Number of threads
#'   to use in training. All cores are used by default.
#' @param verbose boolean (default: FALSE).  If true, prints
#'   additional information about training and architecture.
#' @param return_model boolean (default: FALSE).  If true, trained
#'   autoencoder object is returned. See "Returns".
#' @param return_info boolean (default: FALSE).  If true, all
#'   additional parameters of DCA are stored in "adata.obsm" such as dropout
#'   probabilities (obsm["X_dca_dropout"]) and estimated dispersion values
#'   (obsm["X_dca_dispersion"]), in case that autoencoder is of type zinb or
#'   zinb-conddisp.
#' @param copy boolean (default: FALSE). If true, a copy of anndata
#'   is returned.
#'
#' @importFrom reticulate import py_module_available dict
#' @importFrom parallel detectCores
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
dca <- function(exprDat,
                mode = 'denoise',
                ae_type = 'zinb-conddisp',
                normalize_per_cell = TRUE,
                scale = TRUE,
                log1p = TRUE,
                batchnorm = TRUE,
                activation = 'relu',
                init = 'glorot_uniform',
                epochs = 300,               # training args
                reduce_lr = 10,
                early_stop = 15,
                batch_size = 32,
                optimizer = 'rmsprop',
                random_state = 0,
                threads = NULL,
                verbose = TRUE,
                return_model = FALSE,
                return_info = FALSE,
                copy = FALSE){

  if(!py_module_available('dca')){
    stop("The DCA module is unavailable.
         Please activate the appropriate environment or
         install the module.")
  }
  if(!py_module_available('scanpy')){
    stop("The scanpy module is unavailable.
         Please activate the appropriate environment or
         install the module.")
  }

  if(is.null(threads)){
    threads <- detectCores()
  }

  dca.module <- import(module = 'dca.api', delay_load = TRUE)
  scanpy.module <- import(module = 'scanpy.api', delay_load = TRUE)

  cell.names <- rownames(exprDat)
  gene.names <- colnames(exprDat)

  adata <- scanpy.module$AnnData(exprDat,
                                 obs = cell.names,
                                 var = gene.names)
  adata$obs_names <- cell.names
  adata$var_names <- gene.names
  scanpy.module$pp$filter_genes(data = adata,
                                min_counts = as.integer(1))
  scanpy.module$pp$filter_genes(data = adata,
                                min_cells = as.integer(3))

  dca.module$dca(adata = adata,
                 mode = mode,
                 ae_type = ae_type,
                 normalize_per_cell = normalize_per_cell,
                 scale = scale,
                 log1p = log1p,
                 batchnorm = batchnorm,
                 activation = activation,
                 init = init,
                 epochs = as.integer(epochs),
                 reduce_lr = as.integer(reduce_lr),
                 early_stop = as.integer(early_stop),
                 batch_size = as.integer(batch_size),
                 optimizer = optimizer,
                 random_state = as.integer(random_state),
                 threads = as.integer(threads),
                 verbose = verbose,
                 return_model = return_model,
                 return_info = return_info,
                 copy = copy)

  conv <- adata$data %>% t()
  colnames(conv) <- adata$obs_names$values
  rownames(conv) <- glue("DCA_{adata$var_names$values}")
  return(conv)
}
