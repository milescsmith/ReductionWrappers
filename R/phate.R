#' phate
#'
#' An R wrapper around the PHATE Python module found at
#' https://github.com/KrishnaswamyLab/PHATE
#'
#' Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)
#' embeds high dimensional single-cell data into two or three dimensions for
#' visualization of biological progressions as described in Moon et al, 2017
#'
#' @param n_components integer number of dimensions in which the
#' data will be embedded. Default: 3
#' @param k integer number of nearest neighbors on which to build
#' kernel. Default: 5
#' @param a integer sets decay rate of kernel tails. If NULL,
#' alpha decaying kernel is not used.
#' @param n_landmark integer, default: 2000 number of landmarks to use in fast
#' PHATE. Default: 15
#' @param t integer power to which the diffusion operator is
#' powered. This sets the level of diffusion. If "auto", t is selected
#' according to the knee point in the Von Neumann Entropy of the diffusion
#' operator. Default: "auto"
#' @param gamma numeric Informational distance constant between -1
#' and 1. "gamma=1" gives the PHATE log potential, "gamma=0" gives a square
#' root potential. Default: 1
#' @param n_pca integer Number of principal components to use for
#' calculating neighborhoods. For extremely large datasets, using n_pca < 20
#' allows neighborhoods to be calculated in roughly log(n_samples) time. Default: 100)
#' @param knn_dist character recommended values:
#' 'euclidean', 'cosine', 'precomputed' Any metric from
#' "scipy.spatial.distance" can be used distance metric for building kNN
#' graph. If 'precomputed', "data" should be an n_samples x n_samples distance
#' or affinity matrix. Distance matrices are assumed to have zeros down the
#' diagonal, while affinity matrices are assumed to have non-zero values down
#' the diagonal. This is detected automatically using "data[0,0]". You can
#' override this detection with "knn_dist='precomputed_distance'" or
#' "knn_dist='precomputed_affinity'". Default: 'euclidean'
#' @param mds_dist character recommended values:
#' 'euclidean' and 'cosine' Any metric from "scipy.spatial.distance" can be
#' used distance metric for MDS. Default: 'euclidean'
#' @param mds character choose from ['classic', 'metric',
#' 'nonmetric']. Selects which MDS algorithm is used for dimensionality
#' reduction. Default: 'metric'
#' @param n_jobs integer The number of jobs to use for the
#' computation. If -1 all CPUs are used. If 1 is given, no parallel computing
#' code is used at all, which is useful for debugging. For n_jobs below -1,
#' (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are
#' used. Default: 1
#' @param random_state integer or numpy.RandomState The
#' generator used to initialize SMACOF (metric, nonmetric) MDS If an integer
#' is given, it fixes the seed defaults to the global "numpy" random number
#' generator. Default: NULL
#' @param verbose integer or boolean If "TRUE" or "> 0", print
#' status messages. Default: 1
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#'
#' @return
#' @export
#'
#' @examples
phate <- function(n_components = 3,
                  k = 5,
                  a = 15,
                  n_landmark = 2000,
                  t = 'auto',
                  gamma = 1,
                  n_pca = 100,
                  knn_dist = 'euclidean',
                  mds_dist = 'euclidean',
                  mds = 'metric',
                  n_jobs = NULL,
                  random_state = NULL,
                  verbose = 1){
  if(!py_module_available('phate')){
    stop("The phate module is unavailable.  Please activate the appropriate environment or install the module.")
  }
  phate.module <- import(module = 'phate', delay_load = TRUE)
  if (is.null(n_jobs)){
    n_jobs <- detectCores()
  }
  phate.embed <- phate.module$PHATE(n_components = as.integer(n_components),
                                    k = as.integer(k),
                                    a = as.integer(a),
                                    n_landmark = as.integer(n_landmark),
                                    t = as.integer(t),
                                    gamma = as.numeric(gamma),
                                    n_pca = as.integer(n_pca),
                                    knn_dist = knn_dist,
                                    mds_dist = mds_dist,
                                    mds = mds,
                                    n_jobs = as.integer(n_jobs),
                                    random_state = random_state,
                                    verbose = as.integer(verbose)
  )
  phate.df <- phate.embed$fit_transform(r.data.frame)
  return(phate.df)
}
