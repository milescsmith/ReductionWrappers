#' phate
#'
#' An R wrapper around the PHATE Python module found at
#' https://github.com/KrishnaswamyLab/PHATE
#'
#' Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)
#' embeds high dimensional single-cell data into two or three dimensions for
#' visualization of biological progressions as described in Moon et al, 2017
#'
#' @param ce cell embeddings
#' @param n_components integer (default: 3) number of dimensions in which the
#' data will be embedded
#' @param k integer (default: 5) number of nearest neighbors on which to build
#' kernel
#' @param a integer, (default: 15) sets decay rate of kernel tails. If NULL,
#' alpha decaying kernel is not used
#' @param n_landmark integer, default: 2000 number of landmarks to use in fast
#' PHATE
#' @param t integer (default: 'auto') power to which the diffusion operator is
#' powered. This sets the level of diffusion. If 'auto', t is selected
#' according to the knee point in the Von Neumann Entropy of the diffusion
#' operator
#' @param gamma numeric (default: 1) Informational distance constant between -1
#' and 1. "gamma=1" gives the PHATE log potential, "gamma=0" gives a square
#' root potential.
#' @param n_pca integer (default: 100) Number of principal components to use for
#' calculating neighborhoods. For extremely large datasets, using n_pca < 20
#' allows neighborhoods to be calculated in roughly log(n_samples) time.
#' @param knn_dist character (default: 'euclidean') recommended values:
#' 'euclidean', 'cosine', 'precomputed' Any metric from
#' "scipy.spatial.distance" can be used distance metric for building kNN
#' graph. If 'precomputed', "data" should be an n_samples x n_samples distance
#' or affinity matrix. Distance matrices are assumed to have zeros down the
#' diagonal, while affinity matrices are assumed to have non-zero values down
#' the diagonal. This is detected automatically using "data[0,0]". You can
#' override this detection with "knn_dist='precomputed_distance'" or
#' "knn_dist='precomputed_affinity'".
#' @param mds_dist character (default: 'euclidean') recommended values:
#' 'euclidean' and 'cosine' Any metric from "scipy.spatial.distance" can be
#' used distance metric for MDS
#' @param mds character (default: 'metric') choose from ['classic', 'metric',
#' 'nonmetric']. Selects which MDS algorithm is used for dimensionality
#' reduction
#' @param n_jobs integer (default: 1) The number of jobs to use for the
#' computation. If -1 all CPUs are used. If 1 is given, no parallel computing
#' code is used at all, which is useful for debugging. For n_jobs below -1,
#' (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are
#' used
#' @param random_state integer or numpy.RandomState (default: NULL) The
#' generator used to initialize SMACOF (metric, nonmetric) MDS If an integer
#' is given, it fixes the seed defaults to the global "numpy" random number
#' generator
#' @param verbose integer or boolean (default: 1) If "TRUE" or "> 0", print
#' status messages
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#'
#' @return
#' @export
#'
#' @examples
phate <- function(ce,
                  n_components = 3,
                  k = 5,
                  a = 15,
                  n_landmark = 2000,
                  t = "auto",
                  gamma = 1,
                  n_pca = NULL,
                  knn_dist = "euclidean",
                  mds_dist = "euclidean",
                  mds = "metric",
                  n_jobs = NULL,
                  random_state = NULL,
                  verbose = 1){
  if (!py_module_available("phate")){
    stop("The phate module is unavailable.  Please activate the appropriate environment or install the module.")
  }
  phate <- import(module = "phate", delay_load = TRUE)
  if (is.null(n_jobs)){
    n_jobs <- detectCores()
  }


  if (is.null(n_pca)){
    n_pca <- ncol(ce)
  }
  phate_op <- phate$PHATE(n_components = as.integer(n_components),
                          k = as.integer(k),
                          a = as.integer(a),
                          n_landmark = as.integer(n_landmark),
                          t = as.character(t),
                          gamma = as.numeric(gamma),
                          n_pca = as.integer(n_pca),
                          knn_dist = knn_dist,
                          mds_dist = mds_dist,
                          mds = mds,
                          n_jobs = as.integer(n_jobs),
                          random_state = random_state,
                          verbose = as.integer(verbose)
  )
  phate_df <- phate_op$fit_transform(ce)
  return(phate_df)
}
