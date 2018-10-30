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
