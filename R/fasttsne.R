#' fastTSNE
#'
#' An R wrapper for the fastTSNE Python module found at
#' https://github.com/pavlin-policar/fastTSNE/
#'
#' @param n_components: The dimension of the embedding space. Default: 2
#' @param perplexity: Perplexity can be thought of as the continuous :math:`k` number of
#' neighbors to consider for each data point. To avoid confusion, note that
#' perplexity linearly impacts runtime. Default: 30
#' @param learning_rate: The learning rate for the t-SNE optimization steps. Typical values range
#' from 1 to 1000. Setting the learning rate too low or too high may result
#' in the points forming a "ball". This is also known as the crowding
#' problem. Default: 100
#' @param early_exaggeration_iter: The number of iterations to run in the *early exaggeration* phase. Default: 250
#' @param early_exaggeration: The early exaggeration factor. Typical values range from 12 to 32,
#' however larger values have also been found to work well with specific
#' values of learning rate. Default: 12
#' @param n_iter: The number of iterations to run in the normal optimization regime. Default: 750
#' @param exaggeration: The exaggeration factor to be used during the normal optmimization
#' phase. Standard implementation don't use this exaggeration and it
#' typically isn't necessary for smaller data sets, but it has been shown
#' that for larger data sets, using some exaggeration is necessary in order
#' to obtain good embeddings. Default: NULL
#' @param theta: Only used when "negative_gradient_method='bh'" or its other aliases.
#' This is the trade-off parameter between speed and accuracy of the tree
#' approximation method. Typical values range from 0.2 to 0.8. The value 0
#' indicates that no approximation is to be made and produces exact results
#' also producing longer runtime. Default: 0.5
#' @param n_interpolation_points: Only used when "negative_gradient_method='fft'" or its other aliases.
#' The number of interpolation points to use within each grid cell for
#' interpolation based t-SNE. It is highly recommended leaving this value
#' at the default 3 as otherwise the interpolation may suffer from the
#' Runge phenomenon. Theoretically, increasing this number will result in
#' higher approximation accuracy, but practically, this can also be done
#' with the "ints_in_interval" parameter, which does not suffer from the
#' Runge phenomenon and should always be preferred. Default: 3
#' @param min_num_intervals: Only used when "negative_gradient_method='fft'" or its other aliases.
#' The interpolation approximation method splits the embedding space into a
#' grid, where the number of grid cells is governed by
#' "ints_in_interval". Sometimes, especially during early stages of
#' optimization, that number may be too small, and we may want better
#' accuracy. The number of intervals used will thus always be at least the
#' number specified here. Note that larger values will produce more precise
#' approximations but will have longer runtime. Default: 10
#' @param ints_in_interval: Only used when "negative_gradient_method='fft'" or its other aliases.
#' Since the coordinate range of the embedding changes during optimization,
#' this value tells us how many integers should appear in a single e.g.
#' setting this value to 3 means that the intervals will appear as follows:
#' [0, 3), [3, 6), ... Lower values will need more intervals to fill the
#' space, e.g. 1.5 will produce 4 intervals [0, 1.5), [1.5, 3), ...
#' Therefore lower values will produce more intervals, producing more
#' interpolation points which in turn produce better approximation at the
#' cost of longer runtime. Default: 2
#' @param initialization: The initial point positions to be used in the embedding space. Can be a
#' precomputed numpy array, "pca" or "random". Please note that when
#' passing in a precomputed positions, it is highly recommended that the
#' point positions have small variance (var(Y) < 0.0001), otherwise you may
#' get poor embeddings. Default: 'pca'
#' @param metric: The metric to be used to compute affinities between points in the
#' original space. Default: 'euclidean'
#' @param metric_params: Additional keyword arguments for the metric function. Default: NULL
#' @param initial_momentum: t-SNE optimization uses momentum for faster convergence. This value
#' controls the momentum used during the *early optimization* phase. Default: 0.5
#' @param final_momentum: t-SNE optimization uses momentum for faster convergence. This value
#' controls the momentum used during the normal regime and *late exaggeration* phase. Default: 0.5
#' @param n_jobs: The number of jobs to run in parallel. This follows the scikit-learn
#' convention, "-1" meaning all processors, "-2" meaning all but one
#' processor and so on. Default: all cores.
#' @param neighbors: Specifies the nearest neighbor method to use. Can be either "exact" or
#' "approx". "exact" uses space partitioning binary trees from
#' scikit-learn while "approx" makes use of nearest neighbor descent.
#' Note that "approx" has a bit of overhead and will be slower on smaller
#' data sets than exact search. Default: 'exact'
#' @param negative_gradient_method: Specifies the negative gradient approximation method to use. For smaller
#' data sets, the Barnes-Hut approximation is appropriate and can be
#' set using one of the following aliases: "bh", "BH" or
#' "barnes-hut". Note that the time complexity of Barnes-Hut scales as
#' :math:`\mathcal{O}(N \log N)`. For larger data sets, the FFT accelerated
#' interpolation method is more appropriate and can be set using one of the
#' following aliases: "fft", "FFT" or "ìnterpolation". Note that
#' this method scales linearly in the number of points
#' :math:`\mathcal{O}(N)` and its complexity is governed by the number of
#' interpolation points. Default: 'fft'
#' @param callbacks: We can pass callbacks, that will be run every "callbacks_every_iters"
#' iterations. Each callback should accept three parameters, the first is
#' the current iteration number, the second is the current KL divergence
#' error and the last is the current embedding. The callback may return
#' "True" in order to stop the optimization. Default: NULL
#' @param callbacks_every_iters: How many iterations should run between each time a callback is invoked. Default: 50
#' @param random_state: The random state parameter follows the convention used in scikit-learn.
#' If the value is an int, random_state is the seed used by the random
#' number generator. If the value is a RandomState instance, then it will
#' be used as the random number generator. If the value is None, the random
#' number generator is the RandomState instance used by `np.random`. Default: NULL
#'
#' @importFrom reticulate import
#' @importFrom parallel detectCores
#' @importFrom glue glue
#'
#' @return dataframe
#' @export
#'
#' @examples
fastTSNE <- function(r_data_frame,
                     n_components = 2,
                     perplexity = 30,
                     learning_rate = 100,
                     early_exaggeration_iter = 250,
                     early_exaggeration = 12,
                     n_iter = 750,
                     exaggeration = NULL,
                     theta = 0.5,
                     n_interpolation_points = 3,
                     min_num_intervals = 10,
                     ints_in_interval = 2,
                     initialization = 'pca',
                     metric = 'euclidean',
                     metric_params = NULL,
                     initial_momentum = 0.5,
                     final_momentum = 0.8,
                     n_jobs = NULL,
                     neighbors = 'exact',
                     negative_gradient_method = 'fft',
                     callbacks = NULL,
                     callbacks_every_iters = 50,
                     random_state = NULL){
  if(!py_module_available('fastTSNE')){
    stop("The fitsne module is unavailable.
         Please activate the appropriate environment or install the module.")
  }

  fastTSNE.module <- import(module = 'fastTSNE', delay_load = TRUE)
  # numpy.module <- import(module = 'numpy', delay_load = TRUE)
  if (is.null(n_jobs)){
    n_jobs <- detectCores()
  }
  # X = numpy.module$copy(r.data.frame, order = "C")
  tsne <- fastTSNE.module$TSNE(n_components = as.integer(n_components),
                                    perplexity = as.numeric(perplexity),
                                    learning_rate = as.numeric(learning_rate),
                                    early_exaggeration_iter = as.integer(early_exaggeration_iter),
                                    early_exaggeration = as.numeric(early_exaggeration),
                                    n_iter = as.integer(n_iter),
                                    exaggeration = exaggeration,
                                    theta = as.numeric(theta),
                                    n_interpolation_points = as.integer(n_interpolation_points),
                                    min_num_intervals = as.integer(min_num_intervals),
                                    ints_in_interval = as.numeric(ints_in_interval),
                                    initialization = initialization,
                                    metric = as.character(metric),
                                    metric_params = metric_params,
                                    initial_momentum = as.numeric(initial_momentum),
                                    final_momentum = as.numeric(final_momentum),
                                    n_jobs = as.integer(n_jobs),
                                    neighbors = as.character(neighbors),
                                    negative_gradient_method = as.character(negative_gradient_method),
                                    callbacks = callbacks,
                                    callbacks_every_iters = as.integer(callbacks_every_iters),
                                    random_state = random_state)

  fasttsne.df = tsne$fit(X = r_data_frame)
  rownames(fasttsne.df) <- rownames(r_data_frame)
  colnames(fasttsne.df) <- glue('tsne_{1:n_components}')
  return(fasttsne.df)
}
