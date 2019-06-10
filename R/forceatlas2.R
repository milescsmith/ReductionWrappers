#' ForceAtlas2
#'
#' An R wrapper around the Python implementation of ForceAtlas2 found at
#' https://github.com/bhargavchippada/forceatlas2
#'
#' Currently a very restricted version
#'
#' @param init_pos Array of initial positions, i.e. PCA or PAGA embeddings
#' @param adjacencies A graph in 2D format
#' @param iterations Number of times to iterate the main loop
#' @param ... Additional arguments to pass to the forceatlas2 object instantiation
#'
#' @rdname fa2
#'
fa2 <- function(init_pos,
                adjacencies,
                iterations = 500,
                ...){

  required_modules <- c("fa2","numpy")
  for (i in required_modules) {
    if (!py_module_available(i)) {
      stop(glue("The {i} module is unavailable.
         Please activate the appropriate environment or
         install the module."))
    }
  }

  fa2 <- import(module = "fa2", delay_load = TRUE)
  np <- import(module = "numpy", delay_load = TRUE)

  # shamelessly based from the fantastic Scanpy package
  forceatlas2 <- fa2$ForceAtlas2(outboundAttractionDistribution=FALSE,  # Dissuade hubs
                                 linLogMode=FALSE,  # NOT IMPLEMENTED
                                 adjustSizes=FALSE,  # Prevent overlap (NOT IMPLEMENTED)
                                 edgeWeightInfluence=1.0,
                                 # Performance
                                 jitterTolerance=1.0,  # Tolerance
                                 barnesHutOptimize=TRUE,
                                 barnesHutTheta=1.2,
                                 multiThreaded=FALSE,  # NOT IMPLEMENTED
                                 # Tuning
                                 scalingRatio=2.0,
                                 strongGravityMode=FALSE,
                                 gravity=1.0,
                                 # Log
                                 verbose=TRUE)

  fa_results <- forceatlas2$forceatlas2(G = adjacencies,
                                        pos = init_pos,
                                        iterations = as.integer(iterations))
  fa_embeddings <- np$array(fa_results)
  rownames(fa_embeddings) <- rownames(init_pos)
  colnames(fa_embeddings) <- as.character(glue("fa_{1:ncol(fa_embeddings)}"))
  return(fa_embeddings)
}
