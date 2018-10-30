
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
