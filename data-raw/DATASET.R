## code to prepare `DATASET` dataset goes here
library(Seurat)
library(tidyverse)

pbmc3k <- ReadH5AD("data-raw/pbmc3k_raw.h5ad") %>%
  NormalizeData() %>%
  PercentageFeatureSet(object = .,
                       pattern = "^MT-",
                       col.name = "percent_mt") %>%
  SCTransform(object = .,
              vars.to.regress = "percent_mt") %>%
  RunPCA()

usethis::use_data("DATASET")
