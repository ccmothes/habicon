## code to prepare `DATASET` dataset goes here
library(raster)

corr <- raster('data-raw/current_surface.tif')
suit <- raster('data-raw/suitability_surface.tif')
corr_prior <- raster('data-raw/corr_prior.tif')

usethis::use_data(corr, suit, corr_prior, overwrite = TRUE)
