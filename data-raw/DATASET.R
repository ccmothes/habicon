## code to prepare `DATASET` dataset goes here
library(raster)

corr <- raster('data-raw/current_surface.tif')
suit <- raster('data-raw/suitability_surface.tif')
resist <- raster("data-raw/resistance_surface.tif")
patch <- load("data-raw/patch_priorEx.RData")
corridor <- raster("data-raw/corridor.tif")

usethis::use_data(corr, suit, resist, patch, corridor, overwrite = TRUE)
