## code to prepare `DATASET` dataset goes here
library(terra)

corr <- rast('data-raw/current_surface.tif')
suit <- rast('data-raw/suitability_surface.tif')
resist <- rast("data-raw/resistance_surface.tif")
corridor <- rast("data-raw/corridor.tif")
load("data-raw/patch_priorEx.RData")

# change rasterlayers to terra
patch$qwa <- rast(patch$qwa)
patch$btwn <- rast(patch$btwn)
patch$dECA <- rast(patch$dECA)


usethis::use_data(corr, suit, resist, patch, corridor, overwrite = TRUE)
