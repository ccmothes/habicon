## code to prepare `DATASET` dataset goes here
library(terra)

corr <- wrap(rast('data-raw/current_surface.tif'))
suit <- wrap(rast('data-raw/suitability_surface.tif'))
resist <- wrap(rast("data-raw/resistance_surface.tif"))
corridor <- wrap(rast("data-raw/corridor.tif"))
load("data-raw/patch_priorEx.RData")

# change rasterlayers to terra
patch$qwa <- wrap(rast(patch$qwa))
patch$btwn <- wrap(rast(patch$btwn))
patch$dECA <- wrap(rast(patch$dECA))


usethis::use_data(corr, suit, resist, patch, corridor, overwrite = TRUE)
