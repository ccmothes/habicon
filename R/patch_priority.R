#' Calculate habitat patch priority
#'
#' @param suit A SpatRaster object of continuous habitat suitability values, where values range
#'   from 0-1. See \code{\link{rescale_map}} for rescaling values to 0-1.
#' @param suit_bin A SpatRaster object of binary suitable (1) and unsuitable
#'   (0) habitat created using \code{\link{bin_map}}
#' @param corr_bin A SpatRaster object of binary corridor (1) and non-corridor
#'   (0) values created using \code{\link{bin_map}}
#' @param resist A SpatRaster object where values represent resistance to
#'   movement, higher values represent greater resistance.
#'  @param min_area The minimum area to be considered a habitat patch. Value in square meters. Default
#'    is to remove any patches that are only one pixel in size.
#'  @param d Minimum dispersal distance (in meters).
#'
#' @return A list object with one SpatRaster for each of the three connectivity
#'   metrics and a summary table.
#'
#' @export
patch_priority <- function(suit, suit_bin, corr_bin, resist, min_area = res(suit)[1] * res(suit)[2], d) {

  out <- list()
  # make empty rasters to fill results with
  out$qwa <- terra::setValues(suit, NA)
  out$btwn <- terra::setValues(suit, NA)
  out$dECA <- terra::setValues(suit, NA)


  # ID patches
  landscape_suit <- landscapemetrics::get_patches(suit_bin, class = 1, return_raster = TRUE)[[1]][[1]]

  ## remove patches smaller than minimum area

  ### save area metrics for later
  patch_area <- landscapemetrics::lsm_p_area(landscape_suit) %>%
    rename(patch = class, area_ha = value)

  patch_remove <- patch_area %>%
    filter(area_ha < min_area / 10000)

  if(nrow(patch_remove != 0)) {
    landscape_suit[terra::`%in%`(landscape_suit, ls_remove$patch)] <- NA
  }

  ## assign layer name
  names(landscape_suit) <- "patch"


  # ID unique corridors

  ## filter corridor map to matrix (i.e., remove corridor cells overlapping with patches)
  corr_matrix <- matrix_map(suit_bin, corr_bin)

  landscape_corr <-
    landscapemetrics::get_patches(corr_matrix, class = 1, return_raster = TRUE)[[1]][[1]]

  ## assign layer name
  names(landscape_corr) <- "corridor"


  # ID connections

  pc <-
    tabularaster::as_tibble(c(landscape_corr, landscape_suit)) %>%
    tidyr::drop_na() %>%
    distinct()

  ## create df of patch characteristics
  patch_char <-
    terra::zonal(suit, landscape_suit, fun = "mean") %>%
    dplyr::rename(quality = names(suit)) %>%
    dplyr::inner_join(patch_area, by = "patch") %>%
    dplyr::mutate(area_sqkm = area_ha * 0.01,
           quality_area = area_sqkm * quality) %>%
    dplyr::select(-c(layer, level, id, metric))


  # Create edges object of all patch connections
  edges <- pc %>%
    dplyr::group_by(corridor) %>%
    dplyr::mutate(patch2 = patch) %>%
    tidyr::expand(patch, patch2) %>% # get all comb of values
    dplyr::filter(!duplicated(paste0(pmax(patch, patch2), pmin(patch, patch2)), corridor)) %>% # filter unique combinations
    dplyr::ungroup() %>%
    dplyr::filter(!(patch == patch2)) %>%
    dplyr::rename(patch1 = patch) %>%
    dplyr::select(patch1, patch2, corridor)


  # calculate pairwise distances between patches

  ## filter to just patches connected by corridors to speed up computation
  rpoly <- terra::app(landscape_suit, function(x) {
    x[!x %in% pc$patch] <- NA
    return(x)
  })

  ## convert to polygons and calculate distances
  rpoly <- terra::as.polygons(rpoly)

  dis <- terra::distance(rpoly, pairs = TRUE) %>%
    as_tibble() %>%
  ### change ID to their layer (patch) ID
    mutate(from = rpoly$lyr.1[from],
           to = rpoly$lyr.1[to])

  ## add distances to edges df
  edges <- edges %>%
    left_join(dis, by = c("patch1" = "from", "patch2" = "to"))

  ## get all unique patches (nodes)
  nodes <- unique(c(edges$patch1, edges$patch2))


  # create graph object
  patch_network <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = F)


  # weighted betweenness
  igraph::V(patch_network)$betweenness <- igraph::betweenness(patch_network, directed = FALSE, weights = igraph::E(patch_network)$distance)

  centrality <-
    tibble(
      patch = as.numeric(names(igraph::V(patch_network))),
      betweenness = as.numeric(igraph::V(patch_network)$betweenness)
    )

  patch_char <-
    left_join(patch_char, centrality, by = "patch")




  # ECA for entire network


  ## for all patch connections (edges), sum( AiQiAjQj*exp(dij*(log(0.5)/D50)) ) where D50 is given dispersal distance

  patch_edge <-
    landscapemetrics::get_boundaries(landscape_suit,
      consider_boundary = TRUE,
      return_raster = TRUE
    )[[1]] %>%
    terra::classify(matrix(c(-Inf, 0, NA), ncol = 3)) %>%
    terra::mask(landscape_suit, .) # associate edge cell values with their patch ID


  tran <-
    transition(1 / mask(resist, corr_matrix), function(x) { # convert resist to conductance and mask to just corridors
      mean(x)
    }, 8)

  tran <- geoCorrection(tran, type = "c")

  theta <- log(0.5) / d # meters

  edges[, "ECA"] <- NA
  for (i in 1:nrow(edges)) {
    # get edge points for patch i and patch j
    p1 <- as.numeric(edges[i, 1])
    p2 <- as.numeric(edges[i, 2])
    c <- as.numeric(edges[i, 3]) # corridor connecting patches
    patch_i <-
      mask(patch_edge, calc(landscape_corr, function(x) {
        x[x != c] <- NA
        return(x)
      })) %>%
      tabularaster::as_tibble() %>%
      dplyr::filter(cellvalue == p1) %>%
      dplyr::pull(cellindex) %>%
      xyFromCell(patch_edge, .)

    patch_j <-
      mask(patch_edge, calc(landscape_corr, function(x) {
        x[x != c] <- NA
        return(x)
      })) %>%
      tabularaster::as_tibble() %>%
      dplyr::filter(cellvalue == p2) %>%
      dplyr::pull(cellindex) %>%
      xyFromCell(patch_edge, .)

    spatdf <- vector("list", length = nrow(patch_i))

    for (j in (1:nrow(patch_i))) {
      paths <- shortestPath(tran,
        origin = patch_i[j, ],
        goal = patch_j,
        output = "SpatialLines"
      )

      spatdf[[j]] <- data.frame(distance = sapply(1:length(paths), function(x) {
        tryCatch(
          {
            rgeos::gLength(paths[x, ])
          },
          error = function(msg) {
            # message(paste("Error for path number:", x, "\nInvalid geometry"))
            return(NA)
          }
        )
      }))
    }

    spatdf <- do.call("rbind", spatdf)


    dij <- min(spatdf$distance, na.rm = TRUE) # meters
    Ai <- as.numeric(patch_char[patch_char$patch == p1, 4])
    Qi <- as.numeric(patch_char[patch_char$patch == p1, 2])
    Aj <- as.numeric(patch_char[patch_char$patch == p2, 4])
    Qj <- as.numeric(patch_char[patch_char$patch == p2, 2])
    edges$ECA[i] <- Ai * Qi * Aj * Qj * exp(dij * theta)
  }

  ECA <- sqrt(sum(edges$ECA))

  # now calculate ECA for each patch and add to patch char
  # ECAi <- vector("list", length = length(nodes))

  patch_char[, "dECA"] <- NA
  for (i in 1:length(nodes)) {
    ECAi <- filter(edges, patch1 != nodes[i] & patch2 != nodes[i]) %>% pull(ECA)
    ECAi <- sqrt(sum(ECAi))
    dECA <- (ECA - ECAi) / ECA
    patch_char$dECA[patch_char$patch == nodes[i]] <- dECA
  }


  # for each patch cell, assign value that corresponds with patch ID

  # Quality weighted area
  for (i in 1:ncell(landscape_suit)) {
    if (is.na(landscape_suit[i]) | landscape_suit[i] == 0) {
      # keep non-patch cells as NA
      out$qwa[i] <- NA
    } else {
      patch_id <- landscape_suit[i]

      out$qwa[i] <- dplyr::filter(patch_char, patch == patch_id) %>% dplyr::pull(quality_area)
    }
  }
  names(out$qwa) <- "Quality_weighted_area"



  # Weighted betweenness centrality
  for (i in 1:ncell(landscape_suit)) {
    if (is.na(landscape_suit[i]) | landscape_suit[i] == 0) {
      # keep non-patch cells as NA
      out$btwn[i] <- NA
    } else {
      patch_id <- landscape_suit[i]

      out$btwn[i] <- dplyr::filter(patch_char, patch == patch_id) %>% dplyr::pull(betweenness)
    }
  }
  names(out$btwn) <- "Weighted_betweenness_centrality"

  # dECA
  for (i in 1:ncell(landscape_suit)) {
    if (is.na(landscape_suit[i]) | landscape_suit[i] == 0) {
      # keep non-patch cells as NA
      out$dECA[i] <- NA
    } else {
      patch_id <- landscape_suit[i]

      out$dECA[i] <- dplyr::filter(patch_char, patch == patch_id) %>% dplyr::pull(dECA)
    }
  }
  names(out$dECA) <- "dECA"


  out$summary_table <- patch_char


  return(out)
}
