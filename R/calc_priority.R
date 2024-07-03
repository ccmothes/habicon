#' Habitat Patch and Corridor priority calculations
#'
#' @param type Type of habitat to calculate priority for, either "patch", "corridor", or "all" for both.
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
#' @return A SpatRaster object
#'
#' @export
calc_priority <-
  function(type = "all",
           suit,
           suit_bin,
           corr_bin,
           resist,
           min_area = res(suit)[1] * res(suit)[2],
           d) {

    # create output based on type
    if(type == "all") {
      out <- list()
      # make empty rasters to fill results with
      out$qwa <- terra::setValues(suit, NA)
      out$btwn <- terra::setValues(suit, NA)
      out$dECA <- terra::setValues(suit, NA)
      out$corridor <- terra::setValues(suit, NA)

    }

    if(type == "patch") {
      out <- list()
      # make empty rasters to fill results with
      out$qwa <- terra::setValues(suit, NA)
      out$btwn <- terra::setValues(suit, NA)
      out$dECA <- terra::setValues(suit, NA)

    }

    if(type == "corridor") {
      out <- terra::setValues(suit, NA)

    }

    # error handling with type input
    if(!type %in% c("all", "patch", "corridor")) {
      stop("type must be one of 'patch', 'corridor' or 'all'")
    }


    # ID patches
    landscape_suit <- landscapemetrics::get_patches(suit_bin, class = 1, directions = 4, return_raster = TRUE)[[1]][[1]]

    ## remove patches smaller than minimum area

    ### save area metrics for later
    patch_area <- landscapemetrics::lsm_p_area(landscape_suit) %>%
      rename(patch = class, area_ha = value)

    patch_remove <- patch_area %>%
      filter(area_ha < min_area / 10000)

    if (nrow(patch_remove != 0)) {
      landscape_suit[terra::`%in%`(landscape_suit, ls_remove$patch)] <- NA
    }

    ## assign layer name
    names(landscape_suit) <- "patch"


    # ID unique corridors

    ## filter corridor map to matrix (i.e., remove corridor cells overlapping with patches)
    corr_matrix <- matrix_map(suit_bin, corr_bin)

    landscape_corr <-
      landscapemetrics::get_patches(corr_matrix, class = 1, directions = 4, return_raster = TRUE)[[1]][[1]]

    ## assign layer name
    names(landscape_corr) <- "corridor"


    # ID connections

    pc <-
      terra::as.data.frame(c(landscape_corr, landscape_suit)) %>%
      tidyr::drop_na() %>%
      distinct()

    ## create df of patch characteristics
    patch_char <-
      terra::zonal(suit, landscape_suit, fun = "mean") %>%
      dplyr::rename(quality = names(suit)) %>%
      dplyr::inner_join(patch_area, by = "patch") %>%
      dplyr::mutate(
        area_sqkm = area_ha * 0.01,
        quality_area = area_sqkm * quality
      ) %>%
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
      mutate(
        from = rpoly$lyr.1[from],
        to = rpoly$lyr.1[to]
      )

    ## add distances to edges df
    edges <- edges %>%
      left_join(dis, by = c("patch1" = "from", "patch2" = "to")) %>%
      rename(distance = value)

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

    ### identify patch edges
    patch_edge <-
      landscapemetrics::get_boundaries(landscape_suit,
                                       consider_boundary = TRUE,
                                       return_raster = TRUE
      )[[1]] %>%
      terra::classify(matrix(c(-Inf, 0, NA), ncol = 3)) %>%
      # associate edge cell values with their patch ID
      terra::mask(landscape_suit, .) %>%
      # filter only edge cells connected to a corridor
      mask(., terra::app(landscape_corr, function(x) {
        x[!x %in% edges$corridor] <- NA
        return(x)
      })) %>%
      terra::as.data.frame(cells = TRUE, xy = TRUE)

    ### create conductance matrix for shortest path calc and mask to connecting corridors

    ## keep full raster, masking created errors
    conductance_matrix <- 1 / resist

    ### for now, in order to use gdistance must convert terra to raster object
    tran <-
      gdistance::transition(raster::raster(conductance_matrix), function(x) { # convert resist to conductance and mask to just corridors
        mean(x)
      }, 8)

    tran <- gdistance::geoCorrection(tran, type = "c")

    # set up new variable to fill
    edges[, "ECA"] <- NA

    # This is where EC calculation splits base on type ---------------------------






  }
