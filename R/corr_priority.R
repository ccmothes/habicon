#' Corridor priority calculation
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
#'  @param progress Logical, whether or not to show a progress bar.
#'
#' @return A SpatRaster object
#'
#' @export
#'
corr_priority <-
  function(suit,
           suit_bin,
           corr_bin,
           resist,
           min_area = terra::res(suit)[1] * terra::res(suit)[2],
           d,
           progress = TRUE) {

    # all rasters must be able to stack

    out <- terra::setValues(suit, NA)


    # ID patches ----
    landscape_suit <- landscapemetrics::get_patches(suit_bin, class = 1, directions = 4, return_raster = TRUE)[[1]][[1]]

    ## remove patches smaller than minimum area

    ### save area metrics for later
    patch_area <- landscapemetrics::lsm_p_area(landscape_suit) |>
      dplyr::rename(patch = class, area_ha = value)

    patch_remove <- patch_area |>
      dplyr::filter(area_ha < min_area / 10000)

    if (nrow(patch_remove != 0)) {
      landscape_suit[terra::`%in%`(landscape_suit, ls_remove$patch)] <- NA
    }

    ## assign layer name
    names(landscape_suit) <- "patch"


    # ID corridors ----

    ## filter corridor map to matrix (i.e., remove corridor cells overlapping with patches)
    corr_matrix <- matrix_map(suit_bin, corr_bin)

    landscape_corr <-
      landscapemetrics::get_patches(corr_matrix, class = 1, directions = 4, return_raster = TRUE)[[1]][[1]]

    ## assign layer name
    names(landscape_corr) <- "corridor"



    # ID edges ----

    pc <-
      terra::as.data.frame(c(landscape_corr, landscape_suit)) |>
      tidyr::drop_na() |>
      dplyr::distinct()

    ## create df of patch characteristics
    patch_char <-
      terra::zonal(suit, landscape_suit, fun = "mean") |>
      dplyr::rename(quality = names(suit)) |>
      dplyr::inner_join(patch_area, by = "patch") |>
      dplyr::mutate(
        area_sqkm = area_ha * 0.01,
        quality_area = area_sqkm * quality
      ) |>
      dplyr::select(-c(layer, level, id, metric))


    # Create edges object of all patch connections
    edges <- pc |>
      dplyr::group_by(corridor) |>
      dplyr::mutate(patch2 = patch) |>
      tidyr::expand(patch, patch2) |> # get all comb of values
      dplyr::filter(!duplicated(paste0(pmax(patch, patch2), pmin(patch, patch2)), corridor)) |> # filter unique combinations
      dplyr::ungroup() |>
      dplyr::filter(!(patch == patch2)) |>
      dplyr::rename(patch1 = patch) |>
      dplyr::select(patch1, patch2, corridor)


    # calculate pairwise distances between patches

    ## filter to just patches connected by corridors to speed up computation
    rpoly <- terra::app(landscape_suit, function(x) {
      x[!x %in% pc$patch] <- NA
      return(x)
    })


    ## convert to polygons and calculate distances
    rpoly <- terra::as.polygons(rpoly)

    dis <- terra::distance(rpoly, pairs = TRUE) |>
      tidyr::as_tibble() |>
      ### change ID to their layer (patch) ID
      dplyr::mutate(
        from = rpoly$lyr.1[from],
        to = rpoly$lyr.1[to]
      )

    ## add distances to edges df
    edges <- edges |>
      dplyr::left_join(dis, by = c("patch1" = "from", "patch2" = "to")) |>
      dplyr::rename(distance = value)

    ## get all unique patches (nodes)
    nodes <- unique(c(edges$patch1, edges$patch2))


    # create graph object
    patch_network <- igraph::graph_from_data_frame(d = edges,
                                                   vertices = nodes,
                                                   directed = F)

    # weighted betweenness ----
    igraph::V(patch_network)$betweenness <- igraph::betweenness(
      patch_network,
      directed = FALSE,
      weights = (igraph::E(patch_network)$distance) + 1 # add one because function doesn't like 0's
    )

    centrality <-
      tidyr::tibble(
        patch = as.numeric(names(igraph::V(patch_network))),
        betweenness = as.numeric(igraph::V(patch_network)$betweenness)
      )

    patch_char <-
      dplyr::left_join(patch_char, centrality, by = "patch")


    # EC for each corridor pixel ----

    # ID patch edges ----
    patch_edge <-
      landscapemetrics::get_boundaries(landscape_suit,
                                       consider_boundary = TRUE,
                                       return_raster = TRUE
      )[[1]] |>
      terra::classify(matrix(c(-Inf, 0, NA), ncol = 3)) |>
      # associate edge cell values with their patch ID
      terra::mask(landscape_suit, mask = _) |>
      # filter only edge cells connected to a corridor
      terra::mask(terra::app(landscape_corr, function(x) {
        x[!x %in% edges$corridor] <- NA
        return(x)
      })) |>
      terra::as.data.frame(cells = TRUE, xy = TRUE)

    ## create conductance matrix for shortest path calc and mask to connecting corriors

    ## keep full raster, masking created errors
    #conductance_matrix <- 1 / resist

    # try masking to matrix with 2 cell buffer
    conductance_matrix <- 1 / (terra::mask(resist, source_map(terra::buffer(
      source_map(corr_matrix), terra::res(corr_matrix)[1] * 2
    ))))

    ### for now, in order to use gdistance must convert terra to raster object
    tran <-
      gdistance::transition(raster::raster(conductance_matrix), function(x) { # convert resist to conductance and mask to just corridors
        mean(x)
      }, 8)

    tran <- gdistance::geoCorrection(tran, type = "c")

    # set up new variable to fill
    #edges[, "ECA"] <- NA

    ### Progress bar
    if (progress) {
      cli::cli_progress_bar(
        name = "Calculating ECA",
        type = "iterator",
        format = "{cli::pb_name} {cli::pb_bar} {cli::pb_percent} | \\
              ETA: {cli::pb_eta} - {cli::pb_elapsed_clock}",
        total = nrow(edges)
      )
    }


    # for all patch connections in edges
    rastPaths <- vector("list", length = nrow(edges))

    for (i in 1:nrow(edges)) {
      # get edge points for patch i and patch j
      p1 <- as.numeric(edges[i, 1])
      p2 <- as.numeric(edges[i, 2])
      c <- as.numeric(edges[i, 3]) #corridor connecting patches

      patch_i <-
        patch_edge |>
        dplyr::filter(patch == p1) |>
        dplyr::select(x, y)

      patch_j <-
        patch_edge |>
        dplyr::filter(patch == p2) |>
        dplyr::select(x, y)

      min_path_sp1 <-  vector("list", length = nrow(patch_i))

      for (j in (1:nrow(patch_i))) {

        min_path_sp1[[j]] <- tryCatch({
          paths <- gdistance::shortestPath(
            tran,
            origin = as.matrix(patch_i[j, ]),
            goal = as.matrix(patch_j),
            output = "SpatialLines"
          ) |>
            # convert to sf
            sf::st_as_sf()

          paths |>
            # add distance variable
            dplyr::mutate(distance = sf::st_length(paths)) |>
            # keep only shortest path
            dplyr::filter(distance == min(distance, na.rm = TRUE))


        }, error = function(msg) {
          # message(paste("Error for path number:", x, "\nInvalid geometry"))
          return(NULL)
        })



      }
      # now need to do the same for opposite direction patch j to patch i

      min_path_sp2 <- vector("list", length = nrow(patch_j))

      for (k in (1:nrow(patch_j))) {

        min_path_sp2[[k]] <- tryCatch(
          {
            paths2 <- gdistance::shortestPath(
              tran,
              origin = as.matrix(patch_j[k, ]),
              goal = as.matrix(patch_i),
              output = "SpatialLines"
            ) |>
              # convert to sf
              sf::st_as_sf()

            paths2 |>
              # add distance variable
              dplyr::mutate(distance = sf::st_length(paths2)) |>
              # keep only shortest path
              dplyr::filter(distance == min(distance, na.rm = TRUE))


          },
          error = function(msg) {
            # message(paste("Error for path number:", x, "\nInvalid geometry"))
            return(NULL)
          }
        )
      }



      sldf <- dplyr::bind_rows(purrr::compact(min_path_sp1), purrr::compact(min_path_sp2))


      # sldf$CR <- sapply(1:nrow(sldf), function(x) {
      #   sum(values(mask(resist, sldf[x, ])), na.rm = TRUE)
      #})

      # calculate EC values for each path

      Ai <- dplyr::filter(patch_char, patch == p1)$area_sqkm
      Qi <- dplyr::filter(patch_char, patch == p1)$quality
      Aj <- dplyr::filter(patch_char, patch == p2)$area_sqkm
      Qj <- dplyr::filter(patch_char, patch == p2)$quality
      Bi <- dplyr::filter(patch_char, patch == p1)$betweenness
      Bj <- dplyr::filter(patch_char, patch == p2)$betweenness
      sldf <- sldf |>
        dplyr::mutate(EC =  ifelse(nrow(
          sldf) == 0, NA, Ai * Qi * (Bi + 1) * Aj * Qj * (Bj + 1) * exp(distance * (log(0.5) / d))
        ))


      rastPaths[[i]] <-
        # in case sldf does not return any paths
        if (inherits(sldf, "tbl_df")) {
          out
        } else {
          terra::rasterize(
            x = sldf,
            y = out,
            #y = test,
            field = "EC",
            fun = sum
          )

        }

      # Update the progress bar
      if (progress) {
        cli::cli_progress_update()
      }


    }

    out <- terra::sprc(rastPaths) %>% terra::mosaic(fun = sum) %>% terra::app(sqrt)

    return(out)

  }

