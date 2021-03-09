#' Calculate habitat patch priority
#'
#' @param suit A RasterLayer object of habitat suitability values, values range from 0-1
#' @param suit_bin A RasterLayer object of binary suitable (1) and unsuitable (0) habitat created using \code{\link[habicon]{bin_map}}
#' @param corr_bin A RasterLayer object of binary corridor (1) and non-corridor (0) values created using \code{\link[habicon]{bin_map}}
#' @param resist A RasterLayer object where values represent resistance to movement, higher values represent greater resistance
#'
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples

patch_priority <- function(suit, suit_bin, corr_bin, resist, min_area = res(suit), d){
  out <- list()
  #make empty rasters to fill results with
  out$qwa <- setValues(suit, NA)
  out$btwn <- setValues(suit, NA)
  out$dECA <- setValues(suit, NA)



  #id patches
  landscape_suit <-
    landscapemetrics::get_patches(suit_bin, class = 1, return_raster = TRUE)$'1'

  #remove small patches (default min area is resolution of rasters)
  ls_remove <- landscapemetrics::lsm_p_area(landscape_suit) %>% filter(value <= (min_area[1]*min_area[2])/10000)

  #landscape_suit[landscape_suit %in% ls_remove$id] <- NA
  landscape_suit[raster::`%in%`(landscape_suit, ls_remove$id)] <- NA


  #ID unique corridors

  #filter corridor map to matrix
  corr_matrix <-  matrix_map(suit_bin, corr_bin)

  landscape_corr <-
    landscapemetrics::get_patches(corr_matrix, class = 1, return_raster = TRUE)$'1'


  #ID connections
  pc <-
    tabularaster::as_tibble(stack(landscape_corr, landscape_suit)) %>% dplyr::filter(cellvalue != 0 | !is.na(cellvalue)) %>%
    dplyr::filter(duplicated(cellindex) |
                    duplicated(cellindex, fromLast = TRUE))  %>% #keep both duplicated values (i.e. patch and corridor)
    tidyr::spread(dimindex, cellvalue) %>% dplyr::rename(corridor = '1', patch = '2') %>% dplyr::distinct(corridor, patch)


  ## get patch area and quality associate with patch ID
  patch_area <-
    landscapemetrics::lsm_p_area(landscape_suit) %>%   #units are hectares?
    dplyr::select(class, value) %>% rename(patch = class, area_ha = value)

  patch_char <-
    zonal(suit, landscape_suit, fun = 'mean') %>% as_tibble() %>%
    rename(quality = mean, patch = zone) %>%
    left_join(patch_area, by = 'patch') %>%
    mutate(area_sqkm = area_ha * 0.01) %>%
    mutate(quality_area = area_sqkm * quality)


  #Create edges object of all patch connections
  edges <- pc %>% dplyr::group_by(corridor) %>%
    dplyr::mutate(patch2 = patch) %>%
    tidyr::expand(patch, patch2) %>%  #get all comb of values
    dplyr::filter(!duplicated(paste0(pmax(patch, patch2), pmin(patch, patch2)), corridor)) %>% #filter unique combinations
    dplyr::ungroup() %>%
    dplyr::filter(!(patch == patch2)) %>%
    dplyr::rename(patch1 = patch)%>%
    select(patch1, patch2, corridor)


  # calculate pairwise distances between patches
  ## filter to just patches connected by corridors to speed up computation
  rpoly <- calc(landscape_suit, function(x) {
    x[!x %in% pc$patch] <- NA
    return(x)
  })
  rpoly <- rasterToPolygons(rpoly, dissolve=T)

  #rpoly <- subset(rpoly, rpoly$layer %in% pc$patch)

  dis <- rgeos::gDistance(rpoly, byid = TRUE) # these rows and headers are row ID not patch ID
  # change row and column names to their layer (patch) ID

  rownames(dis) <- rpoly$layer
  colnames(dis) <- rpoly$layer

  edges[, "distance"] <- NA
  #now add distances to edges df
  for (i in 1:nrow(edges)){
    edges$distance[i] <- dis[as.character(edges$patch1[i]), as.character(edges$patch2[i])]

  }

  nodes <- unique(c(edges$patch1, edges$patch2))

  #create graph object

  patchNetwork <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)


  #weighted betweenness
  V(patchNetwork)$betweenness <- betweenness(patchNetwork,  directed = FALSE, weights = E(patchNetwork)$distance)

  centrality <-
    data.frame(patch = as.numeric(names(V(patchNetwork))),
               betweenness = as.numeric(V(patchNetwork)$betweenness))

  patch_char <-
    left_join(patch_char, centrality, by = "patch")




  #ECA for entire network



  ##for all patch connections (edges), sum( AiQiAjQj*exp(dij*(log(0.5)/D50)) )
  ##lets say D50 is 5km (convert area to km then)

  patch_edge <-
    landscapemetrics::get_boundaries(landscape_suit,
                                     consider_boundary = TRUE,
                                     return_raster = TRUE)[[1]] %>%
    reclassify(matrix(c(-Inf, 0, NA), ncol = 3)) %>%
    mask(landscape_suit, .) #associate edge cell values with their patch ID


  tran <-
    transition(1/mask(resist, corr_matrix), function(x) #convert resist to conductance and mask to just corridors
      mean(x), 8)
  tran <- geoCorrection(tran, type = "c")

  theta <- log(0.5)/d #meters

  edges[, "ECA"] <- NA
  for (i in 1:nrow(edges)){
    # get edge points for patch i and patch j
    p1 <- as.numeric(edges[i,1])
    p2 <- as.numeric(edges[i,2])
    c <- as.numeric(edges[i,3]) #corridor connecting patches
    patch_i <-
      mask(patch_edge, calc(landscape_corr, function (x) {
        x[x != c] <- NA
        return(x)
      })) %>% tabularaster::as_tibble() %>% dplyr::filter(cellvalue == p1) %>%
      dplyr::pull(cellindex) %>% xyFromCell(patch_edge, .)

    patch_j <-
      mask(patch_edge, calc(landscape_corr, function (x) {
        x[x != c] <- NA
        return(x)
      })) %>% tabularaster::as_tibble() %>% dplyr::filter(cellvalue == p2) %>%
      dplyr::pull(cellindex) %>% xyFromCell(patch_edge, .)

    spatdf <- vector("list", length = nrow(patch_i))

    for (j in (1:nrow(patch_i))) {
      paths <- shortestPath(tran,
                            origin = patch_i[j, ],
                            goal = patch_j,
                            output = "SpatialLines")

      spatdf[[j]] <-  data.frame(distance = sapply(1:length(paths), function(x) {
        tryCatch({
          rgeos::gLength(paths[x,])
        }, error = function(msg) {
          #message(paste("Error for path number:", x, "\nInvalid geometry"))
          return(NA)
        })
      }))

    }

    spatdf <- do.call("rbind", spatdf)


    dij <- min(spatdf$distance, na.rm = TRUE) # meters
    Ai <- as.numeric(patch_char[patch_char$patch == p1, 4])
    Qi <- as.numeric(patch_char[patch_char$patch == p1, 2])
    Aj <- as.numeric(patch_char[patch_char$patch == p2, 4])
    Qj <- as.numeric(patch_char[patch_char$patch == p2, 2])
    edges$ECA[i] <- Ai*Qi*Aj*Qj*exp(dij*theta)

  }
  ECA <- sqrt(sum(edges$ECA))

  # now calculate ECA for each patch and add to patch char
  #ECAi <- vector("list", length = length(nodes))

  patch_char[, "dECA"] <- NA
  for (i in 1:length(nodes)){

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
