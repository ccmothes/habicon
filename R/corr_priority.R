#' Corridor priority calculation
#'
#' @param suit A RasterLayer object of habitat suitability values, values range from 0-1
#' @param suit_bin A RasterLayer object of binary suitable (1) and unsuitable (0) habitat created using \code{\link[habicon]{bin_map}}
#' @param corr_bin A RasterLayer object of binary corridor (1) and non-corridor (0) values created using \code{\link[habicon]{bin_map}}
#'#' @param resist A RasterLayer object where values represent resistance to movement, higher values represent greater resistance

#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
corr_priorityNew <-
  function(suit,
           suit_bin,
           corr_bin,
           resist,
           min_area = res(suit),
           d) {
    #all rasters must be able to stack

    out <- setValues(suit, NA)


    #id patches
    landscape_suit <-
      landscapemetrics::get_patches(suit_bin, class = 1, return_raster = TRUE)$'1'

    #remove small patches (default min area is resolution of rasters)
    ls_remove <-
      landscapemetrics::lsm_p_area(landscape_suit) %>% filter(value == (min_area[1] *
                                                                          min_area[2]) / 10000)

    landscape_suit[landscape_suit %in% ls_remove$id] <- NA

    #ID unique corridors

    #filter corridor map to matrix
    corr_matrix <-  matrix_map(suit_bin, corr_bin)

    landscape_corr <-
      landscapemetrics::get_patches(corr_matrix, class = 1, return_raster = TRUE)$'1'



    #use package tabluraraster to ID overlapping corr/suit patches

    pc <-
      tabularaster::as_tibble(stack(landscape_corr, landscape_suit)) %>% filter(cellvalue != 0 |
                                                                                  !is.na(cellvalue)) %>%
      filter(duplicated(cellindex) |
               duplicated(cellindex, fromLast = TRUE))  %>% #keep both duplicated values (i.e. patch and corridor)
      spread(dimindex, cellvalue) %>% rename(corridor = '1', patch = '2') %>% distinct(corridor, patch)



    ## get patch area and quality associate with patch ID
    patch_area <-
      landscapemetrics::lsm_p_area(landscape_suit) %>%   #units are hectares?
      dplyr::select(class, value) %>% rename(patch = class, area_ha = value)

    patch_char <-
      zonal(suit, landscape_suit, fun = 'mean') %>% as_tibble() %>%
      rename(quality = mean, patch = zone) %>%
      left_join(patch_area, by = 'patch') %>%
      mutate(area_sqkm = area_ha * 0.01)

    patch_edge <-
      landscapemetrics::get_boundaries(landscape_suit,
                                       consider_boundary = TRUE,
                                       return_raster = TRUE)[[1]] %>%
      reclassify(matrix(c(-Inf, 0, NA), ncol = 3)) %>%
      mask(landscape_suit, .) #associate edge cell values with their patch ID

    #calculations for patch centrality
    #Create edges object of all patch connections
    edges <- pc %>% dplyr::group_by(corridor) %>%
      dplyr::mutate(patch2 = patch) %>%
      tidyr::expand(patch, patch2) %>%  #get all comb of values
      dplyr::filter(!duplicated(paste0(
        pmax(patch, patch2), pmin(patch, patch2)
      ), corridor)) %>% #filter unique combinations
      dplyr::ungroup() %>%
      dplyr::filter(!(patch == patch2)) %>%
      dplyr::rename(patch1 = patch) %>%
      select(patch1, patch2, corridor)


    # calculate pairwise distances between patches
    ## filter to just patches connected by corridors to speed up computation
    rpoly <- calc(landscape_suit, function(x) {
      x[!x %in% pc$patch] <- NA
      return(x)
    })
    rpoly <- rasterToPolygons(rpoly, dissolve = T)

    #rpoly <- subset(rpoly, rpoly$layer %in% pc$patch)

    dis <-
      rgeos::gDistance(rpoly, byid = TRUE) # these rows and headers are row ID not patch ID
    # change row and column names to their layer (patch) ID

    rownames(dis) <- rpoly$layer
    colnames(dis) <- rpoly$layer

    edges[, "distance"] <- NA
    #now add distances to edges df
    for (i in 1:nrow(edges)) {
      edges$distance[i] <-
        dis[as.character(edges$patch1[i]), as.character(edges$patch2[i])]

    }

    nodes <- unique(c(edges$patch1, edges$patch2))

    #create graph object

    patchNetwork <-
      graph_from_data_frame(d = edges,
                            vertices = nodes,
                            directed = F)


    #weighted betweenness
    V(patchNetwork)$betweenness <-
      betweenness(patchNetwork,
                  directed = FALSE,
                  weights = E(patchNetwork)$distance)

    centrality <-
      data.frame(patch = as.numeric(names(V(patchNetwork))),
                 betweenness = as.numeric(V(patchNetwork)$betweenness))

    patch_char <-
      left_join(patch_char, centrality, by = "patch")

    # create transition layer for costDistance calc from point to patch
    #corr_bin needs 0 to be NA for the mask
    #corr_bin <- source_map(corr_bin) #adding this made computers crash

    tran <-
      transition(1 / mask(resist, corr_matrix), function(x)
        #convert resist to conductance and mask to just corridors
        mean(x), 8)
    #tran <- transition(pall_current, function(x) mean(x), 8)
    tran <- geoCorrection(tran, type = "c")

    theta <- log(0.5) / d

    # for all patch connections in edges
    rastPaths <- vector("list", length = nrow(edges))

    for (i in 1:nrow(edges)) {
      # get edge points for patch i and patch j
      p1 <- as.numeric(edges[i, 1])
      p2 <- as.numeric(edges[i, 2])
      c <- as.numeric(edges[i, 3]) #corridor connecting patches
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

      spatdf1 <- vector("list", length = nrow(patch_i))
      sldf1 <- vector("list", length = nrow(patch_i))

      for (j in (1:nrow(patch_i))) {
        paths <- shortestPath(tran,
                              origin = patch_i[j, ],
                              goal = patch_j,
                              output = "SpatialLines")

        spatdf1[[j]] <-
          data.frame(distance = sapply(1:length(paths), function(x) {
            tryCatch({
              rgeos::gLength(paths[x,])
            }, error = function(msg) {
              #message(paste("Error for path number:", x, "\nInvalid geometry"))
              return(NA)
            })
          }))

        rownames(spatdf1[[j]]) <-
          sapply(1:length(paths), function(x) {
            paths@lines[[x]]@ID
          })

        sldf1[[j]] <-
          SpatialLinesDataFrame(paths, data = spatdf1[[j]])

        sldf1[[j]] <-
          subset(sldf1[[j]], distance == min(distance))

        #but now need to do for opposite direction patch j to patch i...

      }

      #sldf1 <- do.call("rbind", sldf1)

      spatdf2 <- vector("list", length = nrow(patch_j))
      sldf2 <- vector("list", length = nrow(patch_j))

      for (k in (1:nrow(patch_j))) {
        paths <- shortestPath(tran,
                              origin = patch_j[k, ],
                              goal = patch_i,
                              output = "SpatialLines")

        spatdf2[[k]] <-
          data.frame(distance = sapply(1:length(paths), function(x) {
            tryCatch({
              rgeos::gLength(paths[x,])
            }, error = function(msg) {
              #message(paste("Error for path number:", x, "\nInvalid geometry"))
              return(NA)
            })
          }))

        rownames(spatdf2[[k]]) <-
          sapply(1:length(paths), function(x) {
            paths@lines[[x]]@ID
          })

        sldf2[[k]] <-
          SpatialLinesDataFrame(paths, data = spatdf2[[k]])

        sldf2[[k]] <-
          subset(sldf2[[k]], distance == min(distance))

        #but now need to do for opposite direction patch j to patch i...

      }

      #sldf2 <- do.call("rbind", sldf2)

      sldf <- do.call("rbind", c(sldf1, sldf2))

      # sldf$CR <- sapply(1:nrow(sldf), function(x) {
      #   sum(values(mask(resist, sldf[x, ])), na.rm = TRUE)
      #})




      Ai <- filter(patch_char, patch == p1)$area_sqkm
      Qi <- filter(patch_char, patch == p1)$quality
      Aj <- filter(patch_char, patch == p2)$area_sqkm
      Qj <- filter(patch_char, patch == p2)$quality
      Bi <- filter(patch_char, patch == p1)$betweenness
      Bj <- filter(patch_char, patch == p2)$betweenness
      sldf$EC <- sapply(1:nrow(sldf), function(x) {
        Ai * Qi * (Bi + 1) * Aj * Qj * (Bj + 1) * exp(sldf[x, ]$distance * theta) #* (1 / sldf[x,]$CR)

      }) # what is happening with NA rows? Need to filter those out beforehand?


      rastPaths[[i]] <-
        rasterize(
          x = sldf,
          y = out,
          field = "EC",
          fun = sum
        )



    }
    rastPaths <- do.call(mosaic, c(rastPaths, fun = sum))
    out <- calc(rastPaths, sqrt)

    return(out)
  }

