
#testing
patch <- landscape_suit
conductance <- 1/resist
origin <- as.matrix(patch_i[1, ])
goal <-  as.matrix(patch_j)


shortest_path <- function(patch, conductance, origin, goal){

  cs_rast <- terra::rast(ext(conductance), resolution = res(conductance), crs = crs(conductance))

  from_cell <- terra::cellFromXY(cs_rast, xy = origin)
  to_cell <- terra::cellFromXY(cs_rast, xy = goal)

  # create adjency matrix from conductance layer
  #cs_matrix <- matrix(data = 0, nrow = terra::ncell(r), ncol = terra::ncell(r))
  cs_matrix <- prioritizr::connectivity_matrix(patch, conductance)

  cm_graph <- igraph::graph_from_adjacency_matrix(cs_matrix, mode = "directed", weighted = TRUE)

  igraph::E(cm_graph)$weight <- (1/igraph::E(cm_graph)$weight)

  lcp_graph <- igraph::shortest_paths(cm_graph, from = from_cell, to = to_cell, mode = "out", algorithm = "dijkstra")

  lcps <- lapply(lcp_graph$vpath, FUN = function(i) {

    lcp_xy <- terra::xyFromCell(cs_rast, as.integer(i))
    lcp <- sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(lcp_xy)), crs = crs(conductance))
    return(lcp)
  }
  )

  lcps <- bind_rows(lcps)

}


