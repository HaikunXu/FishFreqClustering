#'
#' \code{packbylatlon} This function ...
#'
#' @export

packbylatlon <- function(datamat, dlat, dlon, nbins) {
  lat = datamat$lat
  lon = datamat$lon
  lats = seq(min(lat), max(lat), dlat) # 9
  lons = seq(min(lon), max(lon), dlon) # 13
  nlat = length(lats) # 9
  nlon = length(lons) # 13
  nrowsize = nlat * nlon
  bggfunc3 = matrix(0, nrowsize, ncol = nbins + 4)
  bggfunc3 = data.frame(bggfunc3)
  ctable = matrix(0, nlat, nlon)
  for (i in 1:nrow(datamat)) {
    nn = (lat[i] - min(lat)) / dlat * nlon + (lon[i] - min(lon)) / dlon + 1
    bggfunc3[nn, 4:(nbins + 3)] <- bggfunc3[nn, 4:(nbins + 3)] + datamat[i, 4:(nbins + 3)]
    ctable[(lat[i] - min(lat)) / dlat + 1, (lon[i] - min(lon)) / dlat + 1] =  
      ctable[(lat[i] - min(lat)) / dlat + 1, (lon[i] - min(lon)) / dlat + 1] + 1
  }

  bggfunc3[, 1] = paste0(min(datamat[1]),"-",max(datamat[1]))
  bggfunc3[, 2] = rep(lats, each = length(lons))
  bggfunc3[, 3] = rep(lons, length(lats))
  bggfunc3[, nbins + 4] = as.vector(t(ctable))
  return(list(table1 = bggfunc3, table2 = ctable))
}