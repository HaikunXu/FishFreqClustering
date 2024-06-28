#'
#' \code{packbylatlon} This function ...
#'
#' @export

packbylatlon <- function(datamat) {
  lat = datamat$lat
  lon = datamat$lon
  lats = seq(min(lat), max(lat), 5) #9
  lons = seq(min(lon), max(lon), 5)#13
  nlat = length(lats) #9
  nlon = length(lons) #13
  nrowsize = nlat * nlon
  bggfunc3 = matrix(0, nrowsize, ncol = 4 + 202 + 1)
  bggfunc3 = data.frame(bggfunc3)
  nrowdatamat = nrow(datamat)
  ctable = matrix(0, nlat, nlon)
  for (i in 1:nrowdatamat) {
    nn = (lat[i] - min(lat)) / 5 * nlon + (lon[i] - min(lon)) / 5 + 1
    bggfunc3[nn, 5:206] <- bggfunc3[nn, 5:206] + datamat[i, 5:206]
    ctable[(lat[i] - min(lat)) / 5 + 1, (lon[i] - min(lon)) / 5 + 1] =  
      ctable[(lat[i] - min(lat)) / 5 + 1, (lon[i] - min(lon)) / 5 + 1] + 1
  }
  bggfunc3[, 1] = 2003
  bggfunc3[, 2] = 1
  bggfunc3[, 3] = rep(lats, each = length(lons))
  bggfunc3[, 4] = rep(lons, length(lats))
  bggfunc3[, 207] = as.vector(t(ctable))
  return(list(table1 = bggfunc3, table2 = ctable))
}