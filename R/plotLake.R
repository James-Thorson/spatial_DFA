plotLake <-
function( mapObj, col, site_set, cex=2, ... ){
  # Load package
  require(maptools)

  # Plot map
  plot(1, type="n", xlim=attr(mapObj,"bbox")['x',], ylim=attr(mapObj,"bbox")['y',], ...)
  plot(mapObj, add=TRUE)

  #load the coordinates for each site--as df
  stMercCoord <- read.table( paste0(DataFile,"bsSiteMercCoord.txt"), head=TRUE)
  LatLong2UTM <- match( site_set, table=as.character(stMercCoord[,'site']))

  #plot site locations
  points(stMercCoord$v[LatLong2UTM], stMercCoord$e[LatLong2UTM], pch=20, cex=cex, col=col)
}
