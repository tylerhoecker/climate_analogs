#' Transform raster to data.table
#' 
#' @param x  Raster* object
#' @param row.names	`NULL` or a character vector giving the row names for the data frame. Missing values are not allowed
#' @param optional	logical. If `TRUE`, setting row names and converting column names (to syntactic names: see make.names) is optional
#' @param xy  logical. If `TRUE`, also return the spatial coordinates
#' @param centroids	logical. If TRUE return the centroids instead of all spatial coordinates (only relevant if xy=TRUE)
#' @param sepNA	logical. If TRUE the parts of the spatial objects are separated by lines that are NA (only if xy=TRUE and, for polygons, if centroids=FALSE
#' @param ...	 Additional arguments (none) passed to `raster::as.data.frame`
#' 
#' @value returns a data.table object
#' @examples
#' logo <- brick(system.file("external/rlogo.grd", package="raster"))
#' v <- as.data.table(logo)
#' @import

as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = canProcessInMemory(x, 2), ...) {
  stopifnot(require("data.table"))
  if(inmem) {
    v <- as.data.table(as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
  } else {
    tr <- blockSize(x, n=2)
    l <- lapply(1:tr$n, function(i) 
      as.data.table(as.data.frame(getValues(x, 
                                            row=tr$row[i], 
                                            nrows=tr$nrows[i]), 
                                  row.names=row.names, optional=optional, xy=xy, ...)))
    v <- rbindlist(l)
  }
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  setnames(v, coln)
  v
}

#' @param xy logical. If TRUE, the coordinates of each raster cell are included
#' @param cells logical. If TRUE, the cell numbers of each raster cell are included
#' @param na.rm logical. If TRUE, cells that have a NA value in at least one layer are removed
#' @param ...	 Additional arguments (none) passed to `terra::as.data.frame`
#' @value returns a data.table object
#' @examples 
#' r <- rast(ncols=2, nrows=2)
#' values(r) <- 1:ncell(r)
#' as.data.table(r, xy = TRUE)
#' @importFrom terra as.data.frame
#' @importFrom data.table as.data.table
as.data.table.SpatRaster <- function(x, optional = FALSE, xy = FALSE, ...) {
  stopifnot(require("data.table"))
  v <- as.data.table(as.data.frame(x, optional = optional, xy = xy, ...))
  
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  setnames(v, coln)
  v
}

if (!isGeneric("as.data.table")) {
  setGeneric("as.data.table", function(x, ...)
    standardGeneric("as.data.table"))
}  

setMethod('as.data.table', signature(x='data.frame'), data.table::as.data.table)
# make sure you have terra or raster loaded (as needed) before
setMethod('as.data.table', signature(x='Raster'), as.data.table.raster)
setMethod('as.data.table', signature(x='SpatRaster'), as.data.table.SpatRaster)
