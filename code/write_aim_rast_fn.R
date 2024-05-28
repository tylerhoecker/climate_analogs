#' @description Convert datatables/frames to rasters and write out
#' @param data The datatable/frame 
#' @param band Which column of the raster/band to write?
#' @param tile Tile ID
#' @return Nothing - to save RAM?

write_aim_rast <- function(data, tile, band) { 
    #print(paste0('Rasterizing ',band,'...'))
    r <- data |> 
      select(x, y, all_of(band))  |> 
      terra::rast(x = _, 
                  type = 'xyz', 
                  crs = crs(lf_bps_220)
                 ) 
        
    writeRaster(r, 
                paste0(out_dir, "/tile_", tile, "_", band, '.tif'), 
                datatype = "INT2S",
                gdal = c("PROJECTION=EPSG:3857",
                         "TILED=YES",
                         "BLOCKXSIZE=128",
                         "BLOCKYSIZE=128",
                         "OVERVIEW-RESAMPLING=NEAREST",
                         "COMPRESS=DEFLATE"),
                overwrite = TRUE)
  }

