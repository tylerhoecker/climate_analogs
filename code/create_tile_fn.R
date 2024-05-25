#' @description Divide focal climate data into tiles, and then data frames (here, focal is future)
#' @param tile adsf 
#' @param tile_grid asfdf
#' @param climate_dir asdf
#' @param annual 
#' @param write_out asdf
#' @return asfds

create_tile <- function(
  tile,
  tile_grid = tile_grid,
  climate_dir, 
  annual = TRUE,
  write_out = FALSE
){

    print(paste0("Beginning ", tile))

    # Define tile in tile_grid
    tile_i  <- tile_grid[[tile]] #198162\

    # Test if empty
    test <- terra::extract(
      terra::rast(
        list.files(climate_dir, full.names = TRUE)[[1]]
      ), 
      terra::vect(tile_i)
    )
    if(all(is.na(test[,-1]))){
      return(NULL)
    }

  # Read in and crop to tile all years of climate data 
  data_rast <- list.files(climate_dir, full.names = TRUE) |>
    purrr::map(\(x){
      terra::crop(terra::rast(x), tile_i) |>
        terra::mask(ecoregs)
    } 
    )

  if(annual != TRUE){
    data_rast <- mean(data_rast)
  }
  # Convert to data table
  data_datatable <- lapply(data_rast, FUN = as.data.table, xy = TRUE) # map doesn't work

  # Write out if specified
  if(write_out != FALSE){
    saveRDS(data_datatable, paste0(write_out, "tile_", tile, ".Rds"))
  }

  return(data_datatable)
}
