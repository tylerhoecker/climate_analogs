#' @description Divide focal climate data into tiles, and then data frames (here, focal is future)
#' @param tile adsf 
#' @param climate_dir asdf
#' @param annual 
#' @param write_out asdf
#' @return asfds

create_tile <- function(
  tile,
  climate_dir, 
  annual,
  write_out = FALSE
){
  # Read in and crop to tile all years of climate data 
  data_rast <- list.files(climate_dir, full.names = TRUE) |>
    purrr::map(\(x) terra::crop(terra::rast(x), tile))

  if(annual != TRUE){
    data_rast <- mean(data_rast)
  }
  # Convert to data table
  data_datatable <- lapply(data_rast, FUN = as.data.table, xy = TRUE) # map doesn't work

  # Write out if specified
  if(write_out == TRUE){
    saveRDS(data_datatable, paste0(write_out, "tile_datatable"))
  }

  return(data_datatable)
}
