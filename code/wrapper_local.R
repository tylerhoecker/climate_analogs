# Packages and functions -------------------------------------------------------
library(terra)
library(sf)
library(tidyr)
library(dplyr)
library(furrr)
library(purrr)
library(data.table)
source("code/as.data.table.R") # https://gist.github.com/etiennebr/9515738#file-as-data-table-r
# My functions 
source("code/calc_mahalanobis_fn.R")
source("code/find_analogs_fn.R")
source("code/calc_sigma_fn.R")
source("code/analog_impact_fn.R")
source("code/create_tile_fn.R")
# ------------------------------------------------------------------------------

# Path to stored climate data
gis_dir <- "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/GIS_Data"
focal_climate_dir <- file.path(gis_dir, "TopoTerra", "Annual futures")
analog_climate_dir <- file.path(gis_dir, "TopoTerra", "Normals")

# Read in analog data - function expects realization (ie, normal/mean vs annual)
# May choose to calculate normals from annuals
# No tiling - global search
norm_hist <- rast(file.path(analog_climate_dir,"topoterra_hist_1961-1990.tif")) |>
  as.data.table(xy = TRUE)

# Create a n*n tile grid across the study area
tile_grid <- terra::rast(paste0(focal_climate_dir, "/topoterra_2C_2000.tif")) |>
  sf::st_make_grid(n = 100)

#st_write(tile_grid, "tile_grid_100.gpkg")

# Here, iterate over grid of tiles ---------------------------------------------
profvis::profvis({

seq_len(length(tile_grid))[8831] |>
  walk(\(tile_i){ print(tile_i)

    # Define tile in tile_grid
    tile  <- tile_grid[[tile_i]] #198162

    # Test if empty
    test <- terra::extract(
      terra::rast(
        list.files(focal_climate_dir, full.names = TRUE)[[1]]
      ), 
      vect(tile)
    )
    if(all(is.na(test[,-1]))){
      return(NULL)
    }

    # Create tile of focal data
    focal_ann <- create_tile(
      tile, 
      climate_dir = focal_climate_dir,
      annual = TRUE
    )

    # Conduct analog search for all cells of tile
    # Calls: calc_mahalanobis and calc_sigma
    find_analogs(
      focal_data_cov = focal_ann,
      focal_data_mean = focal_ann,
      analog_data = norm_hist,
      var_names = c('aet', 'def', 'tmax', 'tmin'), 
      n_analog_pool = 1000000, 
      n_analog_use = 1000, 
      min_dist = 0, # In KM! 
      output_dir = paste0("data/test_",tile),
      use_futures = FALSE,
      n_futures = 1
    )
  })
  
})
