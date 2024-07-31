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
focal_climate_dir <- file.path("data","climate","Annual futures") #file.path(gis_dir, "TopoTerra", "Annual futures")
analog_climate_dir <- file.path("data","climate","Normals") #file.path(gis_dir, "TopoTerra", "Normals")
# Focal ecoregions
ecoregs <- sf::st_read("/home/hoecker/GitHub/climate_analogs/data/focal ecoregions/focal_ecoregions_3857.gpkg") 

# Read in analog data - function expects realization (ie, normal/mean vs annual)
# May choose to calculate normals from annuals
# No tiling - global search
norm_hist <- rast(file.path(analog_climate_dir,"topoterra_hist_1961-1990.tif")) |>
  #mask(ecoregs) |>
  as.data.table(xy = TRUE)

# Create a n*n tile grid across the study area
# tile_grid <- terra::rast(paste0(focal_climate_dir, "/topoterra_2C_2000.tif")) |>
#   sf::st_make_grid(n = 50)

# Create tiles and save them (this failed in parallel so doing it sequentially first, takes ~1hr)
# seq_len(length(tile_grid)) |>
#   walk(\(x){ 
#     # Create tile of focal data
#     create_tile(
#       tile = x, 
#       tile_grid = tile_grid,
#       climate_dir = focal_climate_dir,
#       annual = TRUE,
#       write_out = "data/ecoreg_dt_tiles/"
#     )
#   })

# Find_analogs for every point in every tile
# Set up parallel "futures"
gc()
options(future.globals.maxSize = 5000000000, future.seed = TRUE)
plan(multicore, workers = 30)

unfinished <- list.files("data/ecoreg_dt_tiles", full.names = TRUE) |>
  map(\(x){
    # Check if complete
    if (!file.exists(paste0("data/ecoreg_output/analogs_",basename(x)))){
      return(x)
    } 
  }
  ) |>
  compact()

unfinished |>
  future_walk(\(x){

    # Progress message
    print(paste0("Running ", x))
    
    
    # Read in the datatable version of the tile
    focal_ann <- readRDS(x)

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
      output_dir = paste0(
        "data/ecoreg_output/analogs_",
        sub(".Rds","",basename(x))
      ),
      use_futures = FALSE,
      n_futures = 40
    )
  })

