library(terra)
library(sf)
library(tidyverse)
library(purrr)
source("code/analog_impact_fn.R")

# Paths to input files
gis_dir <- "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/GIS_Data"
tiles_dir <- "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/analogs_transfer/westwide_tiles_reduced_april"

# Load "impact" data, here, Landfire BPS
lf_bps_220 <- terra::rast(file.path(gis_dir, "Landfire/lf_bps_west_220.tif")) |>
  # Reprojecting to the EPSG that Josh Gage requires, which seems good as canonical
  terra::project("EPSG:3857", method = "mode")

# Map over all tiles
tile_dts <- list.files(tiles_dir, full.names = TRUE)[1:10]
tile_dts <- c(
  "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/analogs_transfer/westwide_tiles_reduced_april/analogs_tile_1965.Rds",
  "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/analogs_transfer/westwide_tiles_reduced_april/analogs_tile_288.Rds",
  "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/analogs_transfer/westwide_tiles_reduced_april/analogs_tile_1877.Rds"
)

min_dist = 0
max_dist = Inf

aim_result <- tile_dts |>
    map_dfr(\(tile_i){

        # Read in the dataframe for this tile
        tile_dt <- readRDS(tile_i)

        # Skip empty ones
        if(nrow(tile_dt) == 0){
          return(NULL)
        }

        # Evaluate validation parameters
        tile_dt <- tile_dt |>
          filter(dist_km > min_dist, dist_km < max_dist) |>
          # Could eventually remove if not reducing # of analogs
          group_by(focal_id) |>
          arrange(md) |>
          slice_head(n = n_analog_use) 
        
        impact_results  <- analog_impact_fn(
          analog_data = tile_dt,
          impact_data = lf_bps_220,
          n_analog_use = 1000,
          weight_val = "dist_km",
          n_projections = 5
        )

        return(impact_results)
      },
      .progress = TRUE 
    )  |>
    ungroup()

list("aim_1", "aim_2", "aim_3", "aim_4", "aim_5") %>%  #'current_bps','bps1','bps2','bps3','win_prop_raw','n_analogs','n_bps','mean_sigma_100') %>%  
  walk(\(band) {
    
    print(paste0('Rasterizing ',band,'...'))
    
    r <- aim_result |> 
      ungroup() |>
      select(x, y, all_of(band))  |> 
      terra::rast(x = _, 
                  type = 'xyz', 
                  crs = crs(lf_bps_220)
                 ) 
        
    writeRaster(r, 
                paste0("data/westwide_reduced_aim/",band,'.tif'), 
                datatype = "INT2S",
                gdal = c("PROJECTION=EPSG:3857",
                         "TILED=YES",
                         "BLOCKXSIZE=128",
                         "BLOCKYSIZE=128",
                         "OVERVIEW-RESAMPLING=NEAREST",
                         "COMPRESS=DEFLATE"),
                overwrite = TRUE)

    return(r)
  }) 
