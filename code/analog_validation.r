# ------------------------------------------------------------------------------
# This script performs validation and performance testing of a 
# climate analog impact model for future vegetation

# Packages and functions -------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(furrr)
library(purrr)
library(data.table)
library(here)
library(exactextractr)
library(tmap)
library(ggplot2)

# Convert SpatRast to data.table
# https://gist.github.com/etiennebr/9515738#file-as-data-table-r
source(here("code/as.data.table.r"))

# Data setup -------------------------------------------------------------------
climate_dir <- file.path("data", "climate")

# Read in climate data, crop, convert to data.table 
# 2C annuals
ann_fut <- grep("topoterra_2C_[0-9]{4}.tif",
                list.files(climate_dir),
                value = TRUE) |>
  paste0(climate_dir, "/", ... = _) |>
  map(\(x) rast(x))

# Create validation set ------------------------------------------------------------------
# From the West-wide climate data, sample a subset of points to run
# validation on.
# Create a n*n tile grid across the study area
# tile_grid <- st_make_grid(ann_fut[[1]], n = 300)
# Sample n of these tiles
# tile_grid_samp <- vect(sample(tile_grid, 20))
huc12s <- st_read(here("data/huc12/huc12_western_states.gpkg")) 
random_hucs <- huc12s[sample(1:length(huc12s[['objectid']]), size = 10),] |>
  st_transform(crs = st_crs(ann_fut[[1]]))
#plot(random_hucs$geom)


# Mask each element of ann_fut to these sample tiles
fut_samp <- ann_fut |>
    map(\(x) mask(x, random_hucs)) 

# Convert each element into a data.table
fut_samp <- lapply(fut_samp, FUN = as.data.table, xy = TRUE)

# For contemporary validation, the focal cells will be historical normals, but we will 
# build the covariance matrix based on the future annuals, because that works here
# Historical normals (downscaled from TerraClimate normals)
norm_hist <- grep("topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
                  list.files(climate_dir),
                  value = TRUE) |>
  paste0(climate_dir, "/", ... = _) |>
  rast() |>
  as.data.table(xy = TRUE)

# Create the same random sample of tiles for the historical data
hist_samp <- grep("topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
                  list.files(climate_dir),
                  value = TRUE) |>
  paste0(climate_dir, "/", ... = _) |>
  rast() |>
  mask(random_hucs) |>
  as.data.table(xy = TRUE)


# Remove things
rm(ann_fut, huc12s)
gc()

# Select contemporary analogs for validation set --------------------------------------------------
source("code/mahalanobis_D_fn.R")

# This takes several hours for 17867 points...
#library(profvis)
start = Sys.time()
#profvis({
find_analogs(focal_data_cov = fut_samp,
             focal_data_mean = hist_samp,
             analog_data = norm_hist,
             var_names = c('aet', 'def', 'tmax', 'tmin'), 
             n_analog_pool = 5000000, 
             n_analog_keep = 5000, 
             min_dist = 50000, 
             output_dir = "data/validation/5M-pool_5k-keep_50km-min",
             use_futures = TRUE,
             n_futures = 6)
#})
end = Sys.time()
start-end
# ------------------------------------------------------------------------------------------------


