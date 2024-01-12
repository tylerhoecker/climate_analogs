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
# Convert SpatRast to data.table
# https://gist.github.com/etiennebr/9515738#file-as-data-table-r
source("code/as.data.table.r")

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
tile_grid <- st_make_grid(ann_fut[[1]], n = 300)
# Sample n of these tiles
tile_grid_samp <- vect(sample(tile_grid, 20))

# Mask each element of ann_fut to these sample tiles
fut_samp <- ann_fut |>
    map(\(x) mask(x, tile_grid_samp)) 

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
  mask(tile_grid_samp) |>
  as.data.table(xy = TRUE)


# Remove things
rm(ann_fut, tile_grid)
gc()

# Select analogs for validation set ----------------------------------------------------------
source("code/mahalanobis_D_fn.R")

find_analogs(focal_data_cov = fut_samp,
             focal_data_mean = hist_samp,
             analog_data = norm_hist,
             var_names = c('aet', 'def', 'tmax', 'tmin'), 
             n_analog_pool = 10000000, 
             n_analog_keep = 1000, 
             min_dist = 50000, 
             output_dir = "data/validation",
             use_futures = TRUE,
             n_futures = 4)
