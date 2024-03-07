# ------------------------------------------------------------------------------
# This script performs validation and performance testing of a
# climate analog impact model for future vegetation

# Packages and functions -------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(data.table)

# Convert SpatRast to data.table
# https://gist.github.com/etiennebr/9515738#file-as-data-table-r
source(here("code/as.data.table.r"))
  
# Data setup -------------------------------------------------------------------
gis_dir <- "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/GIS_Data/"

# Move to script/function
# -------------------------------------------------------------------------------
# Read in climate data, crop, convert to data.table
# 2C annuals
# ann_fut <- grep(
#   "topoterra_2C_[0-9]{4}.tif",
#   list.files(
#     paste0(gis_dir, "TopoTerra/Annual futures"),
#     full.names = TRUE
#   ),
#   value = TRUE
#   # We have 31 years of data... drop 1985 to make it an even 30-yr normal
#   )[2:31] |>
#   map(\(x) rast(x))

# Create validation set ------------------------------------------------------------------
# - TILES
# Create a n*n tile grid across the study area
# tile_grid <- st_make_grid(ann_fut[[1]], n = 300)
# Sample n of these tiles
# tile_grid_samp <- vect(sample(tile_grid, 20))
# - HUCS
# huc12s <- st_read(here(gis_dir, "huc12/huc12_western_states_good.gpkg"))
# random_hucs <- huc12s[sample(1:length(huc12s[["objectid"]]), size = 10), ] |>
#   st_transform(crs = st_crs(ann_fut[[1]]))
# plot(random_hucs$geom)
# write_sf(random_hucs, "data/random_hucs.gpkg")
random_hucs <- read_sf("data/random_hucs.gpkg")

# Crop to random_hucs
# Mask each element of ann_fut to these sample tiles
# fut_samp <- ann_fut |>
#   map(\(x) mask(x, random_hucs))
# # Convert each element into a data.table
# fut_samp <- lapply(fut_samp, FUN = as.data.table, xy = TRUE)
# # This is slow, save for convenience for now
# saveRDS(fut_samp, "fut_samp.Rds")
fut_samp <- readRDS("fut_samp.Rds")

# For contemporary validation, the focal cells will be historical normals, but we will
# build the covariance matrix based on the future annuals, because that works here
# Historical normals (downscaled from TerraClimate normals)
norm_hist <- grep(
  "topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
  list.files(
    paste0(gis_dir, "TopoTerra/Normals"),
    full.names = TRUE
  ),
  value = TRUE
) |>
  rast() |>
  as.data.table(xy = TRUE)

# Create the same random sample of tiles for the historical data
hist_samp <- grep(
  "topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
  list.files(
    paste0(gis_dir, "TopoTerra/Normals"),
    full.names = TRUE
  ),
  value = TRUE
) |>
  rast() |>
  mask(random_hucs) |>
  as.data.table(xy = TRUE)
