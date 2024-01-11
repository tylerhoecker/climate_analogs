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
  map(\(x) crop(rast(x), r6))

ann_fut <- lapply(ann_fut, FUN = as.data.table, xy = TRUE)

# Historical normals (downscaled from TerraClimate normals)
norm_hist <- grep("topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
                  list.files(climate_dir),
                  value = TRUE) |>
  paste0(climate_dir, "/", ... = _) |>
  rast() |>
  as.data.table(xy = TRUE)
