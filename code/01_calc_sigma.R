#-------------------------------------------------------------------------------
# Testing
# historical <- rast('../data/climate/topoterra_hist_1961-1990.tif')
# future <- rast('../data/climate/topoterra_2C_1985-2015.tif')
#-------------------------------------------------------------------------------

# Packages and functions -------------------------------------------------------
library(terra)
library(sf)
library(tidyr)
library(dplyr)
library(furrr)
library(purrr)
library(data.table)
# https://gist.github.com/etiennebr/9515738#file-as-data-table-r
source("code/as.data.table.R")
source("code/find_analogs_fn.R")

climate_dir <- file.path("data", "climate")

# Boundary for crop ------------------------------------------------------
new_crs <- crs(rast(paste0(climate_dir, "/topoterra_2C_2000.tif")))

w_states <- read_sf("data/western_states/western_states.shp") |>
  st_transform(new_crs)

r6 <- w_states  |>
  filter(NAME %in% c("Colorado"))

# Read in climate data, crop, convert to data.table ----------------------------
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

rm(new_crs, w_states, r6)

gc()
#-------------------------------------------------------------------------------

# Break the focal dataset into chunks
# ------------------------------------------------------------------------------
# Set the number of chunks
num_chunks <- 300
# Get the total number of rows in each dataframe
total_rows <- nrow(ann_fut[[1]])
# Calculate the number of rows in each chunk
rows_per_chunk <- total_rows %/% num_chunks
# Initialize an empty list to store chunks
chunks <- list()

# Iterate over the range of chunks
for (i in 1:num_chunks) {
  # Calculate the start and end indices for each chunk
  start_index <- (i - 1) * rows_per_chunk + 1
  end_index <- min(i * rows_per_chunk, total_rows)
  # Get the subset of rows for the current chunk
  chunk <- lapply(ann_fut,
                  function(df) df[start_index:end_index, , drop = FALSE])
  # Append the chunk to the list of chunks
  chunks[[i]] <- chunk
}
rm(ann_fut, chunk)
gc()
# ------------------------------------------------------------------------------

# Calculate sigma dissimilarity between all focal locations and analogs
#-------------------------------------------------------------------------------
# How many analogs randomly sampled from global pool?
n_analogs <- 10000000

# Parallelize and run!
options(future.globals.maxSize = 5000000000)
plan(multicore, workers = 46)

chunks |>
  iwalk(\(chunk, chunk_idx) {
     
    # Calculate sigmas to analogs for each focal point
    sigma_dt <- seq_len(nrow(chunk[[1]])) |>
      future_map_dfr(find_analogs, .id = "focal_id", .progress = TRUE)
    
    # Save as RDS
    saveRDS(sigma_dt,
            paste0("data/sigma_output/sigma_dt_", chunk_idx, ".Rds"))
    
    # Explicitly remove and free memory each iteration
    rm(sigma_dt)
    gc()
  })