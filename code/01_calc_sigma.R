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
source("code/as.data.table.r")

# Boundary for crop ------------------------------------------------------
climate_dir <- file.path("data", "climate")

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

# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]] # nolint

md_fun <- function(pt_i) {
  
  # Build cov matrix for pt_i from 30 years of annual projected future data
  cov_i <- map(chunk, ~ .x[pt_i]) |>
    rbindlist() |>
    dplyr::select(-x, -y) |>
    cov()

  # Calculate the mean of the future annuals -----------------------------------
  # This provides different result than calculating the mean here...
  # probably a result of the downscaling process
  fut_mean <- map(chunk, ~ .x[pt_i]) |>
    rbindlist() |>
    summarise(across(c(aet, def, tmax, tmin), mean)) |>
    unlist()
 
  # Analog pool - from historical normals --------------------------------------
  # This *could* be done outside the function, but this is less biased, because
  # random samples are drawn with each iteration,
  # ultimately sampling from (probably)
  # nearly the entire landscape
  
  # This could be integrated with the chunk below, using slice_sample,
  # but saving this index will be useful for extracting coordinates later
  random_pts <- seq_len(nrow(norm_hist)) |>
    sample(size = n_analogs)
  
  # Build reference matrix from random sample of historical normals
  ref_mat <- norm_hist[random_pts, ] |>
    dplyr::select(-x, -y) |>
    as.matrix()
  
  # Sigma dissimilarity between pt_i and analog pool ---------------------------
  # Calculate 'raw' unsquared Mahalanobis distances
  d <- mahalanobis(ref_mat, fut_mean, cov_i)
  
  # Mahoney follows a procedure where D is immediately unsquared,
  # then is fed into a chi distribution function.
  # Mathematically, it is equivalent to feed the squared-D into
  # a chi-squared distribution function, and then unsquare the result. 
  # stats::chisq is an order of magnitude faster than chi::qchi.
  
  # Convert distances to percentiles of chi distribution (mulit-dimensional normal)
  # df = number of dimensions / climate variables
  p <- pchisq(d, df = length(fut_mean))
  # Convert percentiles into quantiles (standard deviations from mean)
  sigma <- qchisq(p, df = 1) # df is now 1...
  sigma <- sqrt(sigma) |>
    round(3)

  # Save output
  out_dt <- data.table("f_x" = chunk[[1]][pt_i, "x"],
                       "f_y" = chunk[[1]][pt_i, "y"],
                       "a_x" = norm_hist[random_pts, "x"],
                       "a_y" = norm_hist[random_pts, "y"],
                       "sigma" = sigma)
  # Not sure why .x/.y get's appended above... data.table thing
  setnames(out_dt,
           old = c("f_x.x", "f_y.y", "a_x.x", "a_y.y"),
           new = c("f_x", "f_y", "a_x", "a_y"))

  out_dt <- out_dt |>
    dplyr::filter(is.finite(sigma), sigma < 2) |>
    mutate("dist_m" = round(sqrt((f_x - a_x)^2 + (f_y - a_y)^2))) |>
    arrange(sigma) |>
    slice_head(n = 1000)

  rm(cov_i, fut_mean, random_pts, ref_mat, d, p, sigma)
  gc()
  return(out_dt)
}

# Parallelize and run!
options(future.globals.maxSize = 5000000000)
plan(multicore, workers = 46)

chunks |>
  iwalk(\(chunk, chunk_idx) {
     
    # Calculate sigmas to analogs for each focal point
    sigma_dt <- seq_len(nrow(chunk[[1]])) |>
      future_map_dfr(md_fun, .id = "focal_id", .progress = TRUE)
    
    # Save as RDS
    saveRDS(sigma_dt,
            paste0("data/sigma_output/sigma_dt_", chunk_idx, ".Rds"))
    
    # Explicitly remove and free memory each iteration
    rm(sigma_dt)
    gc()
  })