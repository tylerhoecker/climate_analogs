#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

# Packages and functions -------------------------------------------------------
library(sf)
library(tidyverse)
library(terra)
source('as.data.table.r') # https://gist.github.com/etiennebr/9515738#file-as-data-table-r
library(chi)
library(furrr)
library(callr)

# OR-WA boundary for crop ------------------------------------------------------
climate_dir <- file.path('..','data','climate')

new_crs <- crs(rast(paste0(climate_dir,'/topoterra_2C_2000.tif')))

w_states <- read_sf('../data/western_states/western_states.shp') %>% 
  st_transform(., new_crs)

r6 <- w_states %>% 
  filter(NAME %in% c('Washington'))

# Read in climate data, crop, convert to data.table ----------------------------
# 2C annuals
ann_fut <- grep('topoterra_2C_[0-9]{4}.tif', 
                list.files(climate_dir), 
                value = T) %>%
  paste0(climate_dir,'/',.) %>% 
  map(\(x) crop(rast(x),r6)) 

ann_fut <- lapply(ann_fut, FUN = as.data.table)

fut_coords <- paste0(climate_dir,'/','topoterra_2C_1985.tif') %>% 
  rast() %>% 
  crop(., r6) %>% 
  crds()
# This is more consistent, but doesn't work
# ann_fut <- ann_fut %>% 
#   map(., as.data.table())

# 2C normals
# Eventually, remove this because the downscaled normals are slightly different than 
# the means/normals calculated from the downscaled annuals 
# norm_fut <- grep('topoterra_2C_[0-9]{4}-[0-9]{4}.tif',
#                  list.files(climate_dir),
#                  value = T) %>%
#   paste0(climate_dir,'/',.) %>%
#   rast(.) %>%
#   crop(., r6) %>%
#   as.data.table(.) %>%

# Historical normals (downscaled from TerraClimate normals)
norm_hist <- grep('topoterra_hist_[0-9]{4}-[0-9]{4}.tif', 
                  list.files(climate_dir), 
                  value = T) %>%
  paste0(climate_dir,'/',.) %>% 
  rast(.) %>% 
  as.data.table(.) 
 
# Save out the centroid coordinates separately
ref_coords <- grep('topoterra_hist_[0-9]{4}-[0-9]{4}.tif', 
                   list.files(climate_dir), 
                   value = T) %>%
  paste0(climate_dir,'/',.) %>% 
  rast(.) %>% 
  crds()

rm(new_crs, w_states, r6)
gc()
#-------------------------------------------------------------------------------

# Calculate sigma dissimilarity between all focal locations and n best analogs
#-------------------------------------------------------------------------------
# How many analogs randomly sampled from global pool?
n_analogs <- 1000000
pt_i = 10000
# The MD/sigma dissimilarity function
md_fun <- function(pt_i){
  print(paste0('Running point ', pt_i, ' of ', nrow(ann_fut[[1]])))
  # Build covariance matrix for pt_i from 30 years of annual projected future data
  cov_i <- map(ann_fut, ~ .x[pt_i]) %>%
    rbindlist() %>%
    cov()
 
  # Calculate the mean of the future annuals -----------------------------------
  # This provides different result than calculating the mean here... 
  # probably a result of the downscaling process
  # fut_mu_vec <- norm_fut[pt_i] %>%
  #   unlist()
  fut_mean <- map(ann_fut, ~ .x[pt_i]) %>%
    rbindlist() %>% 
    summarise(across(c(aet, def, tmax, tmin), mean)) %>% 
    unlist()
  
  # Analog pool - from historical normals --------------------------------------
  # This *could* be done outside the function, but this is less biased, because
  # random samples are drawn with each iteration, ultimately sampling from (probably)
  # nearly the entire landscape
  
  # This could be integrated with the chunk below, using slice_sample, 
  # but saving this index will be useful for extracting coordinates later
  random_pts <- seq(nrow(norm_hist)) %>% 
    sample(size = n_analogs)
  
  ref_mat <- norm_hist[random_pts] %>% 
    as.matrix() 
  
  # Sigma dissimilarity between pt_i and analog pool ---------------------------
  # Calculate 'raw' unsquared Mahalanobis distances
  D <- mahalanobis(ref_mat, fut_mean, cov_i) 
  
  # Mahoney follows a procedure where D is immediately unsquared, then is fed into 
  # a chi distribution function. Mathematically, it is equivalent to feed the squared-D into 
  # a chi-squared distribution function, and then unsquare the values. This is an order of 
  # magnitude faster, because the chi::qchi function is very slow.
  # %>% 
  #   sqrt()
  
  # Convert distances (D) to percentiles of chi distribution (mulit-dimensional normal)
  P <- pchisq(D, df = 2) 
  # Convert percentiles into quantiles (standard deviations from mean)
  sigma <- qchisq(P, df = 1) 
  sigma <- sqrt(sigma) %>% 
    round(., 3)

  # Bind values to coordinates of best analogs
  ref_coords_sample <- ref_coords[random_pts,]

  out_dt <- data.table('fc_x' = fut_coords[pt_i,'x'],
                       'fc_y' = fut_coords[pt_i,'y'],
                       'an_x' = ref_coords_sample[,'x'],
                       'an_y' = ref_coords_sample[,'y'],
                       'md' = D,
                       'sigma' = sigma) %>% 
    filter(is.finite(sigma), sigma < 2) %>% 
    mutate('dist_m' = sqrt((fc_x-an_x)^2 + (fc_y-an_y)^2)) 
  
  rm(cov_i, fut_mean, random_pts, ref_mat, D, P, sigma, ref_coords_sample)
  gc()
  return(out_dt)
}

# Temp plotting during development
# some_pts <- out_dt %>% 
#   slice_sample(n = 10000)
# 
# ggplot(out_dt, aes(x = dist_km, y = sigma)) +
#   geom_point(data = some_pts) +
#   geom_smooth()

# Break the focal dataset into chunks 
# ------------------------------------------------------------------------------
# Set the number of chunks 
num_chunks <- 350
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
  # Slice each dataframe in the list to get the subset of rows for the current chunk
  chunk <- lapply(ann_fut, function(df) df[start_index:end_index, , drop = FALSE])
  # Append the chunk to the list of chunks
  chunks[[i]] <- chunk
}
# ------------------------------------------------------------------------------


# Parallelize and run!
plan(multicore, workers = 40)
chunk_idx <- seq(length(chunks))

chunks %>% 
  iwalk(\(chunk, chunk_idx){
    
    sigma_dt <- seq(nrow(chunk[[1]])) %>% 
      future_map_dfr(md_fun, .id = 'focal_id')
    
    saveRDS(sigma_dt, paste0('../data/sigma_output/sigma_dt_',chunk_idx,'.rds'))  
    
  })











