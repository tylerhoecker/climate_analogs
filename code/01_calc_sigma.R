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
# TEMPORARILY REMOVE TMIN/TMAX DUE TO ERROR
ann_fut <- lapply(ann_fut, function(i) select(i, -tmin, -tmax))

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
#   # TEMPORARILY REMOVE TMIN/TMAX DUE TO ERROR
#   select(-tmin, -tmax)

# Historical normals (downscaled from TerraClimate normals)
norm_hist <- grep('topoterra_hist_[0-9]{4}-[0-9]{4}.tif', 
                  list.files(climate_dir), 
                  value = T) %>%
  paste0(climate_dir,'/',.) %>% 
  rast(.) %>% 
  as.data.table(.) %>% 
  # TEMPORARILY REMOVE TMIN/TMAX DUE TO ERROR
  select(-tmin, -tmax)

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
    summarise(across(c(aet, def), mean)) %>% 
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
  sigma <- sqrt(sigma)

  # Bind values to coordinates of best analogs
  ref_coords_sample <- ref_coords[random_pts,]

  out_dt <- data.table('x' = ref_coords_sample[,'x'],
                       'y' = ref_coords_sample[,'y'],
                       'D' = D,
                       'sigma' = sigma) %>% 
    filter(sigma > 0, is.finite(sigma))
  
  return(out_dt)
}

plan(callr, workers = 5)

sigma_df <- seq(nrow(ann_fut[[1]])) %>% 
  future_map_dfr(md_fun)

saveRDS(sigma_df, '../data/sigma_df_11-9-23.RD')  










