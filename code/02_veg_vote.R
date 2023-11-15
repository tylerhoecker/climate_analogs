library(tidyverse)
library(terra)
library(sf)

# Point to sigma output
simga_path <- '~/Google Drive/My Drive/analogs_transfer/sigma_output'
sigma_output <- list.files(simga_path, full.names = T) 

# OR-WA boundary for crop ------------------------------------------------------
climate_dir <- file.path('..','data','climate')
new_crs <- crs(rast(paste0(climate_dir,'/topoterra_2C_2000.tif')))

# Load landfire data
lf_bps <- rast('../data/landfire/lf_bps_west.tif') 
# w_states <- read_sf('../data/western_states/western_states.shp') %>% 
#   st_transform(., new_crs)
# lf_bps <- crop(lf_bps, w_states)
# writeRaster(lf_bps, '../data/landfire/lf_bps_west.tif')

# Veg extraction and vote fn ---------------------------------------------------

# Veg voting function
chunk <- readRDS(sigma_output[[1]])

# Extract BPS group for each focal cell in this chunk
focal_pts <- unique(chunk[['focal_id']])

focal_pt <- focal_pts[[1]]

veg_vote_fn <- function(focal_pt){
  
  profvis({
    lf_bps <- rast('../data/landfire/lf_bps_west.tif')
    
    focal_df <- chunk %>% 
      # Select focal point from chunk
      filter(focal_id == focal_pt) %>% 
      # Extract Landfire 2020 BPS code
      mutate(bps = extract(lf_bps, 
                           as.matrix(.[,c('an_x','an_y')]))[['lf_bps_32611']],
             # Rescale sigma to 0-1, to use as weight for vote
             sigma_z = scales::rescale(sigma),
             # Record the total number of analogs (< 2 sigma) in pool
             n_analogs = n()) %>% 
      # For each unique BPS code...
      group_by(bps) %>% 
      # Sum the weighted vote (sigma 0-1)
      summarize(weighted_vote = sum(sigma_z),
                # Count the number of analogs in pool w/ this BPS code
                raw_vote = n(),
                # Record the total number of analogs in pool
                n_analogs = first(n_analogs),
                # Calculate the proportion of votes for BPS code
                weighted_prop = weighted_vote/n_analogs,
                raw_prop = raw_vote/n_analogs) %>% 
      # Pick plurality vote
      arrange(desc(weighted_prop)) %>% 
      slice_head(n = 1)
    
  })
  

  return(focal_df)
}

profvis({map_dfr(focal_pts[[1]], veg_vote_fn, .id = 'focal_pt', .progress = TRUE)})

library(furrr)
availableCores()
plan(multisession, workers = 10)
options(future.globals.maxSize= 1500000000)
library(tictoc)
tic()
veg_df <- map_dfr(focal_pts[[1]], veg_vote_fn, .id = 'focal_pt', .progress = TRUE)
toc()
plan(sequential)

#-------------------------------------------------------------------------------
# Simplify to top vote

veg_df_1 <- veg_df %>% 
  group_by(focal_pt) %>% 
  slice_head(n = 1) %>% 
  select()

chunk_coords <- chunk %>% 
  group_by(focal_id) %>% 
  slice_head(n = 1) %>% 
  select(focal_id, x = fc_x, y = fc_y)

out_rast <-  chunk_coords[,c('x','y')] |> 
  cbind(veg_df_1[,'bps']) |> 
  rast(type = 'xyz', crs = new_crs)










  