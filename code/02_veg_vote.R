library(tidyverse)
library(terra)
library(sf)
library(data.table)
# Point to sigma output
simga_path <- '../data/sigma_output'
sigma_output <- list.files(simga_path, full.names = T) 

chunks_coords <- readRDS('../data/coord_chunks.Rds')

# Resample landfire data to the same grid as the climate data, so BPS classes are the mode of all 
# 30-m pixels in the 220-m grid, not a random centroid
# Load landfire data - cropped and tranformed in QGIS
# lf_bps <- rast('../data/landfire/lf_bps_west.tif') 
# clim_grid <- rast('../data/climate/topoterra_2C_1985-2015.tif')
# lf_bps <- resample(lf_bps, clim_grid, method = 'mode')
# writeRaster(lf_bps, '../data/landfire/lf_bps_west_220.tif')
# This is slooow.... read in the save version

# Veg extraction and vote fn ---------------------------------------------------

# For testing
lf_bps <- rast('../data/landfire/lf_bps_west_220.tif')

veg_vote_fn <- function(chunk, chunk_idx){
  
  veg_chunk_df <- chunk %>% 
    readRDS() %>% 
    mutate(focal_id = as.numeric(focal_id)) %>% 
    # Select focal point from chunk
    # filter(focal_id == focal_pt) %>%
    # Extract Landfire 2020 BPS code
    mutate(bps = extract(lf_bps, as.matrix(.[,c('an_x','an_y')]))[['lf_bps_32611']]) %>% 
    group_by(focal_id) %>% 
    # Rescale sigma to 0-1, to use as weight for vote
    mutate(sigma_z = scales::rescale(sigma),
           # Record the total number of analogs (< 2 sigma) in pool
           n_analogs = n()) %>% 
    # For each unique BPS code...
    group_by(focal_id, bps) %>% 
    # Sum the weighted vote (sigma 0-1)
    summarize(wt_vote = sum(sigma_z),
              # Count the number of analogs in pool w/ this BPS code
              raw_vote = n(),
              # Record the total number of analogs in pool
              n_analogs = first(n_analogs),
              # Calculate the proportion of votes for BPS code
              wt_prop = wt_vote/n_analogs,
              raw_prop = raw_vote/n_analogs) %>% 
    # This sequence saves the rows for each focal point, one for each BPS class,
    # that account for 80% of the (unweighted) analogs 
    arrange(desc(raw_prop), .by_group = T) %>% 
    filter(
      cumsum(
        raw_prop == accumulate(
          raw_prop, ~ ifelse(.x <= .80, .x + .y, .y))) == 1) %>% 
    # Could save the winning BPS class by raw vote instead/addition to weighted BPS
    # mutate(raw_bps = first(bps)) %>% 
    # Instead, sort by the weighted votes for BPS classes
    arrange(desc(wt_prop), .by_group = T) %>% 
    mutate(#w_bps = first(bps),
      n_bps = length(bps)) %>% 
    # And save the 
    slice_head(n = 3) %>%  
    summarise(bps1 = bps[1],
              bps2 = bps[2],
              bps3 = bps[3],
              win_prop = wt_prop[1],
              n_analogs = n_analogs[1],
              n_bps = n_bps[1])  
  
  current_bps <- extract(lf_bps, chunks_coords[[chunk_idx]])[['lf_bps_32611']]
  
  veg_chunk_df <- cbind(chunks_coords[[chunk_idx]],
                        veg_chunk_df,
                        'current_bps' = current_bps)
}


veg_votes_df <- sigma_output %>% 
  imap_dfr(veg_vote_fn)
  

#-------------------------------------------------------------------------------
# Map to raster
out_rast <-  veg_votes_df %>% 
  select(x, y, bps1) %>% 
  rast(., type = 'xyz', crs = crs(lf_bps))
  









  