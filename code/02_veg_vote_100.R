library(tidyverse)
library(terra)
library(sf)
library(data.table)
# Point to sigma output
simga_path <- '../data/sigma_output_11-24'
sigma_output <- list.files(simga_path, full.names = T) 

# Resample landfire data to the same grid as the climate data, so BPS classes are the mode of all 
# 30-m pixels in the 220-m grid, not a random centroid
# Load landfire data - cropped and tranformed in QGIS
# lf_bps <- rast('../data/landfire/lf_bps_west.tif') 
# clim_grid <- rast('../data/climate/topoterra_2C_1985-2015.tif')
# lf_bps <- resample(lf_bps, clim_grid, method = 'mode')
# writeRaster(lf_bps, '../data/landfire/lf_bps_west_220.tif')
# This is slooow.... read in the save version

# Veg extraction and vote fn ---------------------------------------------------
#chunk = sigma_output[[100]]
#simlplify_bps = FALSE
# Landfire BPS data
lf_bps <- rast('../data/landfire/lf_bps_west_220.tif')
names(lf_bps) <- 'bps'

veg_vote_fn <- function(chunk, simlplify_bps){
  
  if(simlplify_bps == TRUE){
    # BPS simplification scheme to reclassify these 
    bps_key <- read_csv('../data/bps_simplification.csv')
    bps_rast <- classify(lf_bps, cbind(bps_key$VALUE,bps_key$bps_group))
  } else {
    bps_rast <- lf_bps
  }
  
  veg_chunk_df <- chunk %>% 
    readRDS() %>% 
    mutate(focal_id = as.numeric(focal_id)) %>% 
    # 1000 best analogs for each point... 
    group_by(focal_id) %>% 
    arrange(sigma, .by_group = T) %>%
    slice_head(n = 100) %>% 
    ungroup() %>% 
    # Extract Landfire 2020 BPS code
    mutate(bps = extract(bps_rast, as.matrix(.[,c('an_x','an_y')]))[['bps']]) %>% 
    # For each unique BPS code...
    group_by(focal_id, bps) %>% 
    # Sum the weighted vote (sigma 0-1)
    summarize(# Preserve coordinates
              x = first(x), y = first(y),
              # Count the number of analogs in pool w/ this BPS code
              votes = n(),
              prop = votes/100,
              mean_sigma = mean(sigma),
              dist_bps_km = round(min(dist_m)*0.001, 0)) %>% 
    # This sequence saves the rows for each focal point, one for each BPS class,
    # that account for 80% of the (unweighted) analogs 
    arrange(desc(prop), .by_group = T) %>% 
    filter(
      cumsum(
        prop == accumulate(
          prop, ~ ifelse(.x <= .80, .x + .y, .y))) == 1) %>% 
    mutate(n_bps = length(bps),
           win_prop = first(prop)) %>%
    # Could save the winning BPS class by raw vote instead/addition to weighted BPS
    # mutate(raw_bps = first(bps)) %>% 
    # Instead, sort by the weighted votes for BPS classes
    # And save the 
    slice_head(n = 3) %>%  
    summarise(x = first(x), y = first(y),
              bps1 = bps[1],
              bps2 = bps[2],
              bps3 = bps[3],
              win_prop = first(prop),
              n_bps = first(n_bps),
              mean_sigma = first(mean_sigma),
              dist_bps_km = first(dist_bps_km))  
  
  current_bps <- extract(bps_rast, as.matrix(veg_chunk_df[,c('x','y')]))[['bps']]
  
  veg_chunk_df <- cbind(veg_chunk_df,current_bps)
  
  return(veg_chunk_df)
}

veg_votes_df <- map_dfr(sigma_output, veg_vote_fn, simlplify_bps = FALSE)
  
saveRDS(veg_votes_df, 'wa_veg_proj_fullBPS_fixedN_1000.Rds')
veg_votes_df <- readRDS('wa_veg_proj_fullBPS_fixedN_1000.Rds')

#-------------------------------------------------------------------------------
# Rasterize and write out the TIFFs
bps_key <- read_csv('../data/bps_simplification.csv')

list('current_bps','bps1','bps2','bps3','mean_sigma','dist_bps_km','win_prop','n_bps') %>%  #'current_bps','bps1','bps2','bps3','win_prop_raw','n_analogs','n_bps','mean_sigma_100') %>%  
  walk(function(band) {
    
    print(paste0('Rasterizing ',band,'...'))
    
    r <- veg_votes_df %>% 
      select(x, y, band) %>% 
      rast(., type = 'xyz', crs = crs(lf_bps)) 
    
    # if(band %in% c('current_bps','bps1','bps2','bps3')){
    #   r <- classify(r, cbind(bps_key$VALUE,bps_key$veg_num))
    # }
    
    r <- project(r, 'EPSG:4326', method = 'near')
    
    writeRaster(r, filename = paste0('../output/gee_fixedN/fullBPS/wa_',band,'.tiff'), overwrite = T)
    
    return(r)
  }) 




  