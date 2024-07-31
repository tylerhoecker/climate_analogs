# Some work was done to prepare the Lanfire BPS data for this analysis

# Resample landfire data to the same grid as the climate data, so BPS classes are the mode of all 
# 30-m pixels in the 220-m grid, not a random centroid
# Load landfire data - cropped and tranformed in QGIS
lf_bps <- rast('../data/landfire/lf_bps_west.tif') 
clim_grid <- rast('../data/climate/topoterra_2C_1985-2015.tif')
lf_bps <- resample(lf_bps, clim_grid, method = 'mode')
writeRaster(lf_bps, '../data/landfire/lf_bps_west_220.tif')
# This is slooow.... read in the save version


# Just saving this chunk in case 
if(simlplify_bps == TRUE){
    # BPS simplification scheme to reclassify these 
    bps_key <- read_csv('../data/bps_simplification.csv')
    bps_rast <- classify(lf_bps, cbind(bps_key$VALUE,bps_key$bps_group))
  } else {
    bps_rast <- lf_bps
  }
  