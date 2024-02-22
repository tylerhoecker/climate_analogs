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
library(here)
library(exactextractr)
library(tmap)
library(ggplot2)

# Convert SpatRast to data.table
# https://gist.github.com/etiennebr/9515738#file-as-data-table-r
source(here("code/as.data.table.r"))
  
# Data setup -------------------------------------------------------------------
gis_dir <- "/Users/hoecker/Library/CloudStorage/GoogleDrive-tyler.hoecker@vibrantplanet.net/My Drive/GIS_Data/"

# Create validation set ------------------------------------------------------------------
# From the West-wide climate data, sample a subset of points to run
# validation on.
# Create a n*n tile grid across the study area
# tile_grid <- st_make_grid(ann_fut[[1]], n = 300)
# Sample n of these tiles
# tile_grid_samp <- vect(sample(tile_grid, 20))
huc12s <- st_read(here(gis_dir, "huc12/huc12_western_states_good.gpkg"))
random_hucs <- huc12s[sample(1:length(huc12s[["objectid"]]), size = 10), ] |>
  st_transform(crs = st_crs(ann_fut[[1]]))
# plot(random_hucs$geom)

# Move to script/function
# -------------------------------------------------------------------------------
# Read in climate data, crop, convert to data.table
# 2C annuals
ann_fut <- grep(
  "topoterra_2C_[0-9]{4}.tif",
  list.files(paste0(gis_dir, "TopoTerra")),
  value = TRUE
  # We have 31 years of data... drop 1985 to make it an even 30-yr normal
  )[2:31] |>
  paste0(climate_dir, "/", ... = _) |>
  map(\(x) rast(x))


# Mask each element of ann_fut to these sample tiles
fut_samp <- ann_fut |>
  map(\(x) mask(x, random_hucs))

# Convert each element into a data.table
fut_samp <- lapply(fut_samp, FUN = as.data.table, xy = TRUE)

# This is slow, save for convenience for now
saveRDS(fut_samp, "fut_samp.Rds")

# For contemporary validation, the focal cells will be historical normals, but we will
# build the covariance matrix based on the future annuals, because that works here
# Historical normals (downscaled from TerraClimate normals)
norm_hist <- grep(
  "topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
  list.files(climate_dir),
  value = TRUE
) |>
  paste0(climate_dir, "/", ... = _) |>
  rast() |>
  as.data.table(xy = TRUE)

# Create the same random sample of tiles for the historical data
hist_samp <- grep(
  "topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
  list.files(climate_dir),
  value = TRUE
) |>
  paste0(climate_dir, "/", ... = _) |>
  rast() |>
  mask(random_hucs) |>
  as.data.table(xy = TRUE)

# Remove things
rm(ann_fut, huc12s)
gc()

# Select contemporary analogs for validation set --------------------------------------------------
source("code/mahalanobis_D_fn.R")

# This takes 4.5 hours to 28011 points on 6 cores 
# library(profvis)
start <- Sys.time()
# profvis({
find_analogs(focal_data_cov = fut_samp,
             focal_data_mean = hist_samp,
             analog_data = norm_hist,
             var_names = c('aet', 'def', 'tmax', 'tmin'), 
             n_analog_pool = 5000000, 
             n_analog_use = 5000, 
             min_dist = 50000, 
             output_dir = "data/validation/5M-pool_5k-keep_50km-min",
             use_futures = TRUE,
             n_futures = 6)
# })
end <- Sys.time()
start - end

# Actual validation part starts here.... move to own script eventually
# Vegetation
lf_bps_220 <- terra::rast(paste0(gis_dir, "Landfire/lf_bps_west_220.tif")) |>
  terra::project("EPSG:3857", method = "mode")
# ------------------------------------------------------------------------------
val_datatable <- readRDS("data/validation/5M-pool_5k-keep_50km-min_2.Rds") 
val_datalist <- split(val_datatable, by = "focal_id")

# Quick explore of distance-sigma relationship
plot_dat <- sample(val_datatable, 10000) 

# Table of function arguments / validation parameters
val_params <- expand.grid(
  # "n_analog_pool" = c(5000000, 1000000, 500000, 10000), # HOPEFULLY ADD LATER
  "n_analog_use" = c(5000, 1000, 500),
  "dist_rule" = c(50000, 100000),
  "weighted" = c("TRUE", "FALSE"),
  "n_projections" = 1
  #"val_datatable" = val_datatable
)

val_wrapper <- function(
  n_analog_use,
  dist_rule,
  weighted,
  n_projections
){

  val_result <- val_datalist |>
    map_dfr(\(pt_i){

      # Evaluate validation parameters
      pt_i <- pt_i |>
        filter(dist_m > dist_rule) |>
        arrange(md) |>
        slice_head(n = n_analog_use) 

      veg_vote_fn(
        analog_data = pt_i,
        impact_data = lf_bps_220,
        n_analog_use,
        weighted,
        n_projections
      )
      }) |>
    mutate(n_analog_use = n_analog_use,
           dist_rule = dist_rule,
           weighted = weighted)
}

validation_sets <- purrr::pmap_dfr(val_params, val_wrapper, .progress = TRUE)

chance = 1/length(unique(validation_sets$observed))

validation_stats <- validation_sets |>
  filter(!is.na(aim), !is.na(observed)) |>
  mutate(correct = ifelse(aim == observed, 1, 0)) |>
  group_by(aim) |>
  mutate(occurrence = n()) |>
  ungroup()


lf_levels = freq(lf_bps_220)[,"value"]


validation_total <- validation_stats |>
  #mutate(dist_rule = as.character(dist_rule/1000)) |>
  filter(n_analog_use == 1000, dist_rule == 50000, weighted == TRUE) |>
  group_by(n_analog_use, dist_rule, weighted) |>
  summarise(prop_agree = mean(correct),
            f1 = f1_score(aim, observed),
            kappa = (prop_agree-chance)/(1-chance),
            avg_min_dist = median(min_dist),
            mean_sigma = mean(mean_sigma),
            mean_n90 = mean(n_aim_90))

ggplot(validation_total) +
  geom_col(aes(x = interaction(weighted, n_analog_use), 
               y = kappa, 
               fill = dist_rule), 
               position = 'dodge') +
  theme_bw(base_size = 14)


# ------------------------------------------------------------------------------------------------
# --- 2
val_dt_2 <- readRDS("data/validation/5M-pool_5k-keep_50km-min_2.Rds") 
head(val_dt_2)
val_dt_2_l <- split(val_dt_2, by = "focal_id")

source("code/analog_impact_fn.R")

val_result_2 <- val_dt_2_l |>
  map_dfr(\(x) veg_vote_fn(
    x,
    impact_data = impact_data,
    n_analogs = 500,
    weighted = TRUE 
    ),
    .progress = TRUE
  )

val_stats <- val_result_2 |>
  filter(!is.na(bps_true), !is.na(bps1)) |>
  mutate(correct = ifelse(bps1 == bps_true, 1, 0)) |>
  group_by(bps_true) |>
  mutate(occurrence = n()) |>
  ungroup()

# Correct classifications = 0.404 with 1000 analogs; 0.38 with 100 analogs; ~0.3 higher with wt_vote
sum(val_stats$correct)/length(val_stats$correct)
# Random classification rate 0.014
1/length(unique(val_stats$bps_true))

#-------------------------------------------------------------------------------
# Rasterize and write out COGs
val_rasts <- list('bps_true','bps1','mean_sigma','n_bps_90') %>%  #'current_bps','bps1','bps2','bps3','win_prop_raw','n_analogs','n_bps','mean_sigma_100') %>%  
  map(\(band) {
    
    print(paste0('Rasterizing ',band,'...'))
    
    r <- val_result_2 |> 
      select(x, y, band)  |> 
      terra::rast(x = _, 
                  type = 'xyz', 
                  crs = crs(impact_data)
                 ) 
        
    writeRaster(r, 
                paste0("data/validation/",band,'.tif'), 
                datatype = "INT2S",
                gdal = c("PROJECTION=EPSG:3857",
                         "TILED=YES",
                         "BLOCKXSIZE=128",
                         "BLOCKYSIZE=128",
                         "OVERVIEW-RESAMPLING=NEAREST",
                         "COMPRESS=DEFLATE"),
                overwrite = TRUE)

    return(r)
  }) 

# intersect
random_hucs

# Feather is faster to read/write but is > 2x the size...
feather::write_feather(val_dt_1, "data/validation/5M-pool_5k-keep_50km-min_2.feather")
