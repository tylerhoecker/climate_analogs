# Read in the sample points - a data.table of focal locations and 1000 best analogs > 50 km away 
test_pts <- readRDS(here("data/validation/random_tiles.Rds")) 


# Turn into a spatial object, for plotting purposes
# Get CRS into from climate data
this_crs <- grep("topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
                  list.files(climate_dir),
                  value = TRUE) |>
  paste0(climate_dir, "/", ... = _) |>
  rast()

states <- st_read(here("data/western_states/western_states.shp")) |>
  st_transform(crs = st_crs(this_crs))

# Create a spatial points dataframe of the sampled focal locations
focal_pts_sf <- test_pts |>
  mutate(focal_id = as.numeric(focal_id)) |>
  group_by(focal_id) |>
  summarise(x = first(f_x),
            y = first(f_y)) |>
  st_as_sf(coords = c("x","y"), crs = st_crs(this_crs)) |>
  # Many don't fall in Western state boundaries... crop out for now
  st_intersection(states)

analog_pts <- test_pts |>
  mutate(focal_id = as.numeric(focal_id)) |>
  filter(focal_id %in% focal_pts_sf$focal_id) |>
  (\(x) split(x, x$focal_id))() 

analog_pts_sf <- analog_pts |>
  map(\(x) st_as_sf(x, coords = c("a_x","a_y"), crs = st_crs(this_crs)))

pt_ids <- c(1000,3000,5000,7000)
analog_ids <- names(analog_pts_sf[pt_ids])

focal_pts_sf_eg <- focal_pts_sf[focal_pts_sf$focal_id %in% analog_ids,] |>
  arrange(focal_id)

analog_pts_eg  <- analog_pts[analog_ids]
analog_pts_sf_eg <- analog_pts_sf[analog_ids]

# Map showing random focal locations their 1000 analogs
tm_shape(states) +
  # State borders and fill
  tm_borders(col = "white", lwd = 2) +
    tm_fill(col = "grey10") +
  # Focal set 1
  tm_shape(focal_pts_sf_eg[1,]) +
    tm_dots(size = 1, col = 'blue') +
    tm_shape(analog_pts_sf_eg[[1]]) +
    tm_dots(col = 'sigma', size = 0.25, alpha = 0.5, palette = "Blues") +
  # Focal set 2
  tm_shape(focal_pts_sf_eg[2,]) +
    tm_dots(size = 1, col = 'green4') +
    tm_shape(analog_pts_sf_eg[[2]]) +
    tm_dots(col = 'sigma', size = 0.25, alpha = 0.5, palette = "Greens") +
  # Focal set 3
  tm_shape(focal_pts_sf_eg[3,]) +
    tm_dots(size = 1, col = 'red2') +
    tm_shape(analog_pts_sf_eg[[3]]) +
    tm_dots(col = 'sigma', size = 0.25, alpha = 0.5, palette = "Reds") +
  # Focal set 4
  tm_shape(focal_pts_sf_eg[4,]) +
    tm_dots(size = 1, col = 'purple2') +
    tm_shape(analog_pts_sf_eg[[4]]) +
    tm_dots(col = 'sigma', size = 0.25, alpha = 0.5, palette = "Purples") 

# Extract climate data at points
test_ann <- ann_fut |>
  map(\(x) terra::extract(x, focal_pts_sf_eg))

test_ann_pca <- seq_along(analog_ids) |>
  map(\(id){
    test_anns_i <- map_dfr(test_ann, ~ .x[id,])
    test_pca <- prcomp(test_anns_i[,-1])
    test_pca_ann <- test_pca$x
    return(as.data.frame(test_pca_ann))

    test_ann_norm_hist <- analog_pts_eg[[id]] |>
      rename(x = a_x, y = a_y) |>
      left_join(x = _, y = norm_hist) 

    test_norm_hist <- analog_pts_eg[[id]][1] |>
      rename(x = f_x, y = f_y) |>
      left_join(x = _, y = norm_hist) 

  } )

ggplot() +
  geom_point(data = test_norm_hist, aes(x = aet, y = def), size = 4) +
  geom_point(data = test_ann_norm_hist, aes(x = aet, y = def, color = sigma)) +
  geom_point(data = test_anns_i, aes(x = aet, y = def)) +
  theme_bw(base_size = 14)
