library(tidyverse)
library(rlang)
library(terra)
library(RColorBrewer)
library(viridisLite)
source("code/src/climate_analogs_R/veg_vote_fn.R")
template_r <- rast("data/climate/TopoTerra_1961-1990.tif", lyrs = "aet")
row_cells <- ceiling(dim(template_r)[1] / 5)
col_cells <- ceiling(dim(template_r)[2] / 4)
tile_extents <- getTileExtents(template_r, c(row_cells, col_cells))

##################
# PROCESS FUTURE PROJECTIONS
##################
BPS <- rast("data/veg_data/LC20_BPS_220.tif")
BPS <- crop(BPS, project(template_r, crs(BPS)))
states <- c("washington", "oregon", "california", "idaho", "montana", "nevada", "utah", "wyoming", "colorado", "new_mexico", "arizona") %>%
    str_to_title()
categorical_names <- c(
    "BPS_name_predicted"
)
numerical_names <- c(
    "sigma_score_BPS", "sigma_score_forest", "mean_sigma_all", "minimum_distance_predicted"
)




input_files <- list.files("data/reverse_analogs/outputs", pattern = "\\.csv\\.gz", full.names = TRUE)
# input_files <- input_files[!str_detect(input_files, paste(str_to_lower(states), collapse = "|"))]
for (file in input_files) {
    tile_name <- str_extract(file, "(?<=outputs/reverse_analogs_)[^\\.]+") %>% str_to_lower()
    local_name <- str_extract(tile_name, "^[a-z0-9]+(?:_[a-z0-9]+)*") %>%
        str_replace_all("_", " ") %>%
        str_to_title() %>%
        str_remove(" Tile\\d*$")
    # check if the local_name matches states
    if (!(local_name %in% states)) {
        # if it does, grab the extent of the tile
        local_border <- tile_extents[as.numeric(local_name), ] |>
            ext() |>
            vect(crs = crs(template_r))
    } else {
        local_border <- vect("data/western_states/western_states.shp") |>
            subset(NAME == local_name, NSE = TRUE) |>
            project("EPSG:4326")
    }

    # crop the BPS to the local border with a 500km buffer and project to the template raster
    BPS_local <- crop(BPS, terra::buffer(project(local_border, crs(BPS)), 500 * 1000))
    BPS_local <- project(BPS_local, template_r, method = "near", threads = TRUE)

    # clean outputs of climate analogs for the tile
    if (file.exists(paste0("data/reverse_analogs/outputs/cleaned_data/", tile_name, "_cleaned.csv.gz"))) {
        cleaned_local <- fread(paste0("data/reverse_analogs/outputs/cleaned_data/", tile_name, "_cleaned.csv.gz"))
    } else {
        local_data <- fread(file)
        cleaned_local <- setup_veg_prediction_bps(local_data, local_border, template_r, BPS_local)
        fwrite(cleaned_local, paste0("data/reverse_analogs/outputs/cleaned_data/", tile_name, "_cleaned.csv.gz"), compress = "gzip")
    }
    gc()
    rm(local_data)
    # compute the predictions using the various methods


    # calculate the climate-weighted vote winner for the BPS and Forest_NonForest predictions
    local_sigma_BPS <- calculate_top_sigma(cleaned_local, "BPS_name_actual", "BPS_name_predicted")
    local_sigma_FN <- calculate_top_sigma(cleaned_local, "Forest_NonForest_actual", "Forest_NonForest_predicted")
    reserve <- local_sigma_BPS[, .(BPS_name_actual, BPS_name_predicted, sigma_score)][, sigma_score_BPS := sigma_score][, sigma_score := NULL]


    local_sigma <- cbind(local_sigma_FN, reserve)[, sigma_score_forest := sigma_score][, sigma_score := NULL]

    # find the mean sigma score for all focal pixels
    mean_sigma_all <- cleaned_local[, .(mean_sigma_all = mean(sigma, na.rm = TRUE)), by = .(f_x, f_y)]
    # merge the mean sigma score to the local_sigma dataframe
    local_sigma <- local_sigma[mean_sigma_all, on = .(f_x, f_y)]
    rm(mean_sigma_all)

    # find the minimum distance to the predicted BPS
    ## build filter data.table
    filter_data <- local_sigma[, .(f_x, f_y, BPS_name_predicted)]
    ## filter by filter_data
    filtered_for_dist <- cleaned_local[filter_data, on = .(f_x, f_y, BPS_name_predicted)]
    min_dist <- filtered_for_dist[, .(minimum_distance_predicted = min(dist_km, na.rm = TRUE)), by = .(f_x, f_y)]

    local_sigma <- local_sigma[min_dist, on = .(f_x, f_y)]

    rm(filtered_for_dist, min_dist)


    # create rasters of predicted and observed vegetation using the climate-weighted method
    template_local <- crop(template_r, local_border)
    local_df <- local_sigma
    local_sv <- vect(local_df, geom = c("f_x", "f_y"), crs = "EPSG:4326")
    gc()
    # list of fiels to rasterize
    fields <- c(
        "BPS_name_predicted", "Forest_NonForest_predicted", "sigma_score_BPS", "sigma_score_forest",
        "mean_sigma_all", "minimum_distance_predicted"
    )
    # rasterize the fields to the template raster over a loop
    temp_list <- list()
    forest_nonforest_colors <- data.frame(value = c(1, 2), color = c("darkgreen", "orange"))
    for (i in seq_along(fields)) {
        field <- fields[i]
        if (str_detect(field, "BPS_CODE")) {
            temporary_rast <- rasterize(local_sv, template_local, field, na.rm = TRUE)
        } else if (str_detect(field, "sigma|distance")) {
            temporary_rast <- rasterize(local_sv, template_local, field)
        } else {
            temporary_rast <- rasterize_predicted_veg(local_df, template_local, field)
            if (str_detect(field, "Forest")) {
                temporary_rast <- as.factor(temporary_rast)
                coltab(temporary_rast) <- forest_nonforest_colors
                temporary_rast <- categories(temporary_rast, value = data.frame(ID = c(1, 2), category = c("Correct", "Incorrect")))
            } else {
                temporary_rast <- reclassify_cpal(temporary_rast)
            }
        }
        temp_list[[i]] <- temporary_rast
    }

    # combine the rasters into a single stack
    local_stack <- rast(temp_list)
    names(local_stack) <- fields
    # write the stack to a file
    # split out the categorical and numeric rasters
    # write the rasters to the outputs folder
    writeRaster(local_stack[[categorical_names]], paste0("data/reverse_analogs/outputs/rasters/t", tile_name, "_reverse_analogs_categorical.tif"), overwrite = TRUE)
    writeRaster(local_stack[[numerical_names]], paste0("data/reverse_analogs/outputs/rasters/t", tile_name, "_reverse_analogs_numerical.tif"), overwrite = TRUE)

    rm(
        local_stack, temp_list, local_sv, local_df, cleaned_local, local_sigma, local_sigma_FN, local_sigma_BPS, reserve,
        BPS_local
    )
    gc()
}
# burn in water, barren, etc. using the burnin mask and merge rasters together
# read in all rasters as collection
burnin_mask <- rast("data/veg_data/BPS_burnin_mask.tif")


fields <- c(
    "BPS_name_predicted",
    "sigma_score_BPS",
    "mean_sigma_all",
    "minimum_distance_predicted"
)

rasters <- list.files("data/reverse_analogs/outputs/rasters", full.names = TRUE, pattern = ".tif$")
# rasters <- rasters[!str_detect(rasters, paste(str_to_lower(states), collapse = "|"))]
rasters <- rasters[!str_detect(rasters, paste0(fields, collapse = "|"))]
rasters <- rasters[str_detect(rasters, "[0-9]")]
categorical_rasters <- rasters[str_detect(rasters, "categorical")]
numerical_rasters <- rasters[str_detect(rasters, "numerical")]
raster_collection <- vector(length = length(fields), mode = "list")
# loop through fields and merge rasters together. reclassify categorical rasters as needed and burn in the burnin mask
for (i in seq_along(fields)) {
    field <- fields[i]
    if (field %in% c("BPS_name_predicted", "BPS_name_actual")) {
        raster_collection[[i]] <- lapply(categorical_rasters, rast, lyrs = field) %>%
            lapply(reclassify_cpal) %>%
            sprc() %>%
            terra::merge() %>%
            terra::merge(burnin_mask, first = FALSE)
        names(raster_collection[[i]]) <- field
    } else if (str_detect(field, "BPS_CODE")) {
        raster_collection[[i]] <- lapply(numerical_rasters, rast, lyrs = field) %>%
            sprc() %>%
            terra::merge()
        names(raster_collection[[i]]) <- field
    } else {
        raster_collection[[i]] <- lapply(numerical_rasters, rast, lyrs = field) %>%
            sprc() %>%
            terra::merge()
        names(raster_collection[[i]]) <- field
    }
}

# write each raster out
for (i in seq_along(fields)) {
    field <- fields[i]
    writeRaster(raster_collection[[i]], paste0("data/reverse_analogs/outputs/rasters/", field, "_reverse_analogs.tif"), overwrite = TRUE)
}
