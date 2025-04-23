library(tidyverse)
library(rlang)
library(terra)
library(RColorBrewer)
library(viridisLite)
source("code/src/climate_analogs_R/veg_vote_fn.R")
template_r <- rast("data/climate/topoterra_hist_1961-1990.tif", lyrs = "aet")
row_cells <- ceiling(dim(template_r)[1] / 5)
col_cells <- ceiling(dim(template_r)[2] / 4)
tile_extents <- getTileExtents(template_r, c(row_cells, col_cells))

BPS <- rast("data/veg_data/LC20_BPS_220.tif")
BPS <- crop(BPS, project(template_r, crs(BPS)))
states <- c("washington", "oregon", "california", "idaho", "montana", "nevada", "utah", "wyoming", "colorado", "new_mexico", "arizona") %>%
    str_to_title()
categorical_names <- c("BPS_name_actual", "Forest_NonForest_actual", "BPS_name_predicted", "Forest_NonForest_predicted")
numerical_names <- c("sigma_score_BPS", "sigma_score_forest", "mean_sigma_all", "minimum_distance_predicted")



input_files <- list.files("data/validation/outputs", pattern = "\\.csv\\.gz", full.names = TRUE)
input_files <- input_files[!str_detect(input_files, paste(str_to_lower(states), collapse = "|"))]
for (file in input_files) {
    tile_name <- str_extract(file, "(?<=outputs/validation_)[^\\.]+") %>% str_to_lower()
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
    BPS_local <- crop(BPS, buffer(project(local_border, crs(BPS)), 500 * 1000))
    BPS_local <- project(BPS_local, template_r, method = "near", threads = TRUE)
    gc()
    # if (file.exists(paste0("data/validation/outputs/cleaned_data/", tile_name, "_cleaned.csv.gz"))) {
    #     cleaned_local <- fread(paste0("data/validation/outputs/cleaned_data/", tile_name, "_cleaned.csv.gz"))
    # } else {

    cleaned_local <- setup_veg_prediction_bps(fread(file), local_border, template_r, BPS_local)
    fwrite(cleaned_local, paste0("data/validation/outputs/cleaned_data/", tile_name, "_cleaned.csv.gz"), append = FALSE, compress = "gzip")
    # }
    gc()

    # compute the predictions using the various methods


    local_sigma_BPS <- calculate_top_sigma(cleaned_local, "BPS_name_actual", "BPS_name_predicted")
    local_sigma_FN <- calculate_top_sigma(cleaned_local, "Forest_NonForest_actual", "Forest_NonForest_predicted")
    reserve <- local_sigma_BPS[, .(BPS_name_actual, BPS_name_predicted, sigma_score)][, sigma_score_BPS := sigma_score][, sigma_score := NULL]


    local_sigma <- cbind(local_sigma_FN, reserve)[, sigma_score_forest := sigma_score][, sigma_score := NULL]

    # find the mean sigma score for all focal pixels
    mean_sigma_all <- cleaned_local[, .(mean_sigma_all = mean(sigma)), by = .(f_x, f_y)]
    # merge the mean sigma score to the local_sigma dataframe
    local_sigma <- local_sigma[mean_sigma_all, on = .(f_x, f_y)]
    rm(mean_sigma_all)
    gc()

    # find the minimum distance to the predicted BPS
    ## build filter data.table
    filter_data <- local_sigma[, .(f_x, f_y, BPS_name_predicted)]
    ## filter by filter_data
    filtered_for_dist <- cleaned_local[filter_data, on = .(f_x, f_y, BPS_name_predicted)]
    min_dist <- filtered_for_dist[, .(minimum_distance_predicted = min(dist_km)), by = .(f_x, f_y)]

    local_sigma <- local_sigma[min_dist, on = .(f_x, f_y)]
    rm(min_dist, filtered_for_dist)

    gc()


    # create rasters of predicted and observed vegetation using the sigma method
    template_local <- crop(template_r, buffer(local_border, 500 * 1000))
    local_df <- local_sigma
    local_sv <- vect(local_df, geom = c("f_x", "f_y"), crs = "EPSG:4326")
    # list of fiels to rasterize
    fields <- c(
        "BPS_name_predicted", "BPS_name_actual",
        "Forest_NonForest_predicted", "Forest_NonForest_actual",
        "sigma_score_BPS", "sigma_score_forest",
        "mean_sigma_all", "minimum_distance_predicted"
    )
    # rasterize the fields to the template raster over a loop
    temp_list <- list()
    forest_nonforest_colors <- data.frame(value = c(1, 2), color = c("darkgreen", "orange"))
    for (i in seq_along(fields)) {
        field <- fields[i]
        if (field == "match_BPS" | field == "match_forested") {
            temporary_rast <- local_sv %>%
                rasterize(template_local, field, fun = "sum")
            temporary_rast <- as.factor(temporary_rast)
            temporary_rast <- categories(temporary_rast, value = data.frame(ID = c(0, 1), category = c("Incorrect", "Correct")))
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
    # write the rasters to the outputs folder
    writeRaster(local_stack[[categorical_names]], paste0("data/validation/outputs/rasters/", tile_name, "_validation_categorical.tif"), overwrite = TRUE)
    writeRaster(local_stack[[numerical_names]], paste0("data/validation/outputs/rasters/", tile_name, "_validation_numerical.tif"), overwrite = TRUE)
    rm(
        local_stack, temp_list, local_sv, local_df, cleaned_local, local_sigma, local_sigma_FN, local_sigma_BPS, reserve,
        BPS_local
    )
    gc()
}

burnin_mask <- rast("data/veg_data/BPS_burnin_mask.tif")

fields <- c(
    "BPS_name_predicted", "BPS_name_actual",
    "sigma_score_BPS",
    "mean_sigma_all",
    "minimum_distance_predicted"
)

rasters <- list.files("data/validation/outputs/rasters", full.names = TRUE, pattern = ".tif$")
rasters <- rasters[!str_detect(rasters, paste(str_to_lower(states), collapse = "|"))]
rasters <- rasters[!str_detect(rasters, paste0(fields, collapse = "|"))]
rasters <- rasters[str_detect(rasters, "[0-9]")]
categorical_rasters <- rasters[str_detect(rasters, "categorical")]
numerical_rasters <- rasters[str_detect(rasters, "numerical")]
raster_collection <- vector(length = length(fields), mode = "list")
for (i in seq_along(fields)) {
    field <- fields[i]
    if (field %in% c("BPS_name_predicted", "BPS_name_actual")) {
        raster_collection[[i]] <- lapply(categorical_rasters, rast, lyrs = field) %>%
            sprc() %>%
            mosaic(fun = "max") %>%
            terra::merge(., burnin_mask, first = FALSE) %>%
            reclassify_cpal()
        names(raster_collection[[i]]) <- field
    } else {
        raster_collection[[i]] <- lapply(numerical_rasters, rast, lyrs = field) %>%
            sprc() %>%
            mosaic(fun = "max")
        names(raster_collection[[i]]) <- field
    }
}

for (i in seq_along(fields)) {
    field <- fields[i]
    writeRaster(raster_collection[[i]], paste0("data/validation/outputs/rasters/", field, "_validation.tif"), overwrite = TRUE)
}

washington <- fread("data/validation/validation_washington.csv.gz")

washington_border <- vect("data/western_states/western_states.shp") |>
    subset(NAME == "Washington", NSE = TRUE) |>
    project("EPSG:4326")
template_r <- rast("data/climate/topoterra_hist_1961-1990.tif", lyrs = 1)
BPS <- rast("data/veg_data/LC20_BPS_220.tif")
BPS <- crop(BPS, project(template_r, crs(BPS)))
BPS <- project(BPS, crs(template_r), method = "near", threads = TRUE)
BPS <- resample(BPS, template_r, method = "near", threads = TRUE)

cleaned_wa <- setup_veg_prediction_bps(washington, washington_border, template_r, BPS)
# bitvectors
actual_include <- !cleaned_wa$BPS_name_actual %in% c("Pacific riparian systems", "Interior riparian")
predicted_include <- !cleaned_wa$BPS_name_predicted %in% c("Pacific riparian systems", "Interior riparian")
both_include <- actual_include & predicted_include
cleaned_wa <- cleaned_wa[both_include, ]
rm(washington)

# compute the predictions using the various methods
washington_plurality <- compute_accuracy_plurality_bps(cleaned_wa)
washington_plurality_wsigma <- compute_accuracy_plurality_bps_wsigma(cleaned_wa)
washington_sigma <- compute_accuracy_sigma_bps(cleaned_wa)
washington_distance <- compute_accuracy_distance_bps(cleaned_wa)
washington_combined <- compute_accuracy_combined_bps(cleaned_wa)

# create a list of the results
results_list <- list(
    washington_plurality = washington_plurality,
    washington_plurality_wsigma = washington_plurality_wsigma,
    washington_sigma = washington_sigma,
    washington_distance = washington_distance,
    washington_combined = washington_combined
)

# now combine accuracy results into a single data frame
accuracy_df <- build_accuracy_tables(results_list) %>%
    mutate(
        dataset = "validation",
        state = "Washington"
    )

## plot the accuracy metrics
BPS_accuracy_plot <- accuracy_df %>%
    ggplot(aes(x = Method, y = accuracy, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Prediction) +
    labs(
        title = "Accuracy of BPS predictions",
        x = "Dataset",
        y = "Accuracy",
        fill = "Method"
    ) +
    ylim(c(0, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/BPS_accuracy_plot.jpg", BPS_accuracy_plot, width = 10, height = 10)

# kappa
BPS_kappa_plot <- accuracy_df %>%
    ggplot(aes(x = Method, y = kappa, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Prediction) +
    labs(
        title = "Kappa of BPS predictions",
        x = "Dataset",
        y = "Kappa",
        fill = "Method"
    ) +
    ylim(c(0, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/BPS_kappa_plot.jpg", BPS_kappa_plot, width = 10, height = 10)

# create rasters of predicted and observed vegetation using the sigma method
template_r <- rast("data/climate/topoterra_hist_1961-1990.tif", 1) %>%
    crop(washington_border)
wa_df <- washington_sigma$simplified_dt %>%
    mutate(
        match_BPS = ifelse(BPS_name_predicted == BPS_name_actual, 1, 0),
        match_forested = ifelse(Forest_NonForest_predicted == Forest_NonForest_actual, 1, 0),
    )

wa_sv <- vect(wa_df, geom = c("f_x", "f_y"), crs = "EPSG:4326")
# list of fiels to rasterize
fields <- c(
    "match_BPS",
    "BPS_name_predicted", "BPS_name_actual",
    "sigma_score_BPS",
    "average_sigma",
    "minimum_distance"
)
# rasterize the fields to the template raster over a loop
temp_list <- list()
forest_nonforest_colors <- data.frame(value = c(1, 2), color = c("darkgreen", "orange"))
for (i in seq_along(fields)) {
    field <- fields[i]
    if (field == "match_BPS" | field == "match_forested") {
        temporary_rast <- wa_sv %>%
            rasterize(template_r, field, fun = "sum")
        temporary_rast <- as.factor(temporary_rast)
        temporary_rast <- categories(temporary_rast, value = data.frame(ID = c(0, 1), category = c("Incorrect", "Correct")))
    } else if (str_detect(field, "sigma")) {
        temporary_rast <- rasterize(wa_sv, template_r, field)
    } else {
        temporary_rast <- rasterize_predicted_veg(wa_df, template_r, field)
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
wa_stack <- rast(temp_list)
names(wa_stack) <- fields
# write the stack to a file
writeRaster(wa_stack, "figures/wa_validation_stack.tif", overwrite = TRUE, datatype = "INT2U")




# histograms of each statistic
mean_hist <- sigma_df %>%
    ggplot(aes(x = mean, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(count / sum(count)))) +
    facet_wrap(~dataset) +
    labs(
        title = "Mean sigma dissimilarity",
        x = expression(Log[10] ~ Mean ~ sigma),
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(trans = scales::log_trans(base = 10)) +
    theme_minimal()


sd_hist <- sigma_df %>%
    ggplot(aes(x = sd, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(count / sum(count)))) +
    facet_wrap(~dataset) +
    labs(
        title = "SD sigma dissimilarity",
        x = expression(Log[10] ~ SD ~ sigma),
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(trans = scales::log_trans(base = 10)) +
    theme_minimal()

ggsave("figures/washington/sd_hist.jpg", sd_hist, width = 10, height = 10)

min_hist <- sigma_df %>%
    ggplot(aes(x = min, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(count / sum(count)))) +
    facet_wrap(~dataset) +
    labs(
        title = "Min sigma dissimilarity",
        x = "Min sigma",
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()

ggsave("figures/washington/min_hist.jpg", min_hist, width = 10, height = 10)

# scatter plots of the factorials
# mean x sd
mean_x_sd <- sigma_df %>%
    ggplot(aes(x = mean, y = sd, color = dataset)) +
    geom_point(alpha = 0.25) +
    labs(
        title = "mean by sd sigma dissimilarity",
        x = "mean sigma",
        y = "sd sigma",
        color = "dataset"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    facet_wrap(~dataset)

ggsave("figures/washington/mean_x_sd.jpg", mean_x_sd, width = 10, height = 10)
# min x mean

min_x_mean <- sigma_df %>%
    ggplot(aes(x = min, y = mean, color = dataset)) +
    geom_point(alpha = 0.25) +
    labs(
        title = "min by mean sigma dissimilarity",
        x = "min sigma",
        y = "mean sigma",
        color = "dataset"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    facet_wrap(~dataset)

ggsave("figures/washington/min_x_mean.jpg", min_x_mean, width = 10, height = 10)

# min x sd

min_x_sd <- sigma_df %>%
    ggplot(aes(x = min, y = sd, color = dataset)) +
    geom_point(alpha = 0.25) +
    labs(
        title = "min by sd sigma dissimilarity",
        x = "min sigma",
        y = "sd sigma",
        color = "dataset"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    facet_wrap(~dataset)

ggsave("figures/washington/min_x_sd.jpg", min_x_sd, width = 10, height = 10)


# distance to best analog
distance_2_best <- wa_df_full %>%
    group_by(dataset, f_x, f_y) %>%
    slice_min(sigma, n = 1, with_ties = FALSE) %>%
    ungroup()

# visualize the distance to best analog
best_analog_hist <- distance_2_best %>%
    ggplot(aes(x = dist_km, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5) +
    facet_wrap(~dataset) +
    labs(
        title = "Distance to best analog",
        x = expression("Distance to best analog (pseudo-log"[10] * ")"),
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    theme_minimal()

ggsave("figures/washington/best_analog_hist.jpg", best_analog_hist, width = 10, height = 10)

# distance to closest analog

distance_2_closest <- wa_df_full %>%
    group_by(dataset, f_x, f_y) %>%
    slice_min(dist_km, n = 1, with_ties = TRUE) %>%
    ungroup()

# visualize the distance to closest analog
closest_analog_hist <- distance_2_closest %>%
    ggplot(aes(x = dist_km, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5) +
    facet_wrap(~dataset) +
    labs(
        title = "Distance to closest analog",
        x = expression("Distance to closest analog (pseudo-log"[10] * ")"),
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()

ggsave("figures/washington/closest_analog_hist.jpg", closest_analog_hist, width = 10, height = 10)
