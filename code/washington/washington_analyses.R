library(tidyverse)
library(caret)
library(rlang)
library(terra)
library(RColorBrewer)
library(viridisLite)
library(data.table)
source("code/src/climate_analogs_R/veg_vote_fn.R")
source("code/src/veg_vote_fn.R")


# Load the Washington state boundary
wa <- vect("data/western_states/western_states.shp") |>
    subset(NAME == "Washington", NSE = TRUE) |>
    project("EPSG:4326")


wa_5p_df <- fread("data/washington/500km_0min_5p.csv")[, dataset := "5p"]
wa_025p_df <- fread("data/washington/500km_0min_025p.csv")[, dataset := "025p"]
wa_001p_df <- fread("data/washington/500km_0min_001p.csv")[, dataset := "001p"]
wa_5p_1p_df <- fread("data/washington/500km_0min_5p_1p.csv")[, dataset := "5p_1p"]
wa_10p_1p_df <- fread("data/washington/500km_0min_10p_1p.csv")[, dataset := "10p_1p"]
wa_20p_1p_df <- fread("data/washington/500km_0min_20p_1p.csv")[, dataset := "20p_1p"]
wa_5p_0min_df <- fread("data/washington/500km_0min_5p_0min.csv")[, dataset := "5p_0min"]

wa_df_unfiltered <- rbindlist(list(wa_5p_df, wa_025p_df, wa_001p_df, wa_5p_1p_df, wa_10p_1p_df, wa_20p_1p_df, wa_5p_0min_df))
wa_sv_unfiltered <- wa_df_unfiltered %>%
    vect(geom = c("f_x", "f_y"), crs = "EPSG:4326")

wa_df <- wa_df_unfiltered[sigma < 2.0, ]
wa_sv <- unique(wa_df, by = c("f_x", "f_y"))
template_r <- rast("data/washington/template_washington.tif") %>%
    crop(., project(buffer(wa, 500 * 1000), crs(.)))


analogs_df <- unique(wa_df[, .(a_x, a_y)])

analogs_v <- vect(analogs_df, geom = c("a_x", "a_y"), crs = "EPSG:4326")

# Load the land cover data
BPS <- rast("data/veg_data/LC20_BPS_220.tif") %>%
    crop(., project(buffer(wa, 500 * 1000), crs(.))) %>%
    project(., crs(template_r), method = "near") %>%
    resample(., template_r, method = "near")
activeCat(BPS) <- "BPS_CODE"
# extract BPS and EVT to analogs
BPS_analogs <- extract(BPS, analogs_v,
    na.rm = TRUE,
    xy = TRUE
)



# bind  the extracted values to the analogs dataframe
analogs_df_full <- data.table::data.table(analogs_df,
    BPS_CODE = as.numeric(as.character(BPS_analogs$BPS_CODE)) # ,
    #    EVT_NAME = EVT_analogs$EVT_NAME,
    #    ESP_NAME = ESP_analogs$ESP_NAME
)

# join with the 5p dataframe
wa_df_full <- wa_df[analogs_df_full, on = .(a_x = a_x, a_y = a_y)][BPS_CODE > 31]
# left_join(analogs_df_full, by = c("a_x", "a_y"))

# adjust vegetation names to simplify the datasets --------
mapping_categorical <- fread("data/veg_data/mapping_categorical_v5.csv")[, .(BPS_CODE, Cluster_name)]
mapping_categorical <- unique(mapping_categorical, by = "BPS_CODE")
# rename cluster name to simplified name
mapping_categorical <- mapping_categorical[, simplified_name := Cluster_name][, Cluster_name := NULL]
# select(BPS_CODE, final_id, Cluster_name)
BPS_csv <- fread("data/veg_data/LF20_BPS_220.csv")[, ":="(
    Custom_simplified_name = map_chr(BPS_NAME, simplify_lf_names),
    Forest_NonForest = map_chr(BPS_NAME, forest_nonforest))][
    , .(BPS_CODE, Custom_simplified_name, Forest_NonForest)
]
BPS_csv <- unique(mapping_categorical[BPS_csv, on = .(BPS_CODE)], by = "BPS_CODE") # mutate(
#     Simplified_BPS_Name = map_chr(BPS_NAME, simplify_lf_names),
#     Forest_NonForest = map_chr(BPS_NAME, forest_nonforest)
# ) %>%
#     select(BPS_CODE, Simplified_BPS_Name, Forest_NonForest) %>%
#     # remove duplicates
#     left_join(mapping_categorical, by = c("BPS_CODE" = "BPS_CODE")) %>%
#     distinct(BPS_CODE, .keep_all = TRUE)
# plurality votes for each -------

join_fn <- function(input_df, input_column, input_csv, by_statement, count_column) {
    # Perform the join using the provided by_statement
    # result <- input_df %>%
    #     select(dataset, f_x, f_y, {{ input_column }}) %>%
    #     left_join(input_csv, by = by_statement) %>%
    #     select({{ input_column }}, {{ count_column }})
    selection <- input_df[, .(dataset, f_x, f_y, input_column), env = list(input_column = input_column)]
    left_join <- input_csv[selection, on = by_statement, env = list(by_statement = by_statement)][
        , .(input_column, count_column),
        env = list(input_column = input_column, count_column = count_column)
    ]

    return(left_join)
}
best_analogs_BPS <- join_fn(wa_df_full, "BPS_CODE", BPS_csv, "BPS_CODE", "simplified_name")
# join the plurality votes to the full dataframe
wa_df_full_plurality <- cbind(wa_df_full, best_analogs_BPS)[, -c("BPS_CODE")]
# wa_df_full_plurality <- wa_df_full %>%
#     cbind(c(best_analogs_BPS, best_analogs_EVT, best_analogs_ESP)) %>%
#     select(-BPS_CODE, -EVT_NAME, -ESP_NAME)


# build raster that contains the name and votes for the highest plurality for each dataset --------
# function to get the highest plurality for each dataset


BPS_analog_plurality <- calculate_top_vote(wa_df_full_plurality, "simplified_name", "BPS_vote", n = 1)
# filter to each dataset then rasterize
# filter function

BPS_analog_5p_veg <- vote_filter_fn(BPS_analog_plurality, "5p", "simplified_name", "BPS_vote")
BPS_analog_025p_veg <- vote_filter_fn(BPS_analog_plurality, "025p", "simplified_name", "BPS_vote")
BPS_analog_001p_veg <- vote_filter_fn(BPS_analog_plurality, "001p", "simplified_name", "BPS_vote")
BPS_analog_5p_1p_veg <- vote_filter_fn(BPS_analog_plurality, "5p_1p", "simplified_name", "BPS_vote")
BPS_analog_10p_1p_veg <- vote_filter_fn(BPS_analog_plurality, "10p_1p", "simplified_name", "BPS_vote")
BPS_analog_20p_1p_veg <- vote_filter_fn(BPS_analog_plurality, "20p_1p", "simplified_name", "BPS_vote")

# function that converts the names into values

BPS_analog_5p_veg_r <- rasterize_predicted_veg(BPS_analog_5p_veg, template_r, simplified_name)
BPS_analog_025p_veg_r <- rasterize_predicted_veg(BPS_analog_025p_veg, template_r, simplified_name)
BPS_analog_001p_veg_r <- rasterize_predicted_veg(BPS_analog_001p_veg, template_r, simplified_name)
BPS_analog_5p_1p_veg_r <- rasterize_predicted_veg(BPS_analog_5p_1p_veg, template_r, simplified_name)
BPS_analog_10p_1p_veg_r <- rasterize_predicted_veg(BPS_analog_10p_1p_veg, template_r, simplified_name)
BPS_analog_20p_1p_veg_r <- rasterize_predicted_veg(BPS_analog_20p_1p_veg, template_r, simplified_name)

BPS_analog_veg_r <- c(
    BPS_analog_5p_veg_r,
    BPS_analog_025p_veg_r,
    BPS_analog_001p_veg_r,
    BPS_analog_5p_1p_veg_r,
    BPS_analog_10p_1p_veg_r,
    BPS_analog_20p_1p_veg_r
) %>%
    crop(wa)
names(BPS_analog_veg_r) <- c("5p", "025p", "001p", "5p_1p", "10p_1p", "20p_1p")
# write out the rasters
writeRaster(BPS_analog_veg_r, "data/washington/BPS_analog_veg.tif", overwrite = TRUE)
# Now we can do the same, but looking at the vote counts for each dataset following the plurality maps above -------
# BPS

# filter to each dataset then rasterize

BPS_analog_5p_vote_r <- rasterize_predicted_veg(BPS_analog_5p_veg, template_r, BPS_vote)
BPS_analog_025p_vote_r <- rasterize_predicted_veg(BPS_analog_025p_veg, template_r, BPS_vote)
BPS_analog_001p_vote_r <- rasterize_predicted_veg(BPS_analog_001p_veg, template_r, BPS_vote)
BPS_analog_5p_1p_vote_r <- rasterize_predicted_veg(BPS_analog_5p_1p_veg, template_r, BPS_vote)
BPS_analog_10p_1p_vote_r <- rasterize_predicted_veg(BPS_analog_10p_1p_veg, template_r, BPS_vote)
BPS_analog_20p_1p_vote_r <- rasterize_predicted_veg(BPS_analog_20p_1p_veg, template_r, BPS_vote)

BPS_analog_vote_r <- c(
    BPS_analog_5p_vote_r,
    BPS_analog_025p_vote_r,
    BPS_analog_001p_vote_r,
    BPS_analog_5p_1p_vote_r,
    BPS_analog_10p_1p_vote_r,
    BPS_analog_20p_1p_vote_r
) %>%
    crop(wa)
names(BPS_analog_vote_r) <- c("5p", "025p", "001p", "5p_1p", "10p_1p", "20p_1p")

# write out the rasters
writeRaster(BPS_analog_vote_r, "data/washington/BPS_analog_vote.tif", overwrite = TRUE)
# compare actual and predicted vegetation
# BPS
## take all focal points and extract the BPS code
focal_points <- unique(wa_df_full[, .(f_x, f_y)]) %>%
    vect(geom = c("f_x", "f_y"), crs = "EPSG:4326")
# focal_points <- wa_df_full %>%
#     select(f_x, f_y) %>%
#     distinct(f_x, f_y) %>%
#     vect(geom = c("f_x", "f_y"), crs = "EPSG:4326")

BPS_actual_extracted <- extract(BPS, focal_points, xy = FALSE, ID = FALSE) %>%
    pull()

# join to the CSV
BPS_actual_df <- focal_points %>%
    as.data.frame(geom = "XY") %>%
    as.data.table()
BPS_actual_df <- BPS_actual_df[, BPS_actual := as.numeric(as.character(BPS_actual_extracted))]
BPS_actual_df <- BPS_csv[BPS_actual_df, on = .(BPS_CODE = BPS_actual)][BPS_CODE > 31]
# rename
BPS_actual_df <- BPS_actual_df[, .(BPS_name_actual = simplified_name, Forest_NonForest_actual = Forest_NonForest, f_x = x, f_y = y)]

# BPS_actual_df <- focal_points %>%
#     as.data.frame(geom = "XY") %>%
#     mutate(BPS_actual = as.numeric(as.character(BPS_actual_extracted))) %>%
#     left_join(BPS_csv, by = c("BPS_actual" = "BPS_CODE")) %>%
#     rename(
#         BPS_name_actual = simplified_name,
#         Forest_NonForest_actual = Forest_NonForest,
#         f_x = x,
#         f_y = y
#     )

# join to the plurality dataframe
BPS_predicted_df <- BPS_analog_plurality[, .(f_x, f_y, dataset, simplified_name)]
# rename to BPS_predicted
BPS_predicted_df <- BPS_predicted_df[, BPS_name_predicted := simplified_name][, simplified_name := NULL]
BPS_predicted_df <- BPS_actual_df[BPS_predicted_df, on = .(f_x, f_y)]
BPS_predicted_df <- BPS_predicted_df[complete.cases(BPS_predicted_df), ]

# BPS_predicted_df <- BPS_analog_plurality %>%
#     select(f_x, f_y, dataset, simplified_name) %>%
#     left_join(BPS_actual_df, by = c("f_x", "f_y")) %>%
#     drop_na()

## binary comparison, if match 1, else 0 then summarize
BPS_accuracy <- build_accuracy_stats_wclass(BPS_predicted_df, "BPS_name_actual", "BPS_name_predicted", "dataset")

# assess the accuracy of detecting forest vs non-forest
# BPS
## compute whether a pixel's most voted analog is forest or non-forest
BPS_predicted_df <- BPS_predicted_df[, Forest_NonForest_predicted := map_chr(BPS_name_predicted, forest_nonforest)]
# BPS_predicted_df <- BPS_predicted_df %>%
#     mutate(Forest_NonForest_predicted = map_chr(Simplified_BPS_Name, forest_nonforest))

# kappa for forested
BPS_accuracy_forest <- build_accuracy_stats_wclass(BPS_predicted_df, "Forest_NonForest_actual", "Forest_NonForest_predicted", "dataset")


# alternative voting scheme. instead of most votes, let's do a score based system
# Rescale sigma to 0-1, low sigmas being 1 and high sigmas being 0. Add up the scores for each analog and take the highest score
# BPS
## rescale sigma
rescaled_wa_df_full <- wa_df_full[, sigma_rescaled := 1 - (sigma - min(sigma)) / (max(sigma) - min(sigma))]
# now group by dataset and focal point and sum the rescaled sigma for each BPS_CODE,
# then take the highest score for each focal point
BPS_analog_scores <- rescaled_wa_df_full[, .(score = sum(sigma_rescaled)), by = .(dataset, f_x, f_y, BPS_CODE)]

BPS_analog_scores <- BPS_analog_scores[BPS_analog_scores[, .I[which.max(score)], by = .(dataset, f_x, f_y)]$V1]
# now evalueate the accuracy of the score based system
BPS_predicted_df <- BPS_csv[BPS_analog_scores, on = "BPS_CODE"]
BPS_predicted_df <- BPS_predicted_df[, ":="(BPS_predicted = simplified_name, simplified_name = NULL)]
BPS_predicted_df <- BPS_actual_df[BPS_predicted_df, on = .(f_x, f_y)]
BPS_predicted_df <- BPS_predicted_df[complete.cases(BPS_predicted_df), ]
# calculate kappa
BPS_accuracy_scores_sigma <- build_accuracy_stats_wclass(BPS_predicted_df, "BPS_name_actual", "BPS_predicted", "dataset")

BPS_predicted_df <- BPS_predicted_df[, Forest_NonForest_predicted := map_chr(BPS_predicted, forest_nonforest)]
BPS_forest_match_score <- BPS_predicted_df[, Forest_match := as.numeric(Forest_NonForest_actual == Forest_NonForest_predicted)][, .(match = mean(Forest_match)), by = dataset]

# kappa for forested
BPS_accuracy_forest_scores <- build_accuracy_stats_wclass(BPS_predicted_df, "Forest_NonForest_actual", "Forest_NonForest_predicted", "dataset")


# rate by distance weighted best analogs
# BPS
## add up distances for each dataset focal pixel and BPS_CODE
BPS_analog_distance <- wa_df_full[, .(distance_rescaled = 1 - (dist_km - min(dist_km)) / (max(dist_km) - min(dist_km))), by = .(dataset, f_x, f_y, BPS_CODE)]
BPS_analog_distances <- BPS_analog_distance[, .(distance = sum(distance_rescaled)), by = .(dataset, f_x, f_y, BPS_CODE)]
# now take the best distance for each focal pixel
BPS_analog_distances <- BPS_analog_distances[BPS_analog_distances[, .I[which.max(distance)], by = .(dataset, f_x, f_y)]$V1]

BPS_predicted_df <- BPS_csv[BPS_analog_distances, on = "BPS_CODE"]
BPS_predicted_df <- BPS_predicted_df[, ":="(BPS_predicted = simplified_name, simplified_name = NULL)]
# join to the actual data
BPS_predicted_df <- BPS_actual_df[BPS_predicted_df, on = .(f_x, f_y)]
BPS_predicted_df <- BPS_predicted_df[complete.cases(BPS_predicted_df), ]


distance_check <- BPS_predicted_df[, BPS_match := as.factor(BPS_name_actual == simplified_name)][, .(sigmas = mean(sigma), distances = mean(dist_km)), by = .(dataset, BPS_match)]


# calculate kappa
BPS_accuracy_distances <- build_accuracy_stats_wclass(BPS_predicted_df, "BPS_name_actual", "BPS_predicted", "dataset")

BPS_predicted_df <- BPS_predicted_df[, Forest_NonForest_predicted := map_chr(BPS_predicted, forest_nonforest)]

# kappa for forested
BPS_accuracy_forest_distances <- build_accuracy_stats_wclass(BPS_predicted_df, "Forest_NonForest_actual", "Forest_NonForest_predicted", "dataset")

# combine sigma and distance
# BPS
## rescale distances and sigma
rescaled_wa_df_full <- wa_df_full[, sigma_rescaled := 1 - (sigma - min(sigma)) / (max(sigma) - min(sigma))][, distance_rescaled := 1 - (dist_km - min(dist_km)) / (max(dist_km) - min(dist_km))]
# now group by dataset and focal point and compute the euclidean distance between the rescaled sigma and distance
# then take the lowest value
BPS_analog_combined <- rescaled_wa_df_full[, .(combined = sqrt(sigma_rescaled^2 + distance_rescaled^2)), by = .(dataset, f_x, f_y, BPS_CODE)]
BPS_analog_combined <- BPS_analog_combined[BPS_analog_combined[, .I[which.max(combined)], by = .(dataset, f_x, f_y)]$V1]

BPS_predicted_df <- BPS_csv[BPS_analog_combined, on = "BPS_CODE"]
BPS_predicted_df <- BPS_predicted_df[, ":="(BPS_predicted = simplified_name, simplified_name = NULL)]
BPS_predicted_df <- BPS_actual_df[BPS_predicted_df, on = .(f_x, f_y)]
BPS_predicted_df <- BPS_predicted_df[complete.cases(BPS_predicted_df), ]


# calculate kappa
BPS_accuracy_combined <- build_accuracy_stats_wclass(BPS_predicted_df, "BPS_name_actual", "BPS_predicted", "dataset")

BPS_predicted_df <- BPS_predicted_df[, Forest_NonForest_predicted := map_chr(BPS_predicted, forest_nonforest)]

# kappa for forested
BPS_accuracy_forest_combined <- build_accuracy_stats_wclass(BPS_predicted_df, "Forest_NonForest_actual", "Forest_NonForest_predicted", "dataset")

# combine all the results into a single dataframe labeling the method used
BPS_results <- rbind(
    BPS_accuracy %>% mutate(method = "Plurality Vote", prediction = "Simplified BPS"),
    BPS_accuracy_forest %>% mutate(method = "Plurality Vote", prediction = "Forested/nonForested"),
    BPS_accuracy_scores_sigma %>% mutate(method = "Best cumulative sigma", prediction = "Simplified BPS"),
    BPS_accuracy_forest_scores %>% mutate(method = "Best cumulative sigma", prediction = "Forested/nonForested"),
    BPS_accuracy_distances %>% mutate(method = "Highest cumulative distance 'score'", prediction = "Simplified BPS"),
    BPS_accuracy_forest_distances %>% mutate(method = "Highest cumulative distance 'score'", prediction = "Forested/nonForested"),
    BPS_accuracy_combined %>% mutate(method = "Highest cumulative combined 'score'", prediction = "Simplified BPS"),
    BPS_accuracy_forest_combined %>% mutate(method = "Highest cumulative combined 'score'", prediction = "Forested/nonForested")
) %>%
    mutate(accuracy = as.numeric(accuracy), kappa = as.numeric(kappa))
write_csv(BPS_results, "data/washington/BPS_predictions.csv")

# do same for kappa results
# visualize accuracy metrics
# BPS
## plot the accuracy metrics
BPS_accuracy_plot <- BPS_results %>%
    ggplot(aes(x = class, y = accuracy, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~prediction) +
    labs(
        title = "Accuracy of BPS predictions",
        x = "Dataset",
        y = "Accuracy",
        fill = "Method"
    ) +
    ylim(c(0, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw()
ggsave("figures/washington/BPS_accuracy_plot.jpg", BPS_accuracy_plot, width = 10, height = 10)

# kappa
BPS_kappa_plot <- BPS_results %>%
    ggplot(aes(x = class, y = kappa, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~prediction) +
    labs(
        title = "Kappa of BPS predictions",
        x = "Dataset",
        y = "Kappa",
        fill = "Method"
    ) +
    ylim(c(0, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw()
ggsave("figures/washington/BPS_kappa_plot.jpg", BPS_kappa_plot, width = 10, height = 10)

# sigma dissimilarity
sigma_df <- wa_df_full[, .(mean = mean(sigma), sd = sd(sigma), min = min(sigma)), by = .(dataset, f_x, f_y)]
# sigma_df <- wa_df_full %>%
#     group_by(dataset, f_x, f_y) %>%
#     summarise(mean = mean(sigma), sd = sd(sigma), min = min(sigma)) %>%
#     ungroup()


# visualize the sigma results

# histograms of each statistic
mean_hist <- sigma_df %>%
    ggplot(aes(x = mean, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(width * density))) +
    facet_wrap(~dataset) +
    labs(
        title = "Mean sigma dissimilarity",
        x = expression(Mean ~ sigma),
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_continuous(trans = scales::log1p_trans()) +
    theme_bw()

ggsave("figures/washington/mean_hist.jpg", mean_hist, width = 10, height = 10)

sd_hist <- sigma_df %>%
    ggplot(aes(x = sd, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(width * density))) +
    facet_wrap(~dataset) +
    labs(
        title = "SD sigma dissimilarity",
        x = expression(SD ~ sigma),
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_continuous(trans = scales::log_trans(base = 10)) +
    theme_bw()

ggsave("figures/washington/sd_hist.jpg", sd_hist, width = 10, height = 10)

min_hist <- sigma_df %>%
    ggplot(aes(x = min, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(width * density))) +
    facet_wrap(~dataset) +
    labs(
        title = "Min sigma dissimilarity",
        x = "Min sigma",
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_continuous(trans = scales::pseudo_log_trans(10)) +
    theme_bw()

ggsave("figures/washington/min_hist.jpg", min_hist, width = 10, height = 10)

# histogram from the best cumulative sigma score
BPS_analog_scores_hist <- BPS_analog_scores %>%
    ggplot(aes(x = score, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity", aes(y = after_stat(width * density))) +
    facet_wrap(~dataset) +
    labs(
        title = "Best cumulative sigma score",
        x = "Best cumulative sigma score",
        y = "Relative Frequency",
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw()
ggsave("figures/washington/BPS_analog_scores_hist.jpg", BPS_analog_scores_hist, width = 10, height = 10)

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
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
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
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
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
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    facet_wrap(~dataset)

ggsave("figures/washington/min_x_sd.jpg", min_x_sd, width = 10, height = 10)


# distance to best analog
dist_2_best <- wa_df_full[wa_df_full[, .I[which.min(sigma)], by = .(dataset, f_x, f_y)]$V1]
# distance_2_best <- wa_df_full %>%
#     group_by(dataset, f_x, f_y) %>%
#     slice_min(sigma, n = 1, with_ties = FALSE) %>%
#     ungroup()

# visualize the distance to best analog
best_analog_hist <- dist_2_best %>%
    group_by(dataset) %>%
    ggplot(aes(x = dist_km, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, aes(y = after_stat(width * density))) +
    facet_wrap(~dataset) +
    labs(
        title = "Distance to best analog",
        x = expression("Distance to best analog (km)"),
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_continuous(trans = scales::log1p_trans()) +
    theme_bw()

ggsave("figures/washington/best_analog_hist.jpg", best_analog_hist, width = 10, height = 10)

# distance to closest analog

distance_2_closest <- wa_df_full[wa_df_full[, .I[which.min(dist_km)], by = .(dataset, f_x, f_y)]$V1]
# distance_2_closest <- wa_df_full %>%
#     group_by(dataset, f_x, f_y) %>%
#     slice_min(dist_km, n = 1, with_ties = TRUE) %>%
#     ungroup()

# visualize the distance to closest analog
closest_analog_hist <- distance_2_closest %>%
    ggplot(aes(x = dist_km, fill = dataset)) +
    geom_histogram(bins = 100, alpha = 0.5, aes(y = after_stat(width * density))) +
    facet_wrap(~dataset) +
    labs(
        title = "Distance to closest analog",
        x = expression("Distance to closest analog (km)"),
        fill = "dataset"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_continuous(trans = scales::log_trans(10)) +
    theme_bw()

ggsave("figures/washington/closest_analog_hist.jpg", closest_analog_hist, width = 10, height = 10)



# other datasets ------
EVT <- rast("data/veg_data/LC23_EVT_240.tif") %>%
    crop(., project(buffer(wa, 500 * 1000), crs(.))) %>%
    project(., crs(template_r), method = "near") %>%
    resample(., template_r, method = "near")
activeCat(EVT) <- "EVT_NAME"

ESP <- rast("data/veg_data/us_140esp.tif") %>%
    crop(., project(buffer(wa, 500 * 1000), crs(.))) %>%
    project(., crs(template_r), method = "near") %>%
    resample(., template_r, method = "near")
activeCat(ESP) <- "ESP_NAME"

EVT_analogs <- extract(EVT, analogs_v,
    na.rm = TRUE,
    xy = TRUE
)
ESP_analogs <- extract(ESP, analogs_v,
    na.rm = TRUE,
    xy = TRUE
)
analogs_df_full <- data.table::data.table(analogs_df,
    BPS_CODE = as.numeric(as.character(BPS_analogs$BPS_CODE)),
    EVT_NAME = EVT_analogs$EVT_NAME,
    ESP_NAME = ESP_analogs$ESP_NAME
)
EVT_csv <- fread("data/veg_data/LF23_EVT_240.csv")[, ":="(
    Simplified_EVT_Name = map_chr(EVT_NAME, simplify_lf_names),
    Forest_NonForest = map_chr(EVT_NAME, forest_nonforest)
)][
    , .(EVT_NAME, Simplified_EVT_Name, Forest_NonForest)
]
EVT_csv <- unique(EVT_csv, by = "EVT_NAME")

# EVT_csv <- read_csv("data/veg_data/LF23_EVT_240.csv") %>%
#     mutate(
#         Simplified_EVT_Name = map_chr(EVT_NAME, simplify_lf_names),
#         Forest_NonForest = map_chr(EVT_NAME, forest_nonforest)
#     ) %>%
#     select(EVT_NAME, Simplified_EVT_Name, Forest_NonForest) %>%
# distinct(EVT_NAME, .keep_all = TRUE)
ESP_csv <- fread("data/veg_data/us_140esp.csv")[, ":="(
    Simplified_ESP_Name = map_chr(ESP_Name, simplify_lf_names),
    Forest_NonForest = map_chr(ESP_Name, forest_nonforest)
)][
    , .(ESP_Name, Simplified_ESP_Name, Forest_NonForest)
][, ESP_NAME := ESP_Name][, ESP_Name := NULL]
ESP_csv <- unique(ESP_csv, by = "ESP_NAME")
# ESP_csv <- read_csv("data/veg_data/us_140esp.csv") %>%
#     mutate(
#         Simplified_ESP_Name = map_chr(ESP_Name, simplify_lf_names),
#         Forest_NonForest = map_chr(ESP_Name, forest_nonforest)
#     ) %>%
#     select(ESP_Name, Simplified_ESP_Name, Forest_NonForest) %>%
#     distinct(ESP_Name, .keep_all = TRUE)

best_analogs_EVT <- join_fn(wa_df_full, "EVT_NAME", EVT_csv, "EVT_NAME", "Simplified_EVT_Name")
best_analogs_ESP <- join_fn(wa_df_full, "ESP_NAME", ESP_csv, "ESP_NAME", "Simplified_ESP_Name")


wa_df_full_plurality <- cbind(wa_df_full, best_analogs_BPS, best_analogs_EVT, best_analogs_ESP)[, -c("BPS_CODE", "EVT_NAME", "ESP_NAME")]
# EVT
EVT_analog_plurality <- calculate_top_vote(wa_df_full_plurality, Simplified_EVT_Name, EVT_vote, n = 1)

# filter to each dataset then rasterize
EVT_analog_5p_veg <- vote_filter_fn(EVT_analog_plurality, "5p", Simplified_EVT_Name)
EVT_analog_025p_veg <- vote_filter_fn(EVT_analog_plurality, "025p", Simplified_EVT_Name)
EVT_analog_001p_veg <- vote_filter_fn(EVT_analog_plurality, "001p", Simplified_EVT_Name)

EVT_analog_5p_veg_r <- rasterize_predicted_veg(EVT_analog_5p_veg, template_r, Simplified_EVT_Name)
EVT_analog_025p_veg_r <- rasterize_predicted_veg(EVT_analog_025p_veg, template_r, Simplified_EVT_Name)
EVT_analog_001p_veg_r <- rasterize_predicted_veg(EVT_analog_001p_veg, template_r, Simplified_EVT_Name)

EVT_analog_veg_r <- c(EVT_analog_5p_veg_r, EVT_analog_025p_veg_r, EVT_analog_001p_veg_r) %>%
    crop(wa)
names(EVT_analog_veg_r) <- c("5p", "025p", "001p")

# ESP
ESP_analog_plurality <- calculate_top_vote(wa_df_full_plurality, Simplified_ESP_Name, ESP_vote, n = 1)

# filter to each dataset then rasterize
ESP_analog_5p_veg <- vote_filter_fn(ESP_analog_plurality, "5p", Simplified_ESP_Name)
ESP_analog_025p_veg <- vote_filter_fn(ESP_analog_plurality, "025p", Simplified_ESP_Name)
ESP_analog_001p_veg <- vote_filter_fn(ESP_analog_plurality, "001p", Simplified_ESP_Name)

ESP_analog_5p_veg_r <- rasterize_predicted_veg(ESP_analog_5p_veg, template_r, Simplified_ESP_Name)
ESP_analog_025p_veg_r <- rasterize_predicted_veg(ESP_analog_025p_veg, template_r, Simplified_ESP_Name)
ESP_analog_001p_veg_r <- rasterize_predicted_veg(ESP_analog_001p_veg, template_r, Simplified_ESP_Name)

ESP_analog_veg_r <- c(ESP_analog_5p_veg_r, ESP_analog_025p_veg_r, ESP_analog_001p_veg_r) %>%
    crop(wa)
names(ESP_analog_veg_r) <- c("5p", "025p", "001p")
writeRaster(EVT_analog_veg_r, "data/washington/EVT_analog_veg.tif", overwrite = TRUE)
writeRaster(ESP_analog_veg_r, "data/washington/ESP_analog_veg.tif", overwrite = TRUE)

# EVT
EVT_analog_5p_vote_r <- rasterize_predicted_veg(EVT_analog_5p_veg, template_r, EVT_vote)
EVT_analog_025p_vote_r <- rasterize_predicted_veg(EVT_analog_025p_veg, template_r, EVT_vote)
EVT_analog_001p_vote_r <- rasterize_predicted_veg(EVT_analog_001p_veg, template_r, EVT_vote)

EVT_analog_vote_r <- c(EVT_analog_5p_vote_r, EVT_analog_025p_vote_r, EVT_analog_001p_vote_r) %>%
    crop(wa)
names(EVT_analog_vote_r) <- c("5p", "025p", "001p")

# ESP
ESP_analog_5p_vote_r <- rasterize_predicted_veg(ESP_analog_5p_veg, template_r, ESP_vote)
ESP_analog_025p_vote_r <- rasterize_predicted_veg(ESP_analog_025p_veg, template_r, ESP_vote)
ESP_analog_001p_vote_r <- rasterize_predicted_veg(ESP_analog_001p_veg, template_r, ESP_vote)

ESP_analog_vote_r <- c(ESP_analog_5p_vote_r, ESP_analog_025p_vote_r, ESP_analog_001p_vote_r) %>%
    crop(wa)
names(ESP_analog_vote_r) <- c("5p", "025p", "001p")

writeRaster(EVT_analog_vote_r, "data/washington/EVT_analog_vote.tif", overwrite = TRUE)
writeRaster(ESP_analog_vote_r, "data/washington/ESP_analog_vote.tif", overwrite = TRUE)

# EVT
EVT_actual_extracted <- extract(EVT, focal_points, xy = FALSE, ID = FALSE) %>%
    pull()

# join to the CSV
EVT_actual_df <- focal_points %>%
    as.data.frame(geom = "XY") %>%
    mutate(EVT_actual = as.character(EVT_actual_extracted)) %>%
    left_join(EVT_csv, by = c("EVT_actual" = "EVT_NAME")) %>%
    rename(
        EVT_name_actual = Simplified_EVT_Name,
        Forest_NonForest_actual = Forest_NonForest,
        f_x = x,
        f_y = y
    )

# join to the plurality dataframe
EVT_predicted_df <- EVT_analog_plurality %>%
    select(f_x, f_y, dataset, Simplified_EVT_Name) %>%
    left_join(EVT_actual_df, by = c("f_x", "f_y"))

## binary comparison, if match 1, else 0 then summarize
EVT_prediction_strength <- EVT_predicted_df %>%
    mutate(EVT_match = as.numeric(EVT_name_actual == Simplified_EVT_Name)) %>%
    group_by(dataset) %>%
    summarize(match = mean(EVT_match), .groups = "drop")

# ESP
ESP_actual_extracted <- extract(ESP, focal_points, xy = FALSE, ID = FALSE) %>%
    pull()

# join to the CSV
ESP_actual_df <- focal_points %>%
    as.data.frame(geom = "XY") %>%
    mutate(ESP_actual = as.character(ESP_actual_extracted)) %>%
    left_join(ESP_csv, by = c("ESP_actual" = "ESP_Name")) %>%
    rename(
        ESP_name_actual = Simplified_ESP_Name,
        Forest_NonForest_actual = Forest_NonForest,
        f_x = x,
        f_y = y
    )

# join to the plurality dataframe
ESP_predicted_df <- ESP_analog_plurality %>%
    select(f_x, f_y, dataset, Simplified_ESP_Name) %>%
    left_join(ESP_actual_df, by = c("f_x", "f_y"))

## binary comparison, if match 1, else 0 then summarize
ESP_prediction_strength <- ESP_predicted_df %>%
    mutate(ESP_match = as.numeric(ESP_name_actual == Simplified_ESP_Name)) %>%
    group_by(dataset) %>%
    summarize(match = mean(ESP_match), .groups = "drop")

# EVT
## compute whether a pixel's most voted analog is forest or non-forest
EVT_predicted_df <- EVT_predicted_df %>%
    mutate(Forest_NonForest_predicted = map_chr(Simplified_EVT_Name, forest_nonforest))

EVT_forest_match <- EVT_predicted_df %>%
    mutate(Forest_match = as.numeric(Forest_NonForest_actual == Forest_NonForest_predicted)) %>%
    group_by(dataset) %>%
    summarize(match = mean(Forest_match), .groups = "drop")

# ESP
## compute whether a pixel's most voted analog is forest or non-forest
ESP_predicted_df <- ESP_predicted_df %>%
    mutate(Forest_NonForest_predicted = map_chr(Simplified_ESP_Name, forest_nonforest))

ESP_forest_match <- ESP_predicted_df %>%
    mutate(Forest_match = as.numeric(Forest_NonForest_actual == Forest_NonForest_predicted)) %>%
    group_by(dataset) %>%
    summarize(match = mean(Forest_match), .groups = "drop")

# write out the results for veg matching and forest matching after combining datasets and BPS EVT ESP
combined_prediction_strength <- bind_rows(
    BPS_prediction_strength %>% mutate(source = "BPS"),
    EVT_prediction_strength %>% mutate(source = "EVT"),
    ESP_prediction_strength %>% mutate(source = "ESP")
)
write_csv(combined_prediction_strength, "data/washington/combined_prediction_strength.csv")
