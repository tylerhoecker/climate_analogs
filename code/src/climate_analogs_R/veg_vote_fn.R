library(data.table)
library(caret)
library(tidyverse)
library(caret)
library(terra)

simplify_lf_names <- function(name) {
  # Convert to lowercase for uniform processing
  simplified_name <- str_to_lower(name)


  # Identify dominant species within the name
  # This example assumes common species names are of interest. Extend the list as needed.
  physiogamy_list <- c(
    "lodgepole pine", "ponderosa pine", "douglas-fir", "englemann spruce",
    "subalpine fir", "quaking aspen", "aspen", "whitebark pine", "blue oak",
    "valley oak", "black oak", "sugar maple", "red maple", "bigleaf maple",
    "green ash", "black ash", "white ash", "cottonwood", "tamarack",
    "paper birch", "yellow birch", "water birch", "creosotebush", "blackbrush",
    "white bursage", "mormon-tea", "sparsely vegetated", "bigtooth maple",
    "sagebrush", "saltbush", "buffalograss", "bluestem", "needlegrass",
    "pinyon pine", "juniper", "willow", "alder", "sedge", "hemlock",
    "cypress", "redwood", "western redcedar", "pinyon-juniper", "western juniper", "alder",
    "herbaceous", "dryas", "shrubland", "grassland", "forest", "wetland", "riparian", "water",
    "ice", "snow", "sand", "woodland", "savanna", "steppe", "chaparral", "scrub", "heath",
    "moor", "mire", "fen", "bog", "swamp", "marsh", "wet meadow", "dry meadow", "tundra",
    "alpine", "subalpine", "montane", "foothill", "lowland", "upland", "highland",
    "lowland", "floodplain", "alluvial", "riparian", "mixed", "conifer", "deciduous",
    "evergreen", "broadleaf", "needleleaf", "hardwood", "softwood", "coniferous",
    "longleaf pine", "post oak", "forest and woodland", "subland-basin",
    "shrubland-semi-desert", "shrubland-upland", "mountain mahogany", "limber-bristlecone pine",
    "pine-oak", "conifer-hardwood", "hardwoods", "forest-hemlock", "pine(-oak)",
    "pine-hemlock-hardwood", "dune", "encinal", "sitka spruce", "beech-maple",
    "tallgrass", "mixedgrass", "praire", "pine", "oak", "maple", "spruce", "fir",
    "cedar", "cypress", "hemlock", "juniper", "dwarf-shrub", "black spuce", "dwarf-tree",
    "peatland", "desert scrub", "desert riparian", "meadow", "praire", "barrens", "savanna"
  )

  # Find the dominant species in the name
  found_species <- str_extract_all(simplified_name, paste(physiogamy_list, collapse = "|")) %>%
    unlist() %>%
    unique() %>%
    str_to_title()

  if (length(found_species) > 0) {
    return(paste(found_species, collapse = "-"))
  } else {
    # If no dominant species found, simplify to ecosystem type
    ecosystem_type <- case_when(
      str_detect(simplified_name, "forest") ~ "Forest",
      str_detect(simplified_name, "shrub") ~ "Shrubland",
      str_detect(simplified_name, "grass") ~ "Grassland",
      str_detect(simplified_name, "water") ~ "Open Water",
      str_detect(simplified_name, "ice|snow") ~ "Ice/Snow",
      str_detect(simplified_name, "rock") ~ "Rock",
      str_detect(simplified_name, "sand") ~ "Sand",
      str_detect(simplified_name, "riparian") ~ "Riparian",
      str_detect(simplified_name, "wetland") ~ "Wetland",
      TRUE ~ "Other"
    )
    return(ecosystem_type)
  }
}
forest_nonforest <- function(group_veg) {
  group_veg <- str_to_lower(group_veg)
  case_when(
    str_detect(group_veg, "forest|woodland|savanna") ~ "Forest",
    TRUE ~ "Non-Forest"
  )
}
join_fn_BPS <- function(input_dt, input_csv) {
  # Perform the join using the provided by_statement
  # result <- input_dt %>%
  #     select(dataset, f_x, f_y, {{ input_column }}) %>%
  #     left_join(input_csv, by = by_statement) %>%
  #     select({{ input_column }}, {{ count_column }})
  selection <- input_dt[, .(f_x, f_y, BPS_CODE)]
  left_join <- input_csv[selection, on = "BPS_CODE"][
    , .(BPS_CODE, simplified_name, Forest_NonForest)
  ]

  return(left_join)
}
join_fn_wclasses <- function(input_dt, input_column, input_csv, by_statement, count_column) {
  # Perform the join using the provided by_statement
  # result <- input_dt %>%
  #     select(dataset, f_x, f_y, {{ input_column }}) %>%
  #     left_join(input_csv, by = by_statement) %>%
  #     select({{ input_column }}, {{ count_column }})
  selection <- input_dt[, .(dataset, f_x, f_y, input_column), env = list(input_column = input_column)]
  left_join <- input_csv[selection, on = by_statement, env = list(by_statement = by_statement)][
    , .(input_column, count_column),
    env = list(input_column = input_column, count_column = count_column)
  ]

  return(left_join)
}


calculate_top_vote <- function(input_dt, input_column, vote_column) {
  top_vote <- input_dt[, .(f_x, f_y, input_column), env = list(input_column = input_column)][
    , .(count = .N),
    by = .(f_x, f_y, input_column),
    env = list(input_column = input_column)
  ][, .(f_x, f_y, input_column, count),
    env = list(input_column = input_column)
  ]
  top_vote <- top_vote[top_vote[, .I[which.max(count)],
    by = .(f_x, f_y, input_column),
    env = list(input_column = input_column)
  ]$V1][, vote_column := count, env = list(vote_column = vote_column)][
    , count := NULL
  ]
  # top_vote <- input_dt %>%
  #     select(f_x, f_y, {{ input_column }},  %>%
  #     group_by(f_x, f_y,  input_column }}) %>%
  #     summarize({{ vote_column }} := n()) %>%
  #     slice_max(order_by = {{ vote_column }}, ..., with_ties = FALSE) %>%
  #     ungroup()
  focal_points <- unique(input_dt, by = c("f_x", "f_y"), cols = c("BPS_name_actual", "Forest_NonForest_actual"))
  combined_dt <- focal_points[top_vote, on = .(f_x, f_y)]
  combined_dt <- unique(combined_dt, by = c("f_x", "f_y"))

  return(combined_dt)
}
calculate_top_vote_wsigma <- function(input_dt, input_column, vote_column) {
  top_vote <- input_dt[, .(f_x, f_y, input_column, sigma), env = list(input_column = input_column)][
    # rescale sigma
    , sigma_rescaled := 1 - (sigma - min(sigma)) / (max(sigma) - min(sigma))
  ][
    , .(count = .N, sigma_score = sum(sigma_rescaled)),
    by = .(f_x, f_y, input_column),
    env = list(input_column = input_column)
  ][, .(f_x, f_y, input_column, count, sigma_score),
    env = list(input_column = input_column)
  ]
  top_vote <- top_vote[top_vote[, .I[which.max(count) & which.max(sigma_score)],
    by = .(f_x, f_y, input_column),
    env = list(input_column = input_column)
  ]$V1][, vote_column := count, env = list(vote_column = vote_column)][
    , count := NULL
  ]
  # top_vote <- input_dt %>%
  #     select(f_x, f_y, {{ input_column }},  %>%
  #     group_by(f_x, f_y,  {{ input_column }}) %>%
  #     summarize({{ vote_column }} := n()) %>%
  #     slice_max(order_by = {{ vote_column }}, ..., with_ties = FALSE) %>%
  #     ungroup()
  focal_points <- unique(input_dt, by = c("f_x", "f_y"), cols = c("BPS_name_actual", "Forest_NonForest_actual"))
  combined_dt <- focal_points[top_vote, on = .(f_x, f_y)]
  combined_dt <- combined_dt[combined_dt[, .I[which.max(sigma_score)], by = .(f_x, f_y)]$V1]

  return(combined_dt)
}
vote_filter_fn <- function(input_dt, dataset_in, count_column, vote_column) {
  input_dt[dataset == dataset_in, .(f_x, f_y, count_column, vote_column), env = list(count_column = count_column, vote_column = vote_column)]
  # input_dt %>%
  #     filter(dataset == dataset) %>%
  #     select(-dataset) %>%
  #     ungroup()
}
convert2Values <- function(input_dt, input_column) {
  factors <- input_dt[, .(input_column), env = list(input_column = input_column)]
  factors <- unique(factors) %>%
    pull() %>%
    as.factor()
  numerics <- as.numeric(factors)
  characters <- as.character(factors)
  from_to_matrix <- data.frame(ID = numerics, category = characters) %>%
    drop_na()
}

rasterize_predicted_veg <- function(input_dt, template, field) {
  cats_levels <- convert2Values(input_dt, field)
  input_dt[complete.cases(field), env = list(field = field)][
    , field := as.numeric(as.factor(field)),
    env = list(field = field)
  ] %>%
    vect(geom = c("f_x", "f_y"), crs = "EPSG:4326") %>%
    rasterize(template, field = field, fun = "first") %>%
    categories(value = cats_levels)
}
# ------------------------------------------------------------------------------------------------
# Project vegetation for validation set ----------------------------------------------------------
# Read in the sample points - a data.table of focal locations and 1000 best analogs > 50 km away
# val_analogs <- readRDS("data/validation/random_tiles.Rds")

# # Test effect of using top 100 analogs instead of 1000
# val_analogs_100 <- val_analogs |>
#   group_by(focal_id) |>
#   arrange(sigma) |>
#   slice_head(n = 100) |>
#   ungroup()

# # Read in Landtire BPS data ("impact data")
# lf_bps <- rast("data/landtire/lf_bps_west_220.tif")

# val_vegproj <-
#   veg_vote_fn(
#     analog_data = val_analogs,
#     # If raster, must provide path, not a SpatRast
#     impact_path = "data/landtire/lf_bps_west_220.tif",
#     n_analog_keep = 1000
#   )

# val_stats <- val_vegproj |>
#   filter(!is.na(bps_true), !is.na(bps)) |>
#   group_by(focal_id) |>
#   # Sort by the weighted votes for BPS classes
#   arrange(desc(wt_vote), .by_group = T) |>
#   slice_head(n = 1) |>
#   mutate(correct = ifelse(bps == bps_true, 1, 0)) |>
#   group_by(bps_true) |>
#   mutate(occurrence = n()) |>
#   ungroup()

# # Correct classifications = 0.404 with 1000 analogs; 0.38 with 100 analogs; ~0.3 higher with wt_vote
# sum(val_stats$correct) / length(val_stats$correct)
# # Random classification rate 0.014
# 1 / length(unique(val_stats$bps_true))

# top_bps <- val_stats |>
#   group_by(bps_true) |>
#   summarise(occurrence = first(occurrence)) |>
#   arrange(desc(occurrence)) |>
#   # slice_head(n = 25, ) |>
#   select(bps_true, occurrence)

# val_stats_bps <- val_vegproj |>
#   full_join(x = _, y = top_bps) |>
#   # filter(bps_true %in% top_bps$bps_true) |>
#   group_by(focal_id) |>
#   # Sort by the weighted votes for BPS classes
#   arrange(desc(wt_vote), .by_group = T) |>
#   slice_head(n = 1) |>
#   group_by(bps_true) |>
#   mutate(correct = ifelse(bps == bps_true, 1, 0)) |>
#   summarise(
#     mean_sigma = mean(mean_sigma),
#     mean_vote = mean(wt_vote),
#     correct_rate = mean(correct, na.rm = T),
#     occurrence = first(occurrence)
#   ) |>
#   mutate(
#     bps_true = as.factor(bps_true),
#     bps_true = forcats::fct_reorder(bps_true, correct_rate)
#   ) |>
#   filter(!is.na(correct_rate))

# ggplot(val_stats_bps) +
#   geom_col(aes(x = bps_true, fill = mean_vote, y = correct_rate)) +
#   coord_flip()

# cor(val_stats_bps$occurrence, val_stats_bps$correct_rate)

# ggplot(val_stats_bps, aes(x = occurrence, y = correct_rate)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw() +
#   coord_cartesian(ylim = c(0, 1))

# kappa validation
# calculate cohens kappa
# build confusion matrix
build_accuracy_stats <- function(input_dt, actual_column, predicted_column) {
  levels <- unique(c(input_dt[[actual_column]], input_dt[[predicted_column]]))
  input_dt[[actual_column]] <- factor(input_dt[[actual_column]], levels = levels)
  input_dt[[predicted_column]] <- factor(input_dt[[predicted_column]], levels = levels)
  confusion_mat <- confusionMatrix(input_dt[[predicted_column]], input_dt[[actual_column]], mode = "prec_recall")
  accuracy <- confusion_mat$overall["Accuracy"]
  kappa <- confusion_mat$overall["Kappa"]
  output <- data.frame(accuracy = accuracy, kappa = kappa)
  return(output)
}

build_accuracy_stats_wclass <- function(input_dt, actual_column, predicted_column, class_col) {
  classes <- unique(input_dt[[class_col]])
  length_classes <- length(classes)
  output_list <- data.frame(
    class = length_classes, accuracy = length_classes,
    kappa = length_classes
  )
  for (class_i in seq_along(classes)) {
    class <- classes[[class_i]]
    filtered_dt <- input_dt[class_col == class, env = list(class_col = class_col)]
    # convert to factor
    filtered_dt[[actual_column]] <- as.factor(filtered_dt[[actual_column]])
    filtered_dt[[predicted_column]] <- as.factor(filtered_dt[[predicted_column]])
    # ensure the levels are the same
    ## create levels that cover each possible value
    levels <- unique(c(filtered_dt[[actual_column]], filtered_dt[[predicted_column]]))
    # ensure the levels are the same
    filtered_dt[[actual_column]] <- factor(filtered_dt[[actual_column]], levels = levels)
    filtered_dt[[predicted_column]] <- factor(filtered_dt[[predicted_column]], levels = levels)
    confusion_mat <- confusionMatrix(filtered_dt[[predicted_column]], filtered_dt[[actual_column]], mode = "prec_recall")
    accuracy <- confusion_mat$overall["Accuracy"]
    kappa <- confusion_mat$overall["Kappa"]
    output_list[class_i, ] <- c(class, accuracy, kappa)
  }
  return(output_list)
}

filter_by_sigma <- function(input_dt, sigma_threshold) {
  input_dt[sigma < sigma_threshold, env = list(sigma_threshold = sigma_threshold)]
}

setup_veg_prediction_bps <- function(input_csv, input_outline_sv, template_r = NULL, BPS = NULL) {
  # filter the input data by sigma
  input_csv <- filter_by_sigma(input_csv, 2)

  # read in a tempalte raster for climate
  if (is.null(template_r)) {
    template_local <- rast("data/climate/topoterra_hist_1961-1990.tif", 1) %>%
      crop(., project(buffer(input_outline_sv, 500 * 1000), crs(.)))
  } else {
    template_local <- template_r %>%
      crop(., project(buffer(input_outline_sv, 500 * 1000), crs(.)))
  }
  # replace all valid values with 1, else NA
  rcl <- matrix(c(0, Inf, 1, -Inf, 0, NA), byrow = TRUE, ncol = 3)
  template_local <- classify(template_local, rcl, others = NA)

  # grab all analogs
  analogs_dt <- unique(input_csv[, .(a_x, a_y)])
  analogs_v <- vect(analogs_dt, geom = c("a_x", "a_y"), crs = "EPSG:4326")

  if (is.null(BPS)) {
    BPS_local <- rast("data/veg_data/LC20_BPS_220.tif") %>%
      crop(., project(buffer(input_outline_sv, 500 * 1000), crs(.))) %>%
      project(., crs(template_local), method = "near") %>%
      resample(., template_local, method = "near")
  } else {
    # check if BPS matches template
    if (ext(BPS) == ext(template_r) & crs(BPS) == crs(template_r)) {
      BPS_local <- BPS %>%
        crop(., buffer(input_outline_sv, 500 * 1000))
    } else {
      BPS_local <- BPS %>%
        crop(., project(buffer(input_outline_sv, 500 * 1000), crs(.))) %>%
        project(., crs(template_local), method = "near") %>%
        resample(., template_local, method = "near")
    }
  }
  activeCat(BPS_local) <- "BPS_CODE"
  # extract BPS and EVT to analogs
  BPS_analogs <- terra::extract(BPS_local, analogs_v,
    na.rm = TRUE,
    xy = FALSE
  )

  # bind  the extracted values to the analogs dataframe
  analogs_dt_full <- data.table::data.table(analogs_dt,
    BPS_CODE = as.numeric(as.character(BPS_analogs$BPS_CODE))
  )

  input_csv_full <- input_csv[analogs_dt_full, on = .(a_x = a_x, a_y = a_y)][BPS_CODE > 31]
  rm(input_csv)
  gc()

  # adjust vegetation names to simplify the datasets --------
  mapping_categorical <- fread("data/veg_data/mapping_categorical_v5.csv")[, .(BPS_CODE, Cluster_name)]
  mapping_categorical <- unique(mapping_categorical, by = "BPS_CODE")
  # rename cluster name to simplified name
  mapping_categorical <- mapping_categorical[, simplified_name := Cluster_name][, Cluster_name := NULL]
  # select(BPS_CODE, final_id, Cluster_name)
  BPS_csv <- fread("data/veg_data/LF20_BPS_220.csv")[, ":="(
    Forest_NonForest = map_chr(BPS_NAME, forest_nonforest))][
    , .(BPS_CODE, Forest_NonForest)
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
  best_analogs_BPS <- join_fn_BPS(input_csv_full, BPS_csv)
  # join the plurality votes to the full dataframe
  input_csv_w_analogs <- cbind(input_csv_full, best_analogs_BPS)
  input_csv_w_analogs <- input_csv_w_analogs[, ":="(BPS_name_predicted = simplified_name, simplified_name = NULL,
    Forest_NonForest_predicted = Forest_NonForest, Forest_NonForest = NULL,
    BPS_CODE_predicted = BPS_CODE, BPS_CODE = NULL)]
  gc()
  # wa_dt_full_plurality <- wa_dt_full %>%
  #     cbind(c(best_analogs_BPS, best_analogs_EVT, best_analogs_ESP)) %>%
  #     select(-BPS_CODE, -EVT_NAME, -ESP_NAME)

  # Compute actual BPS
  ## take all focal points and extract the BPS code
  focal_points <- unique(input_csv_full[, .(f_x, f_y)]) %>%
    vect(geom = c("f_x", "f_y"), crs = "EPSG:4326")
  # focal_points <- input_csv_full %>%
  #     select(f_x, f_y) %>%
  #     distinct(f_x, f_y) %>%
  #     vect(geom = c("f_x", "f_y"), crs = "EPSG:4326")

  BPS_actual_extracted <- terra::extract(BPS_local, focal_points, xy = FALSE, ID = FALSE) %>%
    pull()

  # join to the CSV
  BPS_actual_dt <- focal_points %>%
    as.data.frame(geom = "XY") %>%
    as.data.table()
  BPS_actual_dt <- BPS_actual_dt[, BPS_actual := as.numeric(as.character(BPS_actual_extracted))]
  BPS_actual_dt <- BPS_csv[BPS_actual_dt, on = .(BPS_CODE = BPS_actual)][BPS_CODE > 31]
  # rename
  BPS_actual_dt <- BPS_actual_dt[, .(BPS_name_actual = simplified_name, Forest_NonForest_actual = Forest_NonForest, f_x = x, f_y = y, BPS_CODE_actual = BPS_CODE)]

  # BPS_actual_dt <- focal_points %>%
  #     as.data.frame(geom = "XY") %>%
  #     mutate(BPS_actual = as.numeric(as.character(BPS_actual_extracted))) %>%
  #     left_join(BPS_csv, by = c("BPS_actual" = "BPS_CODE")) %>%
  #     rename(
  #         BPS_name_actual = simplified_name,
  #         Forest_NonForest_actual = Forest_NonForest,
  #         f_x = x,
  #         f_y = y
  #     )

  # join to the analog dataframe
  BPS_predicted_dt <- input_csv_w_analogs[, .(f_x, f_y, a_x, a_y, sigma, dist_km, BPS_CODE_predicted, BPS_name_predicted, Forest_NonForest_predicted)]
  rm(input_csv_w_analogs)
  gc()
  # rename to BPS_predicted
  BPS_predicted_dt <- BPS_actual_dt[BPS_predicted_dt, on = .(f_x, f_y)]
  rm(BPS_actual_dt)
  gc()
  output_dt <- BPS_predicted_dt[complete.cases(BPS_predicted_dt), ]

  gc()
  return(output_dt)
}

compute_accuracy_plurality_bps <- function(input_dt) {
  combined_dt <- calculate_top_vote(input_dt, "BPS_name_predicted", "BPS_vote")
  reserved <- combined_dt[, .(BPS_name_actual, BPS_name_predicted, BPS_vote)]
  accuracy_stats <- build_accuracy_stats(combined_dt, "BPS_name_actual", "BPS_name_predicted")
  rm(combined_dt)
  gc()
  combined_dt <- calculate_top_vote(input_dt, "Forest_NonForest_predicted", "Forest_NonForest_vote")
  accuracy_forest_plurality <- build_accuracy_stats(combined_dt, "Forest_NonForest_actual", "Forest_NonForest_predicted")
  combined_dt <- cbind(combined_dt, reserved)[, forest_vote := BPS_vote][, BPS_vote := NULL]

  return(list(
    method = "plurality",
    summarized_dt = combined_dt,
    accuracy_simplified_bps = accuracy_stats,
    accuracy_forest = accuracy_forest_plurality
  ))
}

compute_accuracy_plurality_bps_wsigma <- function(input_dt) {
  combined_dt <- calculate_top_vote_wsigma(input_dt, "BPS_name_predicted", "BPS_vote")
  reserved <- combined_dt[, .(BPS_name_actual, BPS_name_predicted, BPS_vote)]
  accuracy_stats <- build_accuracy_stats(combined_dt, "BPS_name_actual", "BPS_name_predicted")
  rm(combined_dt)
  gc()
  combined_dt <- calculate_top_vote_wsigma(input_dt, "Forest_NonForest_predicted", "Forest_NonForest_vote")
  accuracy_forest_plurality <- build_accuracy_stats(combined_dt, "Forest_NonForest_actual", "Forest_NonForest_predicted")
  combined_dt <- cbind(combined_dt, reserved)[, forest_vote := BPS_vote][, BPS_vote := NULL]
  return(list(
    method = "plurality with sigma",
    summarized_dt = combined_dt,
    accuracy_simplified_bps = accuracy_stats,
    accuracy_forest = accuracy_forest_plurality
  ))
}

calculate_min_ranks <- function(input_csv_full, vegarg1, vegarg2) {
  # rank within focal points
  input_csv_full <- input_csv_full[, rank := frank(sigma), by = .(f_x, f_y)]
  top_ranks <- input_csv_full[, .(rank = sum(rank)), by = .(f_x, f_y, vegarg1, vegarg2), env = list(vegarg1 = vegarg1, vegarg2 = vegarg2)]
  top_ranks <- top_ranks[top_ranks[, .I[which.min(rank)], by = .(f_x, f_y)]$V1]
  return(top_ranks)
}

calculate_top_sigma <- function(input_csv_full, vegarg1, vegarg2) {
  rescaled_input_csv_full <- input_csv_full[, sigma_rescaled := 1 - (sigma - min(sigma)) / (max(sigma) - min(sigma))]
  top_sigma <- rescaled_input_csv_full[, .(sigma_score = sum(sigma_rescaled)), by = .(f_x, f_y, vegarg1, vegarg2), env = list(vegarg1 = vegarg1, vegarg2 = vegarg2)]
  top_sigma <- top_sigma[top_sigma[, .I[which.max(sigma_score)], by = .(f_x, f_y)]$V1]
  return(top_sigma)
}

compute_accuracy_sigma_bps <- function(input_csv_full) {
  BPS_analog_scores <- calculate_top_sigma(input_csv_full, "BPS_name_predicted", "BPS_name_actual")
  reserved <- BPS_analog_scores[, .(BPS_name_actual, BPS_name_predicted, sigma_score)][
    , sigma_score_BPS := sigma_score
  ][
    , sigma_score := NULL
  ]
  accuracy_stats <- build_accuracy_stats(BPS_analog_scores, "BPS_name_actual", "BPS_name_predicted")
  rm(BPS_analog_scores)
  gc()
  BPS_analog_scores <- calculate_top_sigma(input_csv_full, "Forest_NonForest_predicted", "Forest_NonForest_actual")
  accuracy_forest_sigma <- build_accuracy_stats(BPS_analog_scores, "Forest_NonForest_actual", "Forest_NonForest_predicted")
  BPS_analog_scores <- cbind(BPS_analog_scores, reserved)[, sigma_score_forest := sigma_score][, sigma_score := NULL]
  return(list(
    method = "sigma",
    simplified_dt = BPS_analog_scores,
    accuracy_simplified_bps = accuracy_stats,
    accuracy_forest = accuracy_forest_sigma
  ))
}

calculate_top_distance <- function(input_csv_full, vegarg1, vegarg2) {
  rescaled_input_csv_full <- input_csv_full[, dist_rescaled := 1 - (dist_km - min(dist_km)) / (max(dist_km) - min(dist_km))]
  top_distance <- rescaled_input_csv_full[, .(dist = sum(dist_rescaled)), by = .(f_x, f_y, vegarg1, vegarg2), env = list(vegarg1 = vegarg1, vegarg2 = vegarg2)]
  top_distance <- top_distance[top_distance[, .I[which.max(dist)], by = .(f_x, f_y)]$V1]
  return(top_distance)
}
compute_accuracy_distance_bps <- function(input_csv_full) {
  BPS_analog_scores <- calculate_top_distance(input_csv_full, "BPS_name_predicted", "BPS_name_actual")
  reserved <- BPS_analog_scores[, .(BPS_name_actual, BPS_name_predicted, dist)][
    , dist_BPS := dist
  ][
    , dist := NULL
  ]
  accuracy_stats <- build_accuracy_stats(BPS_analog_scores, "BPS_name_actual", "BPS_name_predicted")
  rm(BPS_analog_scores)
  gc()
  BPS_analog_scores <- calculate_top_distance(input_csv_full, "Forest_NonForest_predicted", "Forest_NonForest_actual")
  accuracy_forest_dist <- build_accuracy_stats(BPS_analog_scores, "Forest_NonForest_actual", "Forest_NonForest_predicted")
  BPS_analog_scores <- cbind(BPS_analog_scores, reserved)[, dist_forest := dist][, dist := NULL]
  return(list(
    method = "distance",
    simplified_dt = BPS_analog_scores,
    accuracy_simplified_bps = accuracy_stats,
    accuracy_forest = accuracy_forest_dist
  ))
}

calculate_top_combined <- function(input_csv_full, vegarg1, vegarg2) {
  rescaled_input_csv_full <- input_csv_full[, sigma_rescaled := 1 - (sigma - min(sigma)) / (max(sigma) - min(sigma))]
  rescaled_input_csv_full <- rescaled_input_csv_full[, dist_rescaled := 1 - (dist_km - min(dist_km)) / (max(dist_km) - min(dist_km))]
  BPS_analog_scores <- rescaled_input_csv_full[, .(score = sum(sigma_rescaled) + sum(dist_rescaled)), by = .(f_x, f_y, vegarg1, vegarg2), env = list(vegarg1 = vegarg1, vegarg2 = vegarg2)]
  BPS_analog_scores <- BPS_analog_scores[BPS_analog_scores[, .I[which.max(score)], by = .(f_x, f_y)]$V1]
  return(BPS_analog_scores)
}
compute_accuracy_combined_bps <- function(input_csv_full) {
  BPS_analog_scores <- calculate_top_combined(input_csv_full, "BPS_name_predicted", "BPS_name_actual")
  reserved <- BPS_analog_scores[, .(BPS_name_actual, BPS_name_predicted, score)][
    , score_BPS := score
  ][
    , score := NULL
  ]
  accuracy_stats <- build_accuracy_stats(BPS_analog_scores, "BPS_name_actual", "BPS_name_predicted")
  rm(BPS_analog_scores)
  gc()
  BPS_analog_scores <- calculate_top_combined(input_csv_full, "Forest_NonForest_predicted", "Forest_NonForest_actual")
  accuracy_forest_combined <- build_accuracy_stats(BPS_analog_scores, "Forest_NonForest_actual", "Forest_NonForest_predicted")
  BPS_analog_scores <- cbind(BPS_analog_scores, reserved)[, score_forest := score][, score := NULL]
  return(list(
    method = "combined",
    simplified_dt = BPS_analog_scores,
    accuracy_simplified_bps = accuracy_stats,
    accuracy_forest = accuracy_forest_combined
  ))
}

build_accuracy_tables <- function(input_list) {
  final_table <- purrr::map_dfr(input_list, \(dataset){
    method <- dataset$method
    bps_accuracy <- dataset$accuracy_simplified_bps
    forest_accuracy <- dataset$accuracy_forest
    full_table <- rbind(
      data.frame(Prediction = "BPS", bps_accuracy),
      data.frame(Prediction = "Forest", forest_accuracy)
    ) %>%
      mutate(Method = method)
  })
}

# create reclassify matrix and color palette for a single layer input_spatRast based on mapping file
create_reclass_cpal <- function(input_spatRast) {
  if (!is(input_spatRast, "SpatRaster")) {
    stop("input_spatRast must be a SpatRaster object")
  }
  if (nlyr(input_spatRast) > 1) {
    stop("input_spatRast must have only one layer")
  }
  # read in the mapping file
  mapping <- read_csv("data/veg_data/mapping_palette.csv", show_col_types = FALSE)
  # create a reclassify matrix
  ## grab levels from linput_spatRast
  categorical_df <- cats(input_spatRast) %>%
    as.data.frame()
  if (nrow(categorical_df) != 0) {
    factor_df <- setNames(categorical_df, c("ID", "name"))
    joined_df <- factor_df %>%
      left_join(mapping, by = c("name" = "name"))
  } else {
    factor_df <- data.frame(ID = as.vector(unique(values(input_spatRast)))) %>%
      drop_na()
    joined_df <- factor_df %>%
      left_join(mapping, by = c("ID" = "value")) %>%
      mutate(value = ID)
  }
  ## join mapping to levels based on name

  # create reclassify matrix
  rcl <- matrix(c(joined_df$ID, joined_df$value), ncol = 2, byrow = FALSE)

  # create color palette
  cpal <- data.frame(
    value = joined_df$value,
    color = joined_df$color
  )
  labels <- data.frame(
    ID = joined_df$value,
    category = joined_df$name
  )

  return(list(rcl = rcl, cpal = cpal, labels = labels))
}
# function that reclassifies and sets color palette for a single layer input_spatRast based on mapping file
reclassify_cpal <- function(input_spatRast) {
  # create reclassify matrix and color palette
  reclass_cpal <- create_reclass_cpal(input_spatRast)

  # reclassify
  reclassified_rast <- classify(input_spatRast, reclass_cpal$rcl, others = NA)

  # set labels
  reclassified_rast <- categories(reclassified_rast, value = reclass_cpal$labels)
  # set color palette
  coltab(reclassified_rast) <- reclass_cpal$cpal

  return(reclassified_rast)
}
