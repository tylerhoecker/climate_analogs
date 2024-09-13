library(data.table)
library(tidyverse)
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


calculate_top_vote <- function(input_df, input_column, vote_column, ...) {
  top_vote <- input_df[, .(f_x, f_y, dataset, input_column), env = list(input_column = input_column)][
    , .(count = .N),
    by = .(f_x, f_y, dataset, input_column),
    env = list(input_column = input_column)
  ][, .(f_x, f_y, dataset, input_column, count),
    env = list(input_column = input_column)
  ]
  top_vote <- top_vote[top_vote[, .I[which.max(count)],
    by = .(f_x, f_y, dataset, input_column),
    env = list(input_column = input_column)
  ]$V1][, vote_column := count, env = list(vote_column = vote_column)][
    , count := NULL
  ]
  # top_vote <- input_df %>%
  #     select(f_x, f_y, {{ input_column }}, dataset) %>%
  #     group_by(f_x, f_y, dataset, {{ input_column }}) %>%
  #     summarize({{ vote_column }} := n()) %>%
  #     slice_max(order_by = {{ vote_column }}, ..., with_ties = FALSE) %>%
  #     ungroup()
  return(top_vote)
}
vote_filter_fn <- function(input_df, dataset_in, count_column, vote_column) {
  input_df[dataset == dataset_in, .(f_x, f_y, count_column, vote_column), env = list(count_column = count_column, vote_column = vote_column)]
  # input_df %>%
  #     filter(dataset == dataset) %>%
  #     select(-dataset) %>%
  #     ungroup()
}
convert2Values <- function(input_df, input_column) {
  factors <- input_df %>%
    select({{ input_column }}) %>%
    distinct() %>%
    pull() %>%
    as.factor()
  numerics <- as.numeric(factors)
  characters <- as.character(factors)
  from_to_matrix <- data.frame(ID = numerics, category = characters) %>%
    drop_na()
}

rasterize_predicted_veg <- function(input_df, template, field) {
  cats_levels <- convert2Values(input_df, {{ field }})
  field_str <- as_string(ensym(field))
  input_df %>%
    drop_na({{ field }}) %>%
    mutate({{ field }} := as.numeric(as.factor({{ field }}))) %>%
    vect(geom = c("f_x", "f_y"), crs = "EPSG:4326") %>%
    rasterize(template, field = field_str, fun = "first") %>%
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

# # Read in Landfire BPS data ("impact data")
# lf_bps <- rast("data/landfire/lf_bps_west_220.tif")

# val_vegproj <-
#   veg_vote_fn(
#     analog_data = val_analogs,
#     # If raster, must provide path, not a SpatRast
#     impact_path = "data/landfire/lf_bps_west_220.tif",
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
