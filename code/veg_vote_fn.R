# ------------------------------------------------------------------------------------------------
# Project vegetation for validation set ----------------------------------------------------------
# Read in the sample points - a data.table of focal locations and 1000 best analogs > 50 km away 
val_analogs <- readRDS(here("data/validation/random_tiles.Rds")) 

# Test effect of using top 100 analogs instead of 1000
val_analogs_100  <- val_analogs |>
  group_by(focal_id) |>
  arrange(sigma) |>
  slice_head(n = 100) |>
  ungroup()

# Read in Landfire BPS data ("impact data")
lf_bps <- rast(here('data/landfire/lf_bps_west_220.tif'))

val_vegproj <- 
  veg_vote_fn(analog_data = val_analogs, 
              # If raster, must provide path, not a SpatRast
              impact_path = here('data/landfire/lf_bps_west_220.tif'),
              n_analog_keep = 1000)

val_stats <- val_vegproj |>
  filter(!is.na(bps_true), !is.na(bps)) |>
  group_by(focal_id) |>
    # Sort by the weighted votes for BPS classes
  arrange(desc(wt_vote), .by_group = T) |>
  slice_head(n = 1) |>
  mutate(correct = ifelse(bps == bps_true, 1, 0)) |>
  group_by(bps_true) |>
  mutate(occurrence = n()) |>
  ungroup()

# Correct classifications = 0.404 with 1000 analogs; 0.38 with 100 analogs; ~0.3 higher with wt_vote
sum(val_stats$correct)/length(val_stats$correct)
# Random classification rate 0.014
1/length(unique(val_stats$bps_true))

top_bps <- val_stats |>
  group_by(bps_true) |>
  summarise(occurrence = first(occurrence)) |>
  arrange(desc(occurrence)) |>
  #slice_head(n = 25, ) |>
  select(bps_true, occurrence) 

val_stats_bps <- val_vegproj |>
  full_join(x = _, y = top_bps) |>
  # filter(bps_true %in% top_bps$bps_true) |>
  group_by(focal_id) |>
  # Sort by the weighted votes for BPS classes
  arrange(desc(wt_vote), .by_group = T) |>
  slice_head(n = 1) |>
  group_by(bps_true) |>
  mutate(correct = ifelse(bps == bps_true, 1, 0)) |>
  summarise(mean_sigma = mean(mean_sigma),
            mean_vote = mean(wt_vote),
            correct_rate = mean(correct, na.rm = T),
            occurrence = first(occurrence)) |>
  mutate(bps_true = as.factor(bps_true),
         bps_true = forcats::fct_reorder(bps_true, correct_rate)) |>
  filter(!is.na(correct_rate)) 

ggplot(val_stats_bps) +
  geom_col(aes(x = bps_true, fill = mean_vote, y = correct_rate)) +
  coord_flip()

cor(val_stats_bps$occurrence, val_stats_bps$correct_rate)

ggplot(val_stats_bps, aes(x = occurrence, y = correct_rate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  coord_cartesian(ylim = c(0,1))
