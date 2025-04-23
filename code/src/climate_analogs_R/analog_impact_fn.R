# Veg extraction and vote fn ---------------------------------------------------
#' @description A function that identifies the vegetation associated with a set
#' of analogs and makes a projection based on weighted votes for vegetation classes
#' @param analog_data A set of analogs for a single point in tabular form
#' As returned by mahalanonis_D_fn.R Expects: "f_x","f_y","a_x","a_y","md","dist_km","sigma"   
#' @param impact_data A SpatRast of environmental data for which projections will be made
#' @param n_analogs How many analogs to keep? May be smaller or equal to the total number
#' @param weigh_val Use a distance-weighted vote or a proportional vote?
#' @param n_projections How many projected classes of AIM to return? (ie, top vote, top 3 vote-getters, etc.)
#' @return A table of analog and project environmental attributes
library(dplyr)
library(terra)
library(tidyr)
library(scales)

analog_impact_fn <- function(
  analog_data,
  impact_data,
  n_analog_use,
  weight_val,
  n_projections
){
  
  # FILL IN OTHER WEIGHTING OPTIONS LATER
  
  # Extract impact data
  result_df <- data.frame(
    "aim" = terra::extract(
      impact_data,
      as.matrix(analog_data[,c('a_x','a_y')])
    )[[names(impact_data)]],
    "observed" = terra::extract(
      impact_data,
      as.matrix(analog_data[,c('f_x','f_y')])
    )[[names(impact_data)]])

  result_df <- cbind(result_df, analog_data) |> 
    group_by(focal_id) |>
    # Turn distances into weights
    mutate(weight := scales::rescale(get(!!weight_val), to = c(1,0)),
           n_analogs = n()) |>
    # For each unique AIM code...
    group_by(focal_id, aim) |> 
    reframe(# Preserve coordinates
            x = first(f_x), y = first(f_y),
            # For validation, 'true' aim code
            observed = first(observed),
            # Weighted vote: sum of weights for AIM/BPS code
            wt_vote = sum(weight),
            # Raw vote: count of analogs w/ this AIM/BPS code
            raw_vote = n(),
            raw_prop = raw_vote/first(n_analogs),
            # Statistics for each AIM/BPS code
            mean_sigma = mean(sigma),
            min_sigma = min(sigma),
            min_dist_km = min(dist_km),
            mean_dist_km = mean(dist_km)) |> 
    group_by(focal_id) |>
    # This sequence saves the rows for each focal point, one for each AIM/BPS class,
    # that account for 90% of the (unweighted) analogs 
    arrange(desc(raw_prop), .by_group = T) |>
    filter(
      cumsum(
        raw_prop == accumulate(
          raw_prop, ~ ifelse(.x <= .90, .x + .y, .y))) == 1) |> 
    mutate(n_aim_90 = length(aim)) |>
    # Sort by the weighted or unweighted votes for AIM classes
    arrange(desc(wt_vote), .by_group = T) |> 
    # Select top n projections
    slice_head(n = n_projections) |>
    mutate(rank = row_number()) |>
    pivot_wider(
      id_cols = c(focal_id, x, y, observed),
      names_from = rank,
      values_from = c(aim, wt_vote, raw_vote, raw_prop, mean_sigma, n_aim_90)
    )
  return(result_df)
}