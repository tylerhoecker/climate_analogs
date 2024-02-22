# Veg extraction and vote fn ---------------------------------------------------
#' @description A function that identifies the vegetation associated with a set
#' of analogs and makes a projection based on weighted votes for vegetation classes
#' @param analog_data A set of analogs for a single point in tabular form
#' As returned by mahalanonis_D_fn.R Expects: "f_x","f_y","a_x","a_y","md","dist_m","sigma"   
#' @param impact_data A SpatRast of environmental data for which projections will be made
#' @param n_analogs How many analogs to keep? May be smaller or equal to the total number
#' @param weighted Use a distance-weighted vote or a proportional vote?
#' @param n_projections How many projected classes of AIM to return? (ie, top vote, top 3 vote-getters, etc.)
#' @return A table of analog and project environmental attributes

veg_vote_fn <- function(
  analog_data,
  impact_data,
  n_analog_use,
  weighted,
  n_projections
){
  
  if(weighted == TRUE){
    sort_val <- "wt_vote"
  } else {
    sort_val  <- "raw_prop"
  }

  if(n_analog_use > nrow(analog_data)){
    print(
      paste0(
        "Number of analogs to use (", n_analog_use, ") ",
        "exceeds number available (", nrow(analog_data), ")"
      )
    )
    return(NULL)
  }
  
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
    # Turn distances into weights
    mutate(weight = scales::rescale(dist_m, to = c(1,0))) |>
    # For each unique AIM code...
    group_by(focal_id, aim) |> 
    summarize(# Preserve coordinates
              x = first(f_x), y = first(f_y),
              # For validation, 'true' aim code
              observed = first(observed),
              # Weighted vote: sum of weights for BPS code
              wt_vote = sum(weight),
              # Raw vote: count of analogs w/ this BPS code
              raw_vote = n(),
              raw_prop = raw_vote/n_analog_use,
              # Statistics for each BPS code
              mean_sigma = mean(sigma),
              min_dist = min(dist_m),
              mean_dist = mean(dist_m),
              .groups = "drop_last") |> 
    # This sequence saves the rows for each focal point, one for each BPS class,
    # that account for 90% of the (unweighted) analogs 
    arrange(desc(raw_prop), .by_group = T) |> 
    filter(
      cumsum(
        raw_prop == accumulate(
          raw_prop, ~ ifelse(.x <= .90, .x + .y, .y))) == 1) |> 
    mutate(n_aim_90 = length(aim)) |>
    # Sort by the weighted or unweighted votes for AIM classes
    arrange(desc(sort_val), .by_group = T) |> 
    ungroup() |>
    # Select top n projections
    slice_head(n = n_projections) 

  return(result_df)
}




  