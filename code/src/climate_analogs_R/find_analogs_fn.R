#' @description Wrapper function for a given dataframe of focal locations and climate data
#' @param focal_data_cov Repeated observations from focal locations (ie, annualized futures) to calculate covariance matrix
#' @param focal_data_mean Mean of repeated observations from focal locations (not required, usually calculated from focal_data_cov)
#' @param analog_data Mean observations from potential analogs (ie historical means from across West)
#' @param var_names Climate variable names; expected to be names for columns in data.tables
#' @param n_analog_pool Size of global analog pool, sampled randomly from extent of analog_data 
#' @param n_analog_use Number of good analogs to keep
#' @param min_dist Minimum distance to analogs in meters (for contemporary validation, exlcude nearby points)
#' @param max_dist Maximum distance to search for analogs in meters
#' @param output_dir Output directory to save resulting data.table as Rds
#' @param use_futures Do parallel processing with futures/furrr?
#' @param n_futures Number of workers to futures
#' #' equal to the number of variables used in the MD calculation
#' @return A vector of "sigmas" - unsquared standard deviations from a chi-squared distribution
source("code/calc_mahalanobis_fn.R")
source("code/calc_sigma_fn.R")
library(future)
library(furrr)
library(data.table)
library(dplyr)
library(purrr)
library(progressr)

max_distance_coordinates <- function(lat, lon, max_dist) {
  # Earth's radius in meters
  R <- 6378.137
  
  # Calculate the change in latitude
  delta_lat <- max_dist / R
  new_lat_north <- lat + (delta_lat * (180 / pi))
  new_lat_south <- lat - (delta_lat * (180 / pi))
  
  # Calculate the change in longitude
  delta_lon <- max_dist / (R * cos(pi * lat / 180))
  new_lon_east <- lon + (delta_lon * (180 / pi))
  new_lon_west <- lon - (delta_lon * (180 / pi))
  
  # Return the new coordinates as a list
  return(list(north = c(new_lat_north, lon), 
              south = c(new_lat_south, lon), 
              east = c(lat, new_lon_east), 
              west = c(lat, new_lon_west)))
}



find_analogs <- function(
    focal_data_cov,
    focal_data_mean,
    analog_data,
    var_names,
    n_analog_pool,
    n_analog_use,
    min_dist,
    max_dist,
    output_dir,
    use_futures = FALSE,
    n_futures
){
    # Map function over all points in supplied dataset

    if (use_futures == TRUE){
        # Set up futures
        options(future.globals.maxSize = object.size(focal_data_cov) + object.size(analog_data) + object.size(var_names) + object.size(focal_data_mean) + object.size(min_dist) + object.size(n_analog_pool) + object.size(n_analog_use) + object.size(output_dir) + 1000000000)
        plan(multisession, workers = n_futures)
        # Run over all points
        with_progress({
        p <- progressor(
          steps = nrow(focal_data_cov[[1]])
          )

        sigma_dt <- seq_len(nrow(focal_data_cov[[1]])) |>
            future_map_dfr(
                \(x){
                     coords <- max_distance_coordinates(
                      lat = focal_data_cov[[1]][x][["y"]],
                      lon = focal_data_cov[[1]][x][["x"]],
                      max_dist = max_dist
                  )
                  
                  ##use the coordinates to filter the analog data using data.table
                  filtered_analog_data <- analog_data[
                      y >= coords$south[1] & y <= coords$north[1] & 
                      x >= coords$west[2] & x <= coords$east[2]
                      ]
               result <- calc_mahalanobis(
                    pt_i = x, 
                    focal_data_cov,
                    focal_data_mean,
                    filtered_analog_data,
                    var_names,
                    n_analog_pool,
                    n_analog_use,
                    min_dist
                ) 
                rm(coords, filtered_analog_data)
                gc()   
                p()
                return(result)
            },  
            .id = "focal_id",
            .options = furrr_options(seed = 3768)
            )
        })
        plan(sequential)
    } else {
        with_progress({
            p <- progressor(
                steps = nrow(focal_data_cov[[1]]),
            )
        sigma_dt <- seq_len(nrow(focal_data_cov[[1]])) |>
            map_dfr(
                \(x){
                     coords <- max_distance_coordinates(
                      lat = focal_data_cov[[1]][x][["y"]],
                      lon = focal_data_cov[[1]][x][["x"]],
                      max_dist = max_dist
                  )
                  
                  ##use the coordinates to filter the analog data using data.table
                  filtered_analog_data <- analog_data[
                      y >= coords$south[1] & y <= coords$north[1] & 
                      x >= coords$west[2] & x <= coords$east[2]
                      ]
                result <- calc_mahalanobis(
                    pt_i = x, 
                    focal_data_cov,
                    focal_data_mean,
                    filtered_analog_data,
                    var_names,
                    n_analog_pool,
                    n_analog_use,
                    min_dist
                )    
                rm(coords, filtered_analog_data)
                gc()
                p()
            },   
            .id = "focal_id"
            )
        })
    }
    
    
    # Save as RDS
    saveRDS(sigma_dt, paste0(output_dir, ".Rds"))
    
    # Explicitly remove and free memory each iteration
    rm(sigma_dt)
    gc()
    return(sigma_dt)
}
