<<<<<<< HEAD
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

calculate_analogs <- function(focal_data_cov,
                              focal_data_mean,
                              analog_data,
                              var_names,
                              n_analog_pool,
                              n_analog_use,
                              min_dist,
                              max_dist,
                              x) {
    if (!is.infinite(max_dist)) {
        coords <- max_distance_coordinates(focal_data_cov[[1]][x, "y"], focal_data_cov[[1]][x, "x"], max_dist)

        # Use the coordinates to filter the analog data using DataFrames
        filtered_analog_data <- analog_data %>%
            dplyr::filter(y >= coords$south[1] &
                y <= coords$north[1] &
                x >= coords$west[2] &
                x <= coords$east[2])
    } else {
        filtered_analog_data <- analog_data
    }

    result_i <- calc_mahalanobis(
        x,
        focal_data_cov,
        focal_data_mean,
        filtered_analog_data,
        var_names,
        n_analog_pool,
        n_analog_use,
        min_dist
    )

    # Remove filtered_analog_data from memory
    filtered_analog_data <- NULL

    # Save the result
    # ...
    return(result_i)
}


run_calculate_analogs <- function(focal_data_cov,
                                  focal_data_mean,
                                  analog_data,
                                  var_names,
                                  n_analog_pool,
                                  n_analog_use,
                                  min_dist,
                                  max_dist,
                                  output_csv,
                                  error_file) {
    result <- future_map_dfr(
        seq_len(size(focal_data_cov[[1]])),
        ~ {
            x <- .
            tryCatch(
                {
                    result_i <- calculate_analogs(
                        focal_data_cov,
                        focal_data_mean,
                        analog_data,
                        var_names,
                        n_analog_pool,
                        n_analog_use,
                        min_dist,
                        max_dist,
                        x
                    )

                    if (x %% 100 == 0 && system("awk '/MemAvailable/ {print $2}' /proc/meminfo",
                        intern = TRUE
                    ) |> as.numeric() / (1024^2) < 10) {
                        open(error_file, "a") %>%
                            write(paste0("Broke at ", x, "\n"))
                    }

                    result_i
                },
                error = function(e) {
                    open(error_file, "a") %>%
                        write(paste0(x, "\n"))
                    data.frame()
                }
            )
        },
        .options = furrr_options(seed = 3768)
    )
}

find_analogs <- function(focal_data_cov,
                         focal_data_mean,
                         analog_data,
                         var_names,
                         n_analog_pool,
                         n_analog_use,
                         min_dist,
                         max_dist,
                         output_file) {
    # Map function over all points in supplied dataset
    # write headers to CSV if csv does not exist
    output_csv <- paste0(output_file, ".csv")
    if (!file.exists(output_csv)) {
        write.csv(
            data.frame(
                a_x = numeric(),
                a_y = numeric(),
                md = numeric(),
                dist_km = numeric(),
                sigma = numeric(),
                f_x = numeric(),
                f_y = numeric()
            ),
            file = output_csv,
            row.names = FALSE
        )
    }

    # Replace error file
    error_file <- paste0("/home/jeff/Github/ClimateAnalogues/climate_analogs/data/", basename(output_file), "_error.txt")
    if (!file.exists(error_file)) {
        file.create(error_file)
    }

    # Check if the data fits in memory
    tile_ids <- spatial_partition(focal_data_cov, focal_data_mean, analog_data, n_analog_use)

    if (!is.null(tile_ids)) {
        # Split the data into tiles
        n_tiles <- max(tile_ids$tile_id)
        for (tile in 1:n_tiles) {
            # Filter the analog pool to only include points within max_dist of the focal points in the tile
            filtered_cov_data <- filter_cov_data(focal_data_cov, tile_extents, tile)
            filtered_data_mean <- filter_mean_data(focal_data_mean, tile_extents, tile)
            filtered_analog_data <- filter_analog_pool(analog_data, tile_extents, tile, max_dist)

            # Run the calculate_analogs function
            result_i <- run_calculate_analogs(
                filtered_cov_data,
                filtered_data_mean,
                filtered_analog_data,
                var_names,
                n_analog_pool,
                n_analog_use,
                min_dist,
                max_dist,
                paste0(output_csv, "_tile", tile),
                paste0(error_file, "_tile", tile)
            )

            write.csv(result_i, file = paste0(output_csv, "_tile", tile), append = TRUE, row.names = FALSE)
            cat("Finished tile", tile, "\n")
        }
        final_result <- "Finished all tiles"
    } else {
        result <- run_calculate_analogs(
            focal_data_cov,
            focal_data_mean,
            analog_data,
            var_names,
            n_analog_pool,
            n_analog_use,
            min_dist,
            max_dist,
            output_csv,
            error_file
        )
        write.csv(result, file = output_csv, append = TRUE, row.names = FALSE)
        final_result <- "Finished all points"
    }

    return(final_result)
}
||||||| 8bdf5d4
=======
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
>>>>>>> 9bace6d764915991d850da432164d016ac4f4903
