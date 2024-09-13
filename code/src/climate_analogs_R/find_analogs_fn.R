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
source("code/src/climate_analogs_R/calc_mahalanobis_fn.R")
source("code/src/climate_analogs_R/calc_sigma_fn.R")
source("code/src/climate_analogs_R/spatial_partition_fn.R")
source("code/src/climate_analogs_R/geography_fns.R")
library(future)
library(furrr)
library(data.table)
library(dplyr)
library(purrr)
library(progressr)

# Function to sample analogs
sample_analogs <- function(filtered_analog_data, n_analog_pool) {
    # Create a random sample of n_analog_pool indices from 1:nrow(filtered_analog_data)
    indices <- sample(1:nrow(filtered_analog_data), n_analog_pool, replace = FALSE)
    return(filtered_analog_data[indices, ])
}

# Function to create a bit vector
create_bitVector <- function(analog_data, coords) {
    # Extract coordinates from the data frame
    dfx <- analog_data$x
    dfy <- analog_data$y

    # Extract bounds from the coordinates dictionary
    south <- coords[["south"]][1]
    north <- coords[["north"]][1]
    west <- coords[["west"]][1]
    east <- coords[["east"]][1]

    # Create a bit vector using vectorized operations
    bit_vector <- (dfy >= south & dfy <= north) & (dfx >= west & dfx <= east)

    return(bit_vector)
}

calculate_analogs <- function(
    focal_data_cov,
    focal_data_mean,
    analog_data,
    var_names,
    n_analog_pool,
    n_analog_use,
    min_dist,
    max_dist,
    x) {
    if (max_dist != Inf) {
        coords <- max_distance_coordinates(focal_data_cov[[1]][x, "y"], focal_data_cov[[1]][x, "x"], max_dist)
        # 1.4 ms
        # Use the coordinates to filter the analog data
        bit_vector <- create_bitVector(analog_data, coords)
        # 236 ms
        # alternative to create_bitVector


        indices <- which(bit_vector)
        # 15 ms
        filtered_analog_data <- analog_data[indices, ] # 17 ms
    } else {
        filtered_analog_data <- analog_data
    }

    analog_data <- NULL

    sampled_analog_data <- sample_analogs(filtered_analog_data, n_analog_pool) # 50 ms

    # Preallocate
    insert_df <- data.table(
        a_x = sampled_analog_data$x,
        a_y = sampled_analog_data$y,
        md = rep(0, n_analog_pool),
        dist_km = rep(0, n_analog_pool)
    )

    result_i <- suppressWarnings(calc_mahalanobis(
        x,
        focal_data_cov,
        focal_data_mean,
        var_names,
        n_analog_use,
        min_dist,
        max_dist,
        sampled_analog_data,
        insert_df
    )) # 490 ms

    # Remove filtered_analog_data from memory
    filtered_analog_data <- NULL

    # Save the result
    # ...
    return(result_i)
}
# Function to calculate analogs distributed
calculate_analogs_distributed <- function(
    focal_data_cov,
    focal_data_mean,
    analog_data,
    var_names,
    n_analog_pool,
    n_analog_use,
    min_dist,
    max_dist,
    error_file,
    futures_nprocs = 1) {
    globals_size <- sum(sapply(.GlobalEnv, object.size))
    furrr_options(
        seed = TRUE
    )
    options(future.globals.maxSize = globals_size + 5 * 1024^3)
    plan(multisession, workers = futures_nprocs) # Set up parallel processing

    # Set future.globals.maxSize to the size of globals plus 5 GB


    # Progress bar setup
    handlers(global = TRUE)
    with_progress({
        p <- progressor(along = 1:nrow(focal_data_cov[[1]]))

        results <- future_map_dfr(1:nrow(focal_data_cov[[1]]), function(x) {
            p()
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
                    ) # 500 ms

                    if (x %% 100 == 0 && available_memory() < 10 * 1024^2) {
                        write(paste("Broke at", x), file = error_file, append = TRUE)
                    }

                    result_i
                },
                error = function(e) {
                    write(as.character(x), file = error_file, append = TRUE)
                    data.frame()
                }
            )
        }, .options = furrr_options(seed = TRUE))
    })
    plan(sequential) # Reset to sequential processing
    return(results)
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
    output_file,
    futures_nprocs = 1) {
    # Define output CSV file
    output_csv <- paste0(output_file, ".csv")

    # Write headers to CSV if it does not exist
    if (!file.exists(output_csv)) {
        write_csv(data.frame(
            a_x = numeric(),
            a_y = numeric(),
            md = numeric(),
            dist_km = numeric(),
            sigma = numeric(),
            f_x = numeric(),
            f_y = numeric()
        ), output_csv)
    }

    # Replace error file
    error_file <- file.path("/project/umontana_climate_analogs/climate_analogs/data", paste0(basename(output_file), "_error.txt"))
    write("", error_file)

    # Check if the data fits in memory
    spatial_partition_result <- spatial_partition(
        focal_data_cov,
        focal_data_mean,
        analog_data,
        n_analog_use
    )

    tile_ids <- spatial_partition_result$tile_ids
    tile_extents <- spatial_partition_result$tile_extents

    if (!is.null(tile_ids)) {
        cat("splitting tiles:", max(tile_ids$tile_id), "\n")

        # Split the data into tiles
        n_tiles <- max(tile_ids$tile_id)
        for (tile in 1:n_tiles) {
            # Filter the analog pool to only include points within max_dist of the focal points in the tile
            filtered_cov_data <- filter_cov_data(
                focal_data_cov,
                tile_extents,
                tile
            )
            filtered_data_mean <- filter_mean_data(
                focal_data_mean,
                tile_extents,
                tile
            )
            filtered_analog_data <- filter_analog_pool(
                analog_data,
                tile_extents,
                tile,
                max_dist
            )

            # Run the calculate_analogs function
            result_i <- calculate_analogs_distributed(
                filtered_cov_data,
                filtered_data_mean,
                filtered_analog_data,
                var_names,
                n_analog_pool,
                n_analog_use,
                min_dist,
                max_dist,
                paste0(error_file, "_tile", tile),
                futures_nprocs
            )

            write_csv(result_i, paste0(output_csv, "_tile", tile))

            result_i <- NULL
            filtered_analog_data <- NULL
            filtered_cov_data <- NULL
            filtered_data_mean <- NULL

            cat("Finished tile", tile, "\n")
        }
        final_result <- "Finished all tiles"
    } else {
        result <- calculate_analogs_distributed(
            focal_data_cov,
            focal_data_mean,
            analog_data,
            var_names,
            n_analog_pool,
            n_analog_use,
            min_dist,
            max_dist,
            error_file,
            futures_nprocs
        )

        write_csv(result, output_csv, append = TRUE)
        final_result <- "Finished all points"
    }

    return(final_result)
}
