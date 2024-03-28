#' @description Wrapper function for a given dataframe of focal locations and climate data
#' @param focal_data_cov Repeated observations from focal locations (ie, annualized futures) to calculate covariance matrix
#' @param focal_data_mean Mean of repeated observations from focal locations (not required, usually calculated from focal_data_cov)
#' @param analog_data Mean observations from potential analogs (ie historical means from across West)
#' @param var_names Climate variable names; expected to be names for columns in data.tables
#' @param n_analog_pool Size of global analog pool, sampled randomly from extent of analog_data 
#' @param n_analog_use Number of good analogs to keep
#' @param min_dist Minimum distance to analogs in meters (for contemporary validation, exlcude nearby points)
#' @param output_dir Output directory to save resulting data.table as Rds
#' @param use_futures Do parallel processing with futures/furrr?
#' @param n_futures Number of workers to futures
#' #' equal to the number of variables used in the MD calculation
#' @return A vector of "sigmas" - unsquared standard deviations from a chi-squared distribution

find_analogs <- function(
    focal_data_cov,
    focal_data_mean,
    analog_data,
    var_names,
    n_analog_pool,
    n_analog_use,
    min_dist,
    output_dir,
    use_futures,
    n_futures
){
    # Map function over all points in supplied dataset
    if (use_futures == TRUE){
        # Set up futures
        options(future.globals.maxSize = 5000000000)
        plan(multisession, workers = n_futures)
        # Run over all points
        sigma_dt <- seq_len(nrow(focal_data_cov[[1]])) |>
            future_map_dfr(
                \(x){
                calc_mahalanobis(
                    pt_i = x, 
                    focal_data_cov,
                    focal_data_mean,
                    analog_data,
                    var_names,
                    n_analog_pool,
                    n_analog_use,
                    min_dist
                )    
            },  
            .id = "focal_id"
            )
        plan(sequential)
    } else {
        sigma_dt <- seq_len(nrow(focal_data_cov[[1]])) |>
            map_dfr(
                \(x){
                calc_mahalanobis(
                    pt_i = x, 
                    focal_data_cov,
                    focal_data_mean,
                    analog_data,
                    var_names,
                    n_analog_pool,
                    n_analog_use,
                    min_dist
                )    
            },  
            .id = "focal_id"
            )
    }
    
    # Save as RDS
    saveRDS(sigma_dt, paste0(output_dir, ".Rds"))
    
    # Explicitly remove and free memory each iteration
    rm(sigma_dt)
    gc()
}
