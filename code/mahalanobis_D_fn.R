# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]] 
# Core Mahalanobis distance -> sigma calculation
md_fun <- function(pt_i, 
                   focal_data_cov,
                   focal_data_mean, 
                   analog_data, 
                   var_names, 
                   n_analog_pool, 
                   n_analog_keep, 
                   min_dist) {
        
        # Build cov matrix for pt_i from 30 years of annual projected future data
        cov_i <- map(focal_data_cov, ~ .x[pt_i]) |>
            rbindlist() |>
            dplyr::select(-x, -y) |>
            cov()

        print(cov_i)

        # Calculate the mean of the future annuals -----------------------------------
        # This provides different result than calculating the mean here...
        # probably a result of the downscaling process
        # Option to supply these mean focal data rather than derive (in the case of contemporary validation)
        if(is.data.table(focal_data_mean)){
            focal_mean <- focal_data_mean[pt_i] |>
                select(all_of(var_names)) |>
                unlist()
        } else {
            focal_mean <- map(focal_data_cov, ~ .x[pt_i]) |>
            rbindlist() |>
            summarise(across(all_of(var_names), mean)) |>
            unlist()
        }

        print(focal_mean)

        # Analog pool (from historical normals) --------------------------------------
        # This could be integrated with the chunk below, using slice_sample,
        # but saving this index will be useful for extracting coordinates later, which we
        # don't want in the analog_mat
        random_pts <- seq_len(nrow(analog_data)) |>
            sample(size = n_analog_pool)
        
        # Build reference matrix from random sample of historical normals
        analog_mat <- analog_data[random_pts, ] |>
            dplyr::select(-x, -y) |>
            as.matrix()
        
        print(analog_mat)

        # Sigma dissimilarity between pt_i and analog pool ---------------------------
        # Calculate Mahalanobis distances
        d <- mahalanobis(analog_mat, focal_mean, cov_i)
        
        # Mahoney follows a procedure where D is immediately unsquared,
        # then is fed into a chi distribution function.
        # Mathematically, it is equivalent to feed the squared-D into
        # a chi-squared distribution function, and then unsquare the result. 
        # stats::chisq is an order of magnitude faster than chi::qchi.
        
        # Convert distances to percentiles of chi distribution (mulit-dimensional normal)
        # df = number of dimensions / climate variables
        p <- pchisq(d, df = length(focal_mean))
        # Convert percentiles into quantiles (standard deviations from mean)
        sigma <- qchisq(p, df = 1) # df is now 1...
        # Here, take square root to unsquare the distances. Round to 3 sig figs, for space. 
        sigma <- sqrt(sigma) |>
            round(4)

        # Save output
        out_dt <- data.table("f_x" = focal_data_mean[pt_i][["x"]],
                            "f_y" = focal_data_mean[pt_i][["y"]],
                            "a_x" = analog_data[random_pts][["x"]],
                            "a_y" = analog_data[random_pts][["y"]],
                            "sigma" = sigma)

        out_dt <- out_dt |>
            dplyr::filter(is.finite(sigma), sigma < 2) |>
            # Euclidean distance between focal point and analogs
            mutate("dist_m" = round(sqrt((f_x - a_x)^2 + (f_y - a_y)^2))) |>
            filter(dist_m > min_dist) |>
            # Sort all analogs by sigma
            arrange(sigma) |>
            # Pick n analogs from pool with lowest sigma value
            slice_head(n = n_analog_keep)

        rm(cov_i, focal_mean, random_pts, analog_mat, d, p, sigma)
        gc()
        print(out_dt)
        return(out_dt)
        }    
# Wrapper function for a given dataframe of focal locations and climate data
find_analogs <- function(focal_data_cov, # Repeated observations from focal locations (ie, annualized futures) to calculate covariance matrix
                         focal_data_mean, # Mean of repeated observations from focal locations (not required, usually calculated from focal_data_cov)
                         analog_data, # Mean observations from potential analogs (ie historical means from across West)
                         var_names, # Climate variable names; expected to be names for columns in data.tables
                         n_analog_pool, # Size of global analog pool, sampled randomly from extent of analog_data 
                         n_analog_keep, # Number of good analogs to keep
                         min_dist, # Minimum distance to analogs in meters (for contemporary validation, exlcude nearby points)
                         output_dir, # Output directory to save resulting data.table as Rds
                         use_futures, # Do parallel processing with futures/furrr?
                         n_futures) { # How many workers to run in parallel?
   

    # Map function over all points in supplied dataset
    if (use_futures == TRUE) {
        options(future.globals.maxSize = 5000000000)
        plan(multicore, workers = n_futures)
        sigma_dt <- seq_len(nrow(focal_data_cov[[1]])) |>
            future_map_dfr(md_fun, .id = "focal_id", .progress = TRUE,
                           focal_data_cov = focal_data_cov,
                           focal_data_mean = focal_data_mean,
                           analog_data = analog_data,
                           var_names = var_names,
                           n_analog_pool = n_analog_pool,
                           n_analog_keep = n_analog_keep,
                           min_dist = min_dist)
    } else {
        sigma_dt <- seq_len(nrow(focal_data_cov[[1]])) |>
            map_dfr(md_fun, .id = "focal_id", .progress = TRUE,
                    focal_data_cov = focal_data_cov,
                    focal_data_mean = focal_data_mean,
                    analog_data = analog_data,
                    var_names = var_names,
                    n_analog_pool = n_analog_pool,
                    n_analog_keep = n_analog_keep,
                    min_dist = min_dist)
    }
    
    # Save as RDS
    saveRDS(sigma_dt, paste0(output_dir, "/random_tiles", ".Rds"))
    
    # Explicitly remove and free memory each iteration
    rm(sigma_dt)
    gc()                        
                         }





