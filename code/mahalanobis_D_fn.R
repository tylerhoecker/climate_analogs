# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]] 

# Core Mahalanobis distance -> sigma calculation
md_fun <- function(pt_i) {
  
  # Build cov matrix for pt_i from 30 years of annual projected future data
  cov_i <- map(focal_data, ~ .x[pt_i]) |>
    rbindlist() |>
    dplyr::select(-x, -y) |>
    cov()

  # Calculate the mean of the future annuals -----------------------------------
  # This provides different result than calculating the mean here...
  # probably a result of the downscaling process
  focal_mean <- map(focal_data, ~ .x[pt_i]) |>
    rbindlist() |>
    summarise(across(var_names, mean)) |>
    unlist()
 
  # Analog pool - from historical normals --------------------------------------
  # This *could* be done outside the function, but this is less biased, because
  # random samples are drawn with each iteration, ultimately sampling from (probably)
  # nearly the entire landscape
  
  # This could be integrated with the chunk below, using slice_sample,
  # but saving this index will be useful for extracting coordinates later
  random_pts <- seq_len(nrow(analog_data)) |>
    sample(size = n_analog_pool)
  
  # Build reference matrix from random sample of historical normals
  analog_mat <- analog_data[random_pts, ] |>
    dplyr::select(-x, -y) |>
    as.matrix()
  
  # Sigma dissimilarity between pt_i and analog pool ---------------------------
  # Calculate 'raw' unsquared Mahalanobis distances
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
    round(3)

  # Save output
  out_dt <- data.table("f_x" = focal_data[[1]][pt_i, "x"],
                       "f_y" = focal_data[[1]][pt_i, "y"],
                       "a_x" = analog_data[random_pts, "x"],
                       "a_y" = analog_data[random_pts, "y"],
                       "sigma" = sigma)
  # Not sure why .x/.y get's appended above... data.table thing
  setnames(out_dt,
           old = c("f_x.x", "f_y.y", "a_x.x", "a_y.y"),
           new = c("f_x", "f_y", "a_x", "a_y"))

  out_dt <- out_dt |>
    dplyr::filter(is.finite(sigma), sigma < 2) |>
    # Euclidean distance between focal point and analogs
    mutate("dist_m" = round(sqrt((f_x - a_x)^2 + (f_y - a_y)^2))) |>
    # Sort all analogs by sigma
    arrange(sigma) |>
    # Pick n analogs from pool with lowest sigma value
    slice_head(n = n_analog_keep)

  rm(cov_i, focal_mean, random_pts, ref_mat, d, p, sigma)
  gc()
  return(out_dt)
}

# Wrapper function for a given dataframe of focal locations and climate data
find_analogs <- function(focal_data_list = chunk, 
                         analog_data = norm_hist, 
                         var_names = c('aet', 'def', 'tmax', 'tmin'),
                         n_analog_pool = 10000000,
                         n_analog_keep = 1000,
                         output_dir = "data/",
                         # Arguments for parallel processing
                         use_futures = TRUE,
                         n_futures = 4) {

    # Calculate sigmas to analogs for each focal point
    if(use_futures == TRUE){
        plan(multicore, workers = n_futures)
        sigma_dt <- seq_along(nrow(focal_data[[1]])) |>
            future_map_dfr(md_fun, .id = "focal_id", .progress = TRUE)
    } else {
        sigma_dt <- seq_len(nrow(focal_data[[1]])) |>
            map_dfr(md_fun, .id = "focal_id", .progress = TRUE)
    }
    
    # Save as RDS
    saveRDS(sigma_dt, paste0(output_dir, ".Rds"))
    
    # Explicitly remove and free memory each iteration
    rm(sigma_dt)
    gc()                        
                         }





