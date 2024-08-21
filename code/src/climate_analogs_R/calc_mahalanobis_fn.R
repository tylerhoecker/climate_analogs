<<<<<<< HEAD
# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]]
# Core Mahalanobis distance -> sigma calculation
library(tictoc)
library(purrr)
library(data.table)
library(dplyr)
source("code/calc_sigma_fn.R")

# great circle distance between two points
# Function for inverse covariance matrix
calculate_covariance_matrix <- function(pt_i, focal_data_cov) {
    # TODO
    mat <- as.matrix(data.frame(focal_data_cov[[pt_i]]))
    mat <- mat[, !(colnames(mat) %in% c("x", "y"))]

    # Check if a column is all the same value
    check_cols <- sapply(mat, function(x) all(x == x[1]))

    # Compute the covariance matrix
    cov_i <- cov(mat)

    # If any columns are all the same value, add Tikhonov regularization
    if (any(check_cols)) {
        cov_i <- cov_i + 0.0001 * diag(ncol(cov_i))
    }

    inv_i <- solve(cov_i)

    return(inv_i)
}

# Function for Mahalanobis distance
calculate_sqmahalanobis_distance <- function(cov_i, analog_mat, focal_data_mean_i) {
    # TODO
    d <- purrr::map2_dbl(analog_mat, focal_data_mean_i, function(x, y) {
        mahalanobis(x, y, cov_i, inverted = TRUE)
    })

    return(d)
}

# Function to process Mahalanobis data
process_mahalanobis_data <- function(out_dt, f_x, f_y, min_dist, n_analog_use, focal_data_mean_i) {
    # Calculate great circle distance and round to 1 decimal place
    out_dt$dist_km <- round(great_circle_distance(cbind(out_dt$a_y, out_dt$a_x), c(f_y, f_x)), digits = 1)

    # Subset the DataFrame based on the minimum distance
    out_dt <- out_dt[out_dt$dist_km > min_dist, ]

    # Sort the DataFrame by Mahalanobis distance
    out_dt <- out_dt[order(out_dt$md), ]

    # Select the top n_analog_use rows
    out_dt <- out_dt[1:n_analog_use, ]

    # Calculate sigma and round to 4 decimal places
    out_dt$sigma <- round(calc_sigma(out_dt$md, length(focal_data_mean_i)), digits = 4)

    # Round Mahalanobis distance to 3 decimal places
    out_dt$md <- round(out_dt$md, digits = 3)

    # Repeat f_x and f_y for each row in the DataFrame
    out_dt$f_x <- rep(f_x, nrow(out_dt))
    out_dt$f_y <- rep(f_y, nrow(out_dt))

    return(out_dt)
}

# Function to calculate Mahalanobis
calc_mahalanobis <- function(pt_i, focal_data_cov, focal_data_mean, analog_data, var_names, n_analog_pool, n_analog_use, min_dist) {
    # Build cov matrix for pt_i from 30 years of annual projected future data
    cov_i <- calculate_covariance_matrix(pt_i, focal_data_cov)

    # Calculate the mean of the future annuals
    if (ncol(focal_data_mean) == 1) {
        focal_data_mean_i <- focal_data_mean[pt_i, var_names]
    } else {
        focal_data_mean_i <- colMeans(as.matrix(data.frame(focal_data_cov[[pt_i]]))[, !(colnames(mat) %in% c("x", "y"))])
    }

    # Analog pool (from historical normals)
    indicies <- sample(1:nrow(analog_data), n_analog_pool, replace = FALSE)[1:n_analog_pool]
    random_pts <- analog_data[indicies, ]
    analog_mat <- random_pts[, var_names]

    # Calculate the sigma dissimilarity
    analog_mat <- as.matrix(analog_mat)

    # Sigma dissimilarity between pt_i and analog pool
    d <- calculate_sqmahalanobis_distance(cov_i, analog_mat, focal_data_mean_i)

    f_x <- focal_data_cov[[1]][pt_i, "x"]
    f_y <- focal_data_cov[[1]][pt_i, "y"]

    # Save output
    out_dt <- data.frame(
        a_x = random_pts$x,
        a_y = random_pts$y,
        md = d,
        dist_km = rep(0, n_analog_pool),
        sigma = rep(0, n_analog_pool),
        f_x = rep(0, n_analog_pool),
        f_y = rep(0, n_analog_pool)
    )

    final_out_dt <- process_mahalanobis_data(out_dt, f_x, f_y, min_dist, n_analog_use, focal_data_mean_i)

    return(final_out_dt)
}
||||||| 8bdf5d4
=======
# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]]
# Core Mahalanobis distance -> sigma calculation
library(tictoc)
library(purrr)
library(data.table)
library(dplyr)
source("code/calc_sigma_fn.R")

# great circle distance between two points

great_circle_distance <- function(lat1, lon1, lat2, lon2) {
    # Convert latitude and longitude from degrees to radians
    lat1 <- lat1 * pi / 180
    lon1 <- lon1 * pi / 180
    lat2 <- lat2 * pi / 180
    lon2 <- lon2 * pi / 180

    # Calculate the change in latitude and longitude
    dlat <- lat2 - lat1
    dlon <- lon2 - lon1

    # Calculate the great circle distance
    a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
    c <- 2 * asin(sqrt(a))
    R <- 6371  # Earth's radius in km
    return(R * c)
}

calc_mahalanobis <- function(
    pt_i,
    focal_data_cov,
    focal_data_mean,
    analog_data,
    var_names,
    n_analog_pool,
    n_analog_use,
    min_dist
) {
    # Build cov matrix for pt_i from 30 years of annual projected future data
    cov_i <- purrr::map(focal_data_cov, ~ .x[pt_i]) |>
        rbindlist() |>
        dplyr::select(-x, -y) |>
        cov()

    # Calculate the mean of the future annuals -----------------------------------
    # Option to supply these mean focal data rather than derive (in the case of contemporary validation)
    if (length(focal_data_mean) == 1 ) {
        focal_data_mean <- focal_data_mean[pt_i] |>
            select(all_of(var_names)) |>
            unlist()
    } else {
        focal_data_mean <- purrr::map(focal_data_cov, ~ .x[pt_i]) |>
            rbindlist() |>
            # Calculate mean of 30 years for every point (row in data.table)
            dplyr::summarise(across(all_of(var_names), mean)) |>
            unlist()
    }

   
    
    # Analog pool (from historical normals) ------------------------------------
    # Random global sample == n_analog_pool
    random_pts <- seq_len(nrow(analog_data)) |>
        sample(size = n_analog_pool)

    # Build reference matrix from random sample of historical normals
    analog_mat <- analog_data[random_pts, ] |>
        dplyr::select(-x, -y) |>
        as.matrix()

    # Sigma dissimilarity between pt_i and analog pool ---------------------------
    # Calculate Mahalanobis distances
    d <- mahalanobis(analog_mat, focal_data_mean, cov_i)
    

    f_x = focal_data_cov[[1]][pt_i][["x"]]
    f_y = focal_data_cov[[1]][pt_i][["y"]]
    # Save output
    out_dt <- data.table(
        "a_x" = analog_data[random_pts][["x"]],
        "a_y" = analog_data[random_pts][["y"]],
        "md" = d
    )

    out_dt <- out_dt[, dist_km := round(great_circle_distance(a_y, a_x, f_y, f_x), 1)][
        dist_km > min_dist][
        order(md)][
        out_dt[, .I[seq_len(n_analog_use)]]][
        , sigma := calc_sigma(md, length(focal_data_mean))][
        , sigma := round(sigma, 4)][
        , md := round(md, 3)][
        , f_x := f_x][
        , f_y := f_y
        ]
    # the above is equivalent to:

    # out_dt <- out_dt |>
    #     # Euclidean/geographic distance between focal point and analogs
    #     dplyr::mutate("dist_km" = round(great_circle_distance(a_y,a_x,f_y,f_x), 1)) |>
    #     dplyr::filter(dist_km > min_dist) |>
    #     # Sort all analogs by sigma
    #     dplyr::arrange(md) |>
    #     # Pick n analogs from pool with lowest sigma value
    #     dplyr::slice_head(n = n_analog_use) |>
    #     # Other speed improvement: only calculate sigma for the analogs we're keeping!
    #     # Can order the best 1000 based on D... relative order is the same as sigma.
    #     dplyr::mutate(sigma = calc_sigma(md, length(focal_data_mean))) |>
    #     dplyr::mutate(sigma = round(sigma, 4),
    #                   md = round(md, 3))
    rm(cov_i, focal_data_mean, random_pts, analog_mat, d)
    # This takes ~50% of compute time... not sure if neccesary to prevent memory leakge
    gc()
    return(out_dt)
}
>>>>>>> 9bace6d764915991d850da432164d016ac4f4903
