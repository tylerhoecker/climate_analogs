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
