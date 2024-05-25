# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]]
# Core Mahalanobis distance -> sigma calculation
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

    n_analog_pool = 1000000
    tic()
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
    toc()

    # Save output
    out_dt <- data.table(
        "f_x" = focal_data_cov[[1]][pt_i][["x"]], 
        "f_y" = focal_data_cov[[1]][pt_i][["y"]],
        "a_x" = analog_data[random_pts][["x"]],
        "a_y" = analog_data[random_pts][["y"]],
        "md" = d
    )

    out_dt <- out_dt |>
        # Euclidean/geographic distance between focal point and analogs
        dplyr::mutate("dist_km" = round(sqrt((f_x - a_x)^2 + (f_y - a_y)^2)*0.001, 1)) |>
        dplyr::filter(dist_km > min_dist) |>
        # Sort all analogs by sigma
        dplyr::arrange(md) |>
        # Pick n analogs from pool with lowest sigma value
        dplyr::slice_head(n = n_analog_use) |>
        # Other speed improvement: only calculate sigma for the analogs we're keeping!
        # Can order the best 1000 based on D... relative order is the same as sigma.
        dplyr::mutate(sigma = calc_sigma(md, length(focal_data_mean))) |>
        dplyr::mutate(sigma = round(sigma, 4),
                      md = round(md, 3))

    rm(cov_i, focal_data_mean, random_pts, analog_mat, d)
    # This takes ~50% of compute time... not sure if neccesary to prevent memory leakge
    gc()
    return(out_dt)
}
