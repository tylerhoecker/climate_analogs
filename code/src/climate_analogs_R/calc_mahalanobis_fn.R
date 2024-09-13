# The MD/sigma dissimilarity function
# Testing: pt_i  <- 1000; chunk  <- chunks[[100]]
# Core Mahalanobis distance -> sigma calculation
library(tictoc)
library(purrr)
library(data.table)
library(dplyr)
source("code/src/climate_analogs_R/calc_sigma_fn.R")
source("code/src/climate_analogs_R/geography_fns.R")

calc_mahalanobis <- function(
    pt_i,
    focal_data_cov,
    focal_data_mean,
    var_names,
    n_analog_use,
    min_dist,
    max_dist,
    input_sample,
    insert_df) {
    # Build cov matrix for pt_i from 30 years of annual projected future data
    cov_i <- purrr::map(focal_data_cov, ~ .x[pt_i]) |>
        rbindlist() |>
        dplyr::select(-x, -y) |>
        cov()

    # Calculate the mean of the future annuals -----------------------------------
    # Option to supply these mean focal data rather than derive (in the case of contemporary validation)
    if (length(focal_data_mean) == 1) {
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




    # Build reference matrix from random sample of historical normals
    analog_mat <- input_sample[, var_names, env = list(var_names = as.list(var_names))] |>
        as.matrix()

    # Sigma dissimilarity between pt_i and analog pool ---------------------------
    # Calculate Mahalanobis distances
    d <- mahalanobis(analog_mat, focal_data_mean, cov_i)
    rm(cov_i, analog_mat)

    f_x <- focal_data_cov[[1]][pt_i][["x"]]
    f_y <- focal_data_cov[[1]][pt_i][["y"]]
    # Save output

    out_dt <- insert_df[, dist_km := round(great_circle_distance(a_y, a_x, f_y, f_x), 1)][
        dist_km > min_dist & dist_km < max_dist
    ][
        order(md)
    ][
        insert_df[, .I[seq_len(n_analog_use)]]
    ][
        , sigma := calc_sigma(md, length(focal_data_mean))
    ][
        , sigma := round(sigma, 4)
    ][
        , md := round(md, 3)
    ][
        , f_x := f_x
    ][
        , f_y := f_y
    ]
    # the above is equivalent to:

    # out_dt <- insert_df |>
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
    rm(focal_data_mean, input_sample, analog_mat, d)
    # This takes ~50% of compute time... not sure if neccesary to prevent memory leakge
    gc()
    return(out_dt)
}
