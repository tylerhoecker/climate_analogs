#' @description
#' # Function to transform squared Mahalanobis distances into standard deviations of a chi-squared distribution
#' Adapted from Mahony et al. 2017 GCB: https://doi.org/10.1111/gcb.13645

#' Mahony follows a procedure where D is immediately unsquared,
#' then is fed into a chi distribution function.
#' Mathematically, it is equivalent to feed the squared-D into
#' a chi-squared distribution function, and then unsquare the result.
#' stats::chisq is an order of magnitude faster than chi::qchi.

#' @param d A vector of squared Mahalanobis distances
#' @param dimensions An integer indicating the degrees of freedom in the Chi-Sq distribution,
#' #' equal to the number of variables used in the MD calculation
#' @return A vector of "sigmas" - unsquared standard deviations from a chi-squared distribution


using Distributions
function calc_sigma(d, dimensions)
    # Convert distances to percentiles of chi-squared distribution (mulit-dimensional normal)
    # df = number of dimensions / climate variables
    p = cdf(Chisq(dimensions), d)
    # Convert percentiles into quantiles (standard deviations from mean)
    sigma = quantile(Chisq(1), p) # df is now 1
    # Here, take square root to unsquare the "distances"
    sigma = sqrt(sigma) 

    return sigma
end
