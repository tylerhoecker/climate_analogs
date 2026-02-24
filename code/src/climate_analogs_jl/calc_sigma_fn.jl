using Distributions
"""
  calc_sigma_fn(d::Vector, dimensions::Int)::Vector

Transform squared Mahalanobis distances into standard deviations of a chi-squared distribution.

Adapted from Mahony et al. 2017 GCB: https://doi.org/10.1111/gcb.13645

# Details

Mahony follows a procedure where D is immediately unsquared, then is fed into a chi distribution
function. Mathematically, it is equivalent to feed the squared-D into a chi-squared distribution
function, and then unsquare the result. `Distributions.chisq` is an order of magnitude faster
than alternative chi distribution functions.

# Arguments

- `d::Vector`: A vector of squared Mahalanobis distances
- `dimensions::Int`: Degrees of freedom in the Chi-squared distribution, equal to the number
  of variables used in the Mahalanobis distance calculation

# Returns

- `Vector`: A vector of "sigmas" - unsquared standard deviations from a chi-squared distribution
"""
function calc_sigma(d::AbstractFloat, dimensions::Int)
  # Convert distances to percentiles of chi-squared distribution (mulit-dimensional normal)
  # df = number of dimensions / climate variables
  p = cdf(Chisq(dimensions), d)
  # Convert percentiles into quantiles (standard deviations from mean)
  sigma = quantile(Chisq(1), p) # df is now 1
  # Here, take square root to unsquare the "distances"
  sigma = sqrt(sigma) |> Float32

  return sigma
end
