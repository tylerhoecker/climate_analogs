using DataFrames
using DataFramesMeta
using Distributions
using Distances
using Suppressor
using Random
using CSV
using Statistics
using StatsBase
using LinearAlgebra


#function for inverse covariance matrix
function calculate_covariance_matrix(pt_i::Integer, focal_data_cov::Vector{Any})
    mat::Matrix{Int16} = DataFrame([focal_data_cov[i][pt_i, :] for i in 1:length(focal_data_cov)]) |>
        #select only the variables we care about
        df -> select(df, Not(:x, :y)) |>
        Matrix
    #check if a column is all the same value, find 
    check_cols = [all(mat[:,i] .== mat[1,i]) for i in 1:size(mat,2)]
    

        # compute the covariance matrix
    cov_i::Matrix{Float32} = cov(mat)

    # if any columns are all the same value, add Tikhonov regularization
    if any(check_cols)
        cov_i = cov_i + 0.0001 * I
    end

    inv_i::Matrix{Float32} = inv(cov_i)
    
    return inv_i
end
#function for Mahalanobis distance
function calculate_sqmahalanobis_distance(
    cov_i::AbstractMatrix{T1},
    analog_mat::AbstractMatrix{T2},
    focal_data_mean_i::AbstractMatrix{T3}) where {T1 <: Real, T2 <: Real, T3 <: Real}
    
        d::Vector{Float32} = pairwise(SqMahalanobis(cov_i, skipchecks = true), analog_mat, focal_data_mean_i, dims = 1) |> 
        vec 
    
    return d
end

function process_mahalanobis_data(
    out_dt::DataFrame, 
    f_x::AbstractFloat,
    f_y::AbstractFloat,
    min_dist::Union{AbstractFloat, Integer},
    n_analog_use::Integer,
    focal_data_mean_i::Matrix{Float32})

    # Calculate great circle distance and round to 1 decimal place
    out_dt.dist_km = round.(great_circle_distance.(out_dt.a_y, out_dt.a_x, f_y, f_x), digits=1)
    
    # Subset the DataFrame based on the minimum distance
    subset!(out_dt, :dist_km => ByRow(>(min_dist)))
    
    # Sort the DataFrame by Mahalanobis distance
    sort!(out_dt, :md)
    
    # Select the top n_analog_use rows
    out_dt = out_dt[1:n_analog_use, :]
    
    # Calculate sigma and round to 4 decimal places
    out_dt.sigma = calc_sigma.(out_dt.md, length(focal_data_mean_i)) |> df -> round.(df, digits=4)
    
    # Round Mahalanobis distance to 3 decimal places
    out_dt.md = round.(out_dt.md, digits=3)
    
    # Repeat f_x and f_y for each row in the DataFrame
    out_dt.f_x = repeat([f_x], nrow(out_dt))
    out_dt.f_y = repeat([f_y], nrow(out_dt))
    
    return out_dt
end


function calc_mahalanobis(
    pt_i::Integer,
    focal_data_cov::Vector{Any},
    focal_data_mean::DataFrame,
    analog_data::DataFrame,
    var_names::Vector{String},
    n_analog_pool::Integer,
    n_analog_use::Integer,
    min_dist::Union{AbstractFloat, Integer},
)
    # Build cov matrix for pt_i from 30 years of annual projected future data
    

    cov_i = calculate_covariance_matrix(pt_i, focal_data_cov)
    
    
    

    # Calculate the mean of the future annuals -----------------------------------
    # Option to supply these mean focal data rather than derive (in the case of contemporary validation)
    if size(focal_data_mean,2) == 1
        focal_data_mean_i = focal_data_mean[pt_i][var_names]
    else
        focal_data_mean_i::Matrix{Float32} = DataFrame([focal_data_cov[i][pt_i,:] for i in 1:length(focal_data_cov)]) |>
            df -> select(df, Not(:x, :y)) |>
            df -> mean(Matrix(df), dims=1) 
    end

    # Analog pool (from historical normals) ------------------------------------
    # Random global sample == n_analog_pool
    #create a random sample of n_analog_pool numbers from 1:nrow(analog_data)
    indicies::Vector{UInt32} =sample(1:nrow(analog_data), n_analog_pool, replace=false)[1:n_analog_pool] 
        
    
    random_pts= analog_data[indicies, :]
    analog_mat = random_pts 
    #select only the variables we care about
    analog_mat = analog_mat[:, var_names]
    analog_mat = Matrix(analog_mat) 
   
    # Calculate the sigma dissimilarity -----------------------------------------


    # Sigma dissimilarity between pt_i and analog pool ---------------------------
    # Calculate Mahalanobis distances from analog mat to focal_data_mean
    
    
    
    d = @suppress calculate_sqmahalanobis_distance(cov_i, analog_mat, focal_data_mean_i)
    cov_i = nothing
    analog_mat = nothing

    f_x = focal_data_cov[1][pt_i,"x"]
    f_y = focal_data_cov[1][pt_i,"y"]
    # Save output
    
    out_dt = DataFrame(
        "a_x" => random_pts[:, "x"],
        "a_y" => random_pts[:, "y"],
        "md" => d,
        #preallocate
        "dist_km" => zeros(Float32, n_analog_pool),
        "sigma" => zeros(Float32, n_analog_pool),
        "f_x" => zeros(Float32, n_analog_pool),
        "f_y" => zeros(Float32, n_analog_pool)
    )

   final_out_dt = process_mahalanobis_data(out_dt, f_x, f_y, min_dist, n_analog_use, focal_data_mean_i) 
   f_x = nothing
   f_y = nothing
   min_dist = nothing
   n_analog_use = nothing
   focal_data_mean_i = nothing
    
    return final_out_dt
end

