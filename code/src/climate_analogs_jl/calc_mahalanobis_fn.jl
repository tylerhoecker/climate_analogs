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
using ThreadPools: @bthreads

#fast sort from <https://discourse.julialang.org/t/faster-sorting-of-dataframes-jl-via-background-ordering/45084>
assignnew!(v, idx) = v .= v[idx]
function fsort!(dataframe, by)
	idx = sortperm(dataframe[!, by])

	@bthreads for n in names(dataframe)
		assignnew!(dataframe[!, n], idx)
	end
	dataframe
end


#function for inverse covariance matrix
function calculate_covariance_matrix(pt_i::Integer, focal_data_cov::Union{Vector{Any}, Vector{DataFrame}})
	mat::Matrix{Int16} =
		DataFrame([focal_data_cov[i][pt_i, :] for i in 1:length(focal_data_cov)]) |>
		#select only the variables we care about
		df -> select(df, Not(:x, :y)) |>
			  Matrix
	#check if a column is all the same value, find 
	check_cols = [all(mat[:, i] .== mat[1, i]) for i in 1:size(mat, 2)]


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
	insert_df::DataFrame,
	f_x::AbstractFloat,
	f_y::AbstractFloat,
	min_dist::Union{AbstractFloat, Integer},
	max_dist::Union{AbstractFloat, Integer},
	n_analog_use::Integer,
	focal_data_mean_i::Matrix{Float32},
)

	# Calculate great circle distance and round to 1 decimal place

	insert_df.dist_km = round.(great_circle_distance(insert_df.a_y, insert_df.a_x, f_y, f_x), digits = 1) #24.1ms

	# Subset the DataFrame based on the minimum distance
	subset!(insert_df, :dist_km => ByRow(>(min_dist))) # 5 ms

	# Subset the DataFrame based on maximum distance
	subset!(insert_df, :dist_km => ByRow(<(max_dist))) # 5 ms

	# Sort the DataFrame by Mahalanobis distance
	fsort!(insert_df, :md) #5 ms 17 allocations

	# Select the top n_analog_use rows
	@inbounds result_i_sub = insert_df[1:n_analog_use, :] # 1.9 μs 26 allocations
	insert_df = nothing
	# Calculate sigma and round to 4 decimal places
	result_i_sub.sigma = calc_sigma.(result_i_sub.md, length(focal_data_mean_i)) |> df -> round.(df, digits = 4) #38 μ s 104 allocations

	# Round Mahalanobis distance to 3 decimal places
	result_i_sub.md = round.(result_i_sub.md, digits = 3)

	# Repeat f_x and f_y for each row in the DataFrame
	result_i_sub.f_x = repeat([f_x], nrow(result_i_sub))
	result_i_sub.f_y = repeat([f_y], nrow(result_i_sub))

	return result_i_sub
end


function calc_mahalanobis(
	pt_i::Integer,
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	focal_data_mean::DataFrame,
	var_names::Vector{String},
	n_analog_use::Integer,
	min_dist::Union{AbstractFloat, Integer},
	max_dist::Union{AbstractFloat, Integer},
	input_sample::DataFrame,
	insert_df::DataFrame,
)
	# Build cov matrix for pt_i from 30 years of annual projected future data


	cov_i = calculate_covariance_matrix(pt_i, focal_data_cov) #50 μs 306 allocations




	# Calculate the mean of the future annuals -----------------------------------
	# Option to supply these mean focal data rather than derive (in the case of contemporary validation)
	#44μs 288 allocations
	if size(focal_data_mean, 2) == 1
		focal_data_mean_i = focal_data_mean[pt_i][var_names]
	else
		focal_data_mean_i::Matrix{Float32} = DataFrame([focal_data_cov[i][pt_i, :] for i in 1:length(focal_data_cov)]) |>
											 df -> select(df, Not(:x, :y)) |>
												   df -> mean(Matrix(df), dims = 1)
	end



	analog_mat = Matrix(input_sample[:, var_names]) #910 μs 10 allocations / 1.9 ms 37 allocations

	# Calculate the sigma dissimilarity -----------------------------------------


	# Sigma dissimilarity between pt_i and analog pool ---------------------------
	# Calculate Mahalanobis distances from analog mat to focal_data_mean



	insert_df.md = @suppress calculate_sqmahalanobis_distance(cov_i, analog_mat, focal_data_mean_i) #32 ms 249 allocations
	cov_i = nothing
	analog_mat = nothing

	f_x = focal_data_cov[1][pt_i, "x"] # 100 ns 1 allocation
	f_y = focal_data_cov[1][pt_i, "y"] # 123 ns 1 allocation
	# Save output

	# script to here is 35.188 ms 881 allocations

	result_i_final = process_mahalanobis_data(insert_df, f_x, f_y, min_dist, max_dist, n_analog_use, focal_data_mean_i) # 32 ms 394 allocations

	f_x = nothing
	f_y = nothing
	min_dist = nothing
	n_analog_use = nothing

	return result_i_final
end

