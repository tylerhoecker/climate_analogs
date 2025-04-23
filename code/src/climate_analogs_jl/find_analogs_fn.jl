using DataFrames
using DataFramesMeta
using ProgressMeter
using Suppressor
using CSV
using CodecZlib
using Base.Threads
using Distributed

function sample_analogs(filtered_analog_data::Union{DataFrame, SubDataFrame},
	n_analog_pool::Integer)
	#create a random sample of n_analog_pool numbers from 1:nrow(analog_data)
	indicies = UInt32[]  # Initialize an empty vector of UInt32
	resize!(indicies, n_analog_pool)  # Resize it to the desired length

	sample!(1:nrow(filtered_analog_data), indicies, replace = false) #all steps take 30.1 ms 7 allocations
	return filtered_analog_data[indicies, :] #30.5 ms alone
end

function create_bitVector(analog_data::DataFrame, coords::Dict{String, Vector{Float32}})
	# Ensure that the columns of analog_data are concrete arrays
	dfx::Vector{Float32} = analog_data.x
	dfy::Vector{Float32} = analog_data.y
	# Initialize BitVector with a size equal to the length of dfx
	bit_vector::BitVector = BitVector(undef, length(dfx))

	# Pre-extract bounds for the loop
	south::Float32 = coords["south"][1]
	north::Float32 = coords["north"][1]
	west::Float32 = coords["west"][1]
	east::Float32 = coords["east"][1]

	# Use @inbounds and make sure that types are consistent within the loop

	@inbounds for i::Int64 in eachindex(dfx)
		x::Float32 = dfx[i]
		y::Float32 = dfy[i]
		bit_vector[i]::Bool = (south <= y && y <= north) && (west <= x && x <= east)
	end

	return bit_vector

	# @. (analog_data.y >= coords["south"][1]) &
	#     (analog_data.y <= coords["north"][1]) &
	#     (analog_data.x >= coords["west"][1]) &
	#     (analog_data.x <= coords["east"][1])
end

function calculate_analogs(
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	focal_data_mean::DataFrame,
	analog_data::DataFrame,
	var_names::Vector{String},
	n_analog_pool::Integer,
	n_analog_use::Integer,
	min_dist::Union{AbstractFloat, Integer},
	max_dist::Union{AbstractFloat, Integer},
	x::Integer)
	if max_dist != Inf
		coords = max_distance_coordinates(focal_data_cov[1][x, "y"], focal_data_cov[1][x, "x"], max_dist)

		# Use the coordinates to filter the analog data using DataFrames

		indices = @inbounds findall(
			create_bitVector(analog_data, coords), #134 ms 3 allocations
		) # 196ms with 24 allocations just for this check #with for loop 148 ms 5 allocations



		filtered_analog_data = @view analog_data[indices, :] # this step is 3.884 ms 2 allocations

	#combined 275ms 51 allocations

	else
		filtered_analog_data = analog_data
	end
	analog_data = nothing


	sampled_analog_data = sample_analogs(filtered_analog_data, n_analog_pool) #70 ms 32 allocations

	#preallocate
	insert_df = DataFrame(
		"a_x" => sampled_analog_data[:, "x"],
		"a_y" => sampled_analog_data[:, "y"],
		"md" => zeros(Float32, n_analog_pool),
		"dist_km" => zeros(Float32, n_analog_pool),
	) #3.6ms ms 48 allocations 
	result_i = Suppressor.@suppress calc_mahalanobis(
		x,
		focal_data_cov,
		focal_data_mean,
		var_names,
		n_analog_use,
		min_dist,
		max_dist,
		sampled_analog_data,
		insert_df) #for 55-69: 120 ms 1400 allocations using fsort: sort is 213 ms 83,000 allocations

	# Remove filtered_analog_data from memory
	filtered_analog_data = nothing

	# Save the result
	# ...
	return result_i
end

function calculate_analogs_distributed(
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	focal_data_mean::DataFrame,
	analog_data::DataFrame,
	var_names::Vector{String},
	n_analog_pool::Integer,
	n_analog_use::Integer,
	min_dist::Union{AbstractFloat, Integer},
	max_dist::Union{AbstractFloat, Integer},
	error_file::String)

	result = ProgressMeter.@showprogress "Calculating Climate Analogs.." @distributed (vcat) for x in 1:size(focal_data_cov[1], 1)
		try
			result_i = calculate_analogs(
				focal_data_cov,
				focal_data_mean,
				analog_data,
				var_names,
				n_analog_pool,
				n_analog_use,
				min_dist,
				max_dist,
				x) #344 ms 1956 allocations

			if x % 100 == 0 && Sys.free_memory() < 10e9
				open(error_file, "a") do f
					write(f, "Broke at $x\n")
				end
			end

			result_i
		catch e
			open(error_file, "a") do f
				write(f, "$x\n")
			end
			result_i = DataFrame()
		end
	end
	return result
end

function find_analogs(
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	focal_data_mean::DataFrame,
	analog_data::DataFrame,
	var_names::Vector{String},
	n_analog_pool::Integer,
	n_analog_use::Integer,
	min_dist::Union{AbstractFloat, Integer},
	max_dist::Union{AbstractFloat, Integer},
	output_file::String)
	# Map function over all points in supplied dataset
	# write headers to CSV if csv does not exist
    output_gzip= joinpath(output_file) * ".csv.gz"
		#replace error file
	error_file = joinpath("/project/umontana_climate_analogs/climate_analogs/data", basename(output_file) * "_error.txt")

	open(error_file, "w") do f
		write(f, "")
	end


	# Check if the data fits in memory
	tile_ids, tile_extents = spatial_partition(
		focal_data_cov,
		focal_data_mean,
		analog_data,
		n_analog_use) # 148 μs 1696 allocations
	flush(stdout)

	if !isnothing(tile_ids)
		println("splitting tiles: $(maximum(tile_ids.tile_id))")
		# Split the data into tiles
		n_tiles = maximum(tile_ids.tile_id)
		for tile in 1:n_tiles
			# Filter the analog pool to only include points within max_dist of the focal points in the tile
			filtered_cov_data = filter_cov_data(
				focal_data_cov,
				tile_extents,
				tile)
			filtered_data_mean = filter_mean_data(
				focal_data_mean,
				tile_extents,
				tile)
			filtered_analog_data = filter_analog_pool(
				analog_data,
				tile_extents,
				tile,
				max_dist)
			# Run the calculate_analogs function
			result_i = calculate_analogs_distributed(
				filtered_cov_data,
				filtered_data_mean,
				filtered_analog_data,
				var_names,
				n_analog_pool,
				n_analog_use,
				min_dist,
				max_dist,
				error_file * "_tile$tile")
            open(GzipCompressorStream, output_file * "_tile$tile.csv.gz","w") do stream
    			CSV.write(stream, result_i, append = false)
            end
			result_i = nothing
			filtered_analog_data = nothing
			filtered_cov_data = nothing
			filtered_data_mean = nothing
			println("Finished tile $tile")
		end
		final_result = "Finished all tiles"
	else
		result = calculate_analogs_distributed(
			focal_data_cov,
			focal_data_mean,
			analog_data,
			var_names,
			n_analog_pool,
			n_analog_use,
			min_dist,
			max_dist,
			error_file)
        open(GzipCompressorStream, output_gzip, "w") do stream
            CSV.write(stream, result, append = false)
        end
		final_result = "Finished all points"
	end

	return final_result
end

