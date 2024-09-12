using DataFrames
using DataFramesMeta
using Distributed

struct NoMemory <: Exception
	msg::String
end

function check_memory(
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	focal_data_mean::DataFrame,
	analog_pool::DataFrame,
	n_analog_use::Integer,
)
	# check if SLURM environment variables exist and grab memory
	if haskey(ENV, "SLURM_CPUS_ON_NODE")

		if haskey(ENV, "SLURM_MEM_PER_NODE")
			total_memory = parse(Float64, ENV["SLURM_MEM_PER_NODE"]) *
						   1e-3
		else
			total_memory = parse(Float64, ENV["SLURM_CPUS_ON_NODE"]) *
						   parse(Float64, ENV["SLURM_MEM_PER_CPU"]) *
						   1e-3
		end
	else
		total_memory = Sys.total_memory() * 1e-9
	end
	used_memory = (Base.summarysize(focal_data_cov) +
				   Base.summarysize(focal_data_mean) +
				   Base.summarysize(analog_pool))

	sampling_memory = sizeof(Float32) * 7 * size(focal_data_mean, 1)

	nprocesses = Distributed.nprocs() * 2

	environment_memory = (used_memory + sampling_memory) * nprocesses
	output_memory = sizeof(Float32) * 7 * size(focal_data_mean, 1) * n_analog_use
	necessary_memory = (environment_memory + output_memory) * 1e-9 + 30
	if total_memory <= 30
		throw(NoMemory("Does not contain the minimal memory needed for overhead (30 GB)!"))
	end

	return necessary_memory, total_memory
end

function create_outer_extents(
	focal_data_coords::DataFrame,
	min_n_tiles::Int,
)
	# Get the minimum and maximum x and y coordinates
	min_x = minimum(focal_data_coords[:, "x"])
	max_x = maximum(focal_data_coords[:, "x"])
	min_y = minimum(focal_data_coords[:, "y"])
	max_y = maximum(focal_data_coords[:, "y"])

	# Calculate the width and height of the grid
	width = abs(max_x - min_x)
	height = abs(max_y - min_y)

	# Calculate the number of tiles in the x and y directions
	n_tiles_x = ceil(sqrt(min_n_tiles * width / height)) |> Int
	n_tiles_y = ceil(sqrt(min_n_tiles * height / width)) |> Int

	n_tiles = n_tiles_x * n_tiles_y
	# Calculate the width and height of each tile
	tile_width = width / n_tiles_x
	tile_height = height / n_tiles_y

	# Create a DataFrame to store the extents of each tile
	tile_extents = DataFrame(
		tile_id = 1:n_tiles,
		min_x = zeros(n_tiles),
		max_x = zeros(n_tiles),
		min_y = zeros(n_tiles),
		max_y = zeros(n_tiles),
	)

	# Calculate the extents of each tile
	for i in 1:n_tiles_x
		for j in 1:n_tiles_y
			tile_id = (i - 1) * n_tiles_y + j
			tile_extents[tile_id, "min_x"] = min_x + (i - 1) * tile_width
			tile_extents[tile_id, "max_x"] = min_x + i * tile_width
			tile_extents[tile_id, "min_y"] = min_y + (j - 1) * tile_height
			tile_extents[tile_id, "max_y"] = min_y + j * tile_height
		end
	end


	return tile_extents
end

function find_tile_id(
	x::AbstractFloat,
	y::AbstractFloat,
	tile_extents::DataFrame,
)
	# Find the tile that contains the point (x, y)
	for i in 1:nrow(tile_extents)
		if x >= tile_extents[i, "min_x"] && x <= tile_extents[i, "max_x"] &&
		   y >= tile_extents[i, "min_y"] && y <= tile_extents[i, "max_y"]
			return (tile_extents[i, "tile_id"] |> Int16)
		end
	end
end



function spatial_partition(
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	focal_data_mean::DataFrame,
	analog_pool::DataFrame,
	n_analog_use::Integer,
)
	# get memory needed for all inputs in GB 
	necessary_memory, total_memory = check_memory(
		focal_data_cov,
		focal_data_mean,
		analog_pool,
		n_analog_use,
	)

	if necessary_memory > total_memory
		println("Data does not fit in memory: necessary:$necessary_memory, total_memory: $total_memory")

		# calculate how many tiles we need to split the data into

		min_n_tiles = ceil(necessary_memory / total_memory) |> Int
		# Split the data into n_tiles using a spatial partitioning algorithm
		#Get the coordinates of the focal data
		focal_data_coords = DataFrame(
			x = focal_data_mean[:, "x"],
			y = focal_data_mean[:, "y"],
		)
		tile_extents = create_outer_extents(focal_data_coords, min_n_tiles)
		# Create a DataFrame to store the tile ID for each point
		tile_id = DataFrame(
			x = focal_data_coords[:, "x"],
			y = focal_data_coords[:, "y"],
			tile_id = zeros(Int16, nrow(focal_data_coords)),
		)
		# Assign each point to a tile
		for i in 1:nrow(focal_data_coords)
			x = focal_data_coords[i, "x"]
			y = focal_data_coords[i, "y"]
			tile_id[i, "tile_id"] = find_tile_id(x, y, tile_extents)
		end
		n_tiles = maximum(tile_extents.tile_id)
		println("Splitting data into $n_tiles tiles")




		return tile_id, tile_extents

	else
		println("Data fits in memory: necessary: $necessary_memory, total_memory: $total_memory")
		tile_id = nothing
		tile_extents = nothing
		return tile_id, tile_extents
	end
end

function filter_analog_pool(
	analog_pool::DataFrame,
	tile_extents::DataFrame,
	tile::Int,
	max_dist::Union{AbstractFloat, Integer},
)
	# Filter the analog pool to only include points within max_dist of the focal points in the tile
	min_x = max_distance_coordinates(tile_extents[tile, "min_y"], tile_extents[tile, "min_x"], max_dist)["west"][2]
	max_x = max_distance_coordinates(tile_extents[tile, "max_y"], tile_extents[tile, "max_x"], max_dist)["east"][2]
	min_y = max_distance_coordinates(tile_extents[tile, "min_y"], tile_extents[tile, "min_x"], max_dist)["south"][1]
	max_y = max_distance_coordinates(tile_extents[tile, "max_y"], tile_extents[tile, "max_x"], max_dist)["north"][1]

	analog_extents = DataFrame(
		min_x = min_x,
		max_x = max_x,
		min_y = min_y,
		max_y = max_y,
	)

	filtered_analog_pool = DataFramesMeta.@subset(analog_pool,
		(:y .>= min_y) .&
		(:y .<= max_y) .&
		(:x .>= min_x) .&
		(:x .<= max_x)
	)

	return filtered_analog_pool
end

function filter_cov_data(
	focal_data_cov::Union{Vector{Any}, Vector{DataFrame}},
	tile_id::DataFrame,
	tile::Int,
)
	# Filter the covariance data to only include points in the tile
	focal_data_cov_tile = copy(focal_data_cov)
	for i in 1:length(focal_data_cov)
		focal_data_cov_tile[i] = DataFramesMeta.@subset(focal_data_cov[i],
			(:x .>= tile_id[tile, "min_x"]) .&
			(:x .<= tile_id[tile, "max_x"]) .&
			(:y .>= tile_id[tile, "min_y"]) .&
			(:y .<= tile_id[tile, "max_y"])
		)
	end

	return focal_data_cov_tile
end

function filter_mean_data(
	focal_data_mean::DataFrame,
	tile_id::DataFrame,
	tile::Int,
)
	# Filter the mean data to only include points in the tile
	focal_data_mean_tile = DataFramesMeta.@subset(focal_data_mean,
		(:x .>= tile_id[tile, "min_x"]) .&
		(:x .<= tile_id[tile, "max_x"]) .&
		(:y .>= tile_id[tile, "min_y"]) .&
		(:y .<= tile_id[tile, "max_y"])
	)

	return focal_data_mean_tile
end


