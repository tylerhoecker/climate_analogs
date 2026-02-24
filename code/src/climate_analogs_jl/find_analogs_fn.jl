using DataFrames
using DataFramesMeta
using ProgressMeter
using Suppressor
using CSV
using CodecZlib
using Base.Threads
using Distributed

# sample analogs from the filtered analog pool
function sample_analogs(filtered_analog_data::Union{DataFrame,SubDataFrame},
  n_analog_pool::Integer)
  #create a random sample of n_analog_pool numbers from 1:nrow(analog_data)
  indicies = UInt32[]  # Initialize an empty vector of UInt32
  resize!(indicies, n_analog_pool)  # Resize it to the desired length

  sample!(1:nrow(filtered_analog_data), indicies, replace=false) #all steps take 30.1 ms 7 allocations
  return filtered_analog_data[indicies, :] #30.5 ms alone
end

# fast method to find the indicies of data within a bounding box defined by max distance from the focal point
function create_bitVector(analog_data::DataFrame, coords::Dict{String,Vector{Float32}})
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
  focal_data_cov::Union{Vector{Any},Vector{DataFrame}},
  focal_data_mean::DataFrame,
  analog_data::DataFrame,
  var_names::Vector{String},
  n_analog_pool::Integer,
  n_analog_use::Integer,
  min_dist::Union{AbstractFloat,Integer},
  max_dist::Union{AbstractFloat,Integer},
  x::Integer)
  if max_dist != Inf
    # Get the coordinates of the bounding box around the focal point
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


  # sample n_analog_pool analogs from the filtered analog pool
  sampled_analog_data = sample_analogs(filtered_analog_data, n_analog_pool) #70 ms 32 allocations

  #preallocate
  insert_df = DataFrame(
    "a_x" => sampled_analog_data[:, "x"],
    "a_y" => sampled_analog_data[:, "y"],
    "md" => zeros(Float32, n_analog_pool),
    "dist_km" => zeros(Float32, n_analog_pool),
  ) #3.6ms ms 48 allocations
  # this is the main function that calculates the Mahalanobis distance, sigma, and filters by distance
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

  return result_i
end

# this distributes the calculate_analogs function across all points in the focal dataset
function calculate_analogs_distributed(
  focal_data_cov::Union{Vector{Any},Vector{DataFrame}},
  focal_data_mean::DataFrame,
  analog_data::DataFrame,
  var_names::Vector{String},
  n_analog_pool::Integer,
  n_analog_use::Integer,
  min_dist::Union{AbstractFloat,Integer},
  max_dist::Union{AbstractFloat,Integer},
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

      # if memory is too low, log the point at which it broke to the error file
      if x % 100 == 0 && Sys.free_memory() < 10e9
        open(error_file, "a") do f
          write(f, "Broke at $x\n")
        end
      end

      result_i
    catch e
      # generic error to log any point that causes an error, but continue the loop
      open(error_file, "a") do f
        write(f, "$x\n")
      end
      # empty dataframe if error
      result_i = DataFrame()
    end
  end
  return result
end

"""
  `find_analogs` is the main function that takes in focal data,
    analog data, and parameters for the analog search
    , and returns a dataframe of the best analogs for each focal point.

  # Arguments
  - `focal_data_cov::Union{Vector{Any},Vector{DataFrame}}`: A vector of DataFrames containing the repeated observations from focal locations (ie, annualized futures) to calculate covariance matrix
  - `focal_data_mean::DataFrame`: A DataFrame containing the mean of repeated observations from focal locations
  - `analog_data::DataFrame`: A DataFrame containing the mean observations from potential analog pixels
  - `var_names::Vector{String}`: A vector of strings containing the names of the climate variables to use in the analog search; expected to be names for columns in data.tables
  - `n_analog_pool::Integer`: The size of the sampled global analog pool
  - `n_analog_use::Integer`: The number of best analogs to keep for each focal point
  - `min_dist::Union{AbstractFloat,Integer}`: The minimum distance (in km) that an analog must be from the focal point
  - `max_dist::Union{AbstractFloat,Integer}`: The maximum distance (in km) that an analog can be from the focal point
  - `output_dir::String`: Directory where the output CSV files will be saved
  - `output_file::String`: The base name for the output CSV file
"""
function find_analogs(
  focal_data_cov::Union{Vector{Any},Vector{DataFrame}},
  focal_data_mean::DataFrame,
  analog_data::DataFrame,
  var_names::Vector{String},
  n_analog_pool::Integer,
  n_analog_use::Integer,
  min_dist::Union{AbstractFloat,Integer},
  max_dist::Union{AbstractFloat,Integer},
  output_dir::String,
  output_file::String)
  # Map function over all points in supplied dataset
  # write headers to CSV if csv does not exist
  output_name = joinpath(output_dir, output_file)
  output_gzip = joinpath(output_name) * ".csv.gz"
  #replace error file
  error_file = joinpath(output_dir, output_file * "_error.txt")

  # create error file
  open(error_file, "w") do f
    write(f, "")
  end


  # Check if the data fits in memory then partition it
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
    # if there are tiles, loop through each tile and run the calculate_analogs function on the subset of data in that tile, then write the results separate CSV
    for tile in 1:n_tiles
      # Filter the analog pool to only include points within max_dist of the tile
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
      open(GzipCompressorStream, output_file * "_tile$tile.csv.gz", "w") do stream
        CSV.write(stream, result_i, append=false)
      end
      result_i = nothing
      filtered_analog_data = nothing
      filtered_cov_data = nothing
      filtered_data_mean = nothing
      println("Finished tile $tile")
    end
    final_result = "Finished all tiles"
  else
    # otherwise, run the calculate_analogs function on the full dataset
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
      CSV.write(stream, result, append=false)
    end
    final_result = "Finished all points"
  end

  return final_result
end
