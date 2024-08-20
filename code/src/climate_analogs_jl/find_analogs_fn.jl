using DataFrames
using DataFramesMeta
using ProgressMeter
using Suppressor
using CSV

function calculate_analogs(
    focal_data_cov::Vector{Any},
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
    filtered_analog_data = DataFramesMeta.@subset(analog_data,
                                     (:y .>= coords["south"][1]) .&
                                      (:y .<= coords["north"][1]) .&
                                       (:x .>= coords["west"][2]) .&
                                        (:x .<= coords["east"][2]))
    else
      filtered_analog_data = analog_data
    end
    result_i = Suppressor.@suppress calc_mahalanobis(
        x,
        focal_data_cov,
        focal_data_mean,
        filtered_analog_data,
        var_names,
        n_analog_pool,
        n_analog_use,
        min_dist)
  
    # Remove filtered_analog_data from memory
    filtered_analog_data = nothing
  
    # Save the result
    # ...
    result_i
  end
  
  function run_calculate_analogs(
    focal_data_cov::Vector{Any},
    focal_data_mean::DataFrame,
    analog_data::DataFrame,
    var_names::Vector{String},
    n_analog_pool::Integer,
    n_analog_use::Integer,
    min_dist::Union{AbstractFloat, Integer},
    max_dist::Union{AbstractFloat, Integer},
    output_csv::String,
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
          x)
        
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
    focal_data_cov::Vector{Any},
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
    output_csv = joinpath(output_file) * ".csv"
    if !isfile(output_csv)
      CSV.write(output_csv, DataFrame(
        a_x = Float32[],
        a_y = Float32[],
        md = Float32[],
        dist_km = Float32[],
        sigma = Float32[],
        f_x = Float32[],
        f_y = Float32[],
      ), writeheader = true)
    end
#replace error file
 error_file = joinpath("/project/umontana_climate_analogs/climate_analogs/data",basename(output_file) * "_error.txt")
 if !isfile(error_file)
      open(error_file, "w") do f
        write(f, "")
      end
    end
    
   # Check if the data fits in memory
   tile_ids, tile_extents = spatial_partition(
                              focal_data_cov,
                              focal_data_mean,
                              analog_data,
                              n_analog_use)

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
        result_i = run_calculate_analogs(
            filtered_cov_data,
            filtered_data_mean,
            filtered_analog_data,
            var_names,
            n_analog_pool,
            n_analog_use,
            min_dist,
            max_dist,
            output_csv * "_tile$tile",
            error_file * "_tile$tile")
            CSV.write(output_csv * "_tile$tile", result_i, append = true)
        println("Finished tile $tile")
      end
    final_result = "Finished all tiles"
    else
      result = run_calculate_analogs(
          focal_data_cov,
          focal_data_mean,
          analog_data,
          var_names,
          n_analog_pool,
          n_analog_use,
          min_dist,
          max_dist,
          output_csv,
          error_file)
          CSV.write(output_csv, result, append = true)
          final_result = "Finished all points"
    end

    return final_result
  end
  
