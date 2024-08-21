using DataFrames
using DataFramesMeta
using ProgressMeter
using Suppressor
using CSV
using Base.Threads

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
    analog_data = nothing
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
  
  function calculate_analogs_distributed(
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

  function calculate_analogs_threaded(
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

    # Determine the total number of rows needed in the output DataFrame
    total_rows = nrow(focal_data_cov[1]) * n_analog_use

    # Preallocate the output DataFrame with the required columns
    result_df = DataFrame(
      a_x = Vector{Float32}(undef,total_rows),
      a_y = Vector{Float32}(undef,total_rows),
      md = Vector{Float32}(undef,total_rows),
      dist_km = Vector{Float32}(undef,total_rows),
      sigma = Vector{Float32}(undef,total_rows),
      f_x = Vector{Float32}(undef,total_rows),
      f_y = Vector{Float32}(undef,total_rows)
    
    )

    # Initialize a lock for thread-safe operations
    lk = ReentrantLock()

    

 # Determine the number of threads and split data accordingly
 n_threads = nthreads()
 n_per_thread = ceil(Int, size(focal_data_cov[1], 1) / n_threads)
 p = Progress(nrow(focal_data_cov[1]))

 result_df = Threads.@threads for t in 1:n_threads
     start_idx = (t - 1) * n_per_thread + 1
     end_idx = min(t * n_per_thread, size(focal_data_cov[1], 1))
     focal_data_cov_i = Vector{Any}(undef, size(focal_data_cov,1))
     for i in 1:size(focal_data_cov,1)
      focal_data_cov_i[i] = focal_data_cov[i][start_idx:end_idx, :]
     end
      focal_data_mean_i = focal_data_mean[start_idx:end_idx, :]
    
     #preallocate local DataFrame
     i_rows = (end_idx-start_idx+1) * n_analog_use
    df_i = DataFrame(
      a_x = Vector{Float32}(undef,i_rows),
      a_y = Vector{Float32}(undef,i_rows),
      md = Vector{Float32}(undef,i_rows),
      dist_km = Vector{Float32}(undef,i_rows),
      sigma = Vector{Float32}(undef,i_rows),
      f_x = Vector{Float32}(undef,i_rows),
      f_y = Vector{Float32}(undef,i_rows))

     for x in start_idx:end_idx
        adj_x = x-(start_idx-1)
        df_x_1 = ((adj_x-1)*100) + 1
        df_x_end = df_x_1+100-1
        
         try
             # Calculate analogs for the current focal point
             df_i[df_x_1:df_x_end, :] = calculate_analogs(
                 focal_data_cov_i,
                 focal_data_mean_i,
                 analog_data,
                 var_names,
                 n_analog_pool,
                 n_analog_use,
                 min_dist,
                 max_dist,
                 adj_x
             )
             next!(p)
            GC.gc(false)
             # Check memory and log errors, thread-safely
             if x % 100 == 0 && Sys.free_memory() < 10e9
                 lock(lk) do
                     open(error_file, "a") do f
                         write(f, "Broke at $x\n")
                     end
                 end
             end
         catch e
             lock(lk) do
                 open(error_file, "a") do f
                     write(f, "$x\n")
                 end
             end
         end
     end

 end

  # Return the filled DataFrame
  result_df = reduce(vcat, result_df)
  return result_df
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
    output_file::String,
    compute_function::Function)
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
        result_i = compute_function(
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
      result = compute_function(
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
  
