using Printf
function format_value(prop::Float64)
    if prop < 1
        # Convert to string with 3 decimal places, remove the decimal point, pad with zeros
        return lpad(Int(1000 * prop), 3, "0")
    else
        # Convert to integer and format as a string
        return string(Int(100 * prop))
    end
end
# Define the parameters
distances = [ 500]
proportions = [0.001, 0.025, 0.05]
for distance in distances
    for prop in proportions
    prop_name = format_value(prop)
    project_dir = "/project/umontana_climate_analogs/climate_analogs"
    #project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
    # Define a template for the script
    script_content = """
    using ProgressMeter
    using RCall
    using Distributed
    println(nprocs())

    @everywhere  project_dir = "/project/umontana_climate_analogs/climate_analogs"
    # project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
    @everywhere begin
      using Pkg; Pkg.activate(project_dir)
      Pkg.instantiate(); Pkg.precompile()
      end\n\n\n


    @everywhere begin  
    using DataFrames
    using DataFramesMeta
    using ProgressMeter
    using Suppressor
    using Distributions
    using Distances
    using Random
    using CSV
    using Statistics
    using StatsBase
    using LinearAlgebra
    using ThreadPools: @bthreads
    #load in custom Julia functions
    include(joinpath(project_dir, "code", "src", "climate_analogs_jl", "climate_analogs.jl"))
    end\n\n



    using AlertPushover
    using Alert


    max_dist = $distance \n
    \n
    annuals_contemporary_path = joinpath(project_dir, "data", "washington", "annuals_contemporary_washington.Rds")
    normal_contemporary_path = joinpath(project_dir, "data", "washington", "normal_contemporary_washington.Rds")
    analog_pool_path = joinpath(project_dir, "data", "washington", "analog_pool_washington.Rds")
    @rput annuals_contemporary_path
    @rput normal_contemporary_path
    @rput analog_pool_path\n\n

    annuals_contemporary = rcopy(R"readRDS(annuals_contemporary_path)");\n\n

    # grab historical focal cells
    normal_contemporary = rcopy(R"readRDS(normal_contemporary_path)");\n\n


    # analogs
    analog_pool = rcopy(R"readRDS(analog_pool_path)");\n\n

    try 

    alert("datasets for $(distance)x0x$(prop_name)p loaded")
    \n
    catch e \n
        println("loaded data") \n
    end \n

    #covert data types to lightest possible
    for i in 1:length(annuals_contemporary)
      ## convert x and y to float16
      annuals_contemporary[i][!, :x] = Float32.(annuals_contemporary[i][!, :x])
      annuals_contemporary[i][!, :y] = Float32.(annuals_contemporary[i][!, :y])
      ## convert all other columns to INT16
      annuals_contemporary[i][!, Not(:x, :y)] = Int16.(annuals_contemporary[i][!, Not(:x, :y)])
    end;\n\n\n



    #covert data types to lightest possible
    normal_contemporary[!, :x] = Float32.(normal_contemporary[!, :x]);
    normal_contemporary[!, :y] = Float32.(normal_contemporary[!, :y]);
    ## convert all other columns to INT16 after rounding
    normal_contemporary[!, Not(:x, :y)] = Int16.(round.(normal_contemporary[!, Not(:x, :y)], digits = 0));\n\n\n


    analog_pool[!, :x] = Float32.(analog_pool[!, :x]);
    analog_pool[!, :y] = Float32.(analog_pool[!, :y]);
    analog_pool[!, Not(:x, :y)] = Int16.(round.(analog_pool[!, Not(:x, :y)], digits = 0));\n\n

    # This takes 4.5 hours to ~1.4m analogs across washington points on 25 cores
    # calculate n_analog_use
    proportion_landscape = $prop
    n_analog_pool::Int32 = round(((2 * max_dist) / 0.270)^2 * proportion_landscape, digits = 0) |> Integer \n
    try \n
    alert("starting analogs $(distance)x0x$(prop_name)p")
    \n
    catch e \n
        println("starting analog calculation") \n
    end\n
    # library(profvis)
    #start = time()
    # profvis({
      # test first 12
      
      #  for i in 1:size(annuals_contemporary,1)
      #    annuals_contemporary[i] = annuals_contemporary[i][1:12, :]
      #  end
      
      analog_results = ClimateAnalogs.find_analogs( annuals_contemporary,
      normal_contemporary,
      analog_pool,
      ["aet", "def", "tmax", "tmin"],
      n_analog_pool,
      100,
      50, # In KM! 
      max_dist, # In KM!
      (project_dir * "/data/washington/$(distance)km_0min_$(prop_name)p"));\n 
        alert("analog results $(distance)x0x$(prop_name) complete")  

    """
             # Write the script to a file
            filename = project_dir * "/code/washington/$(distance)_0m_$(prop_name)p.jl"
            open(filename, "w") do f
                write(f, script_content)
            end
        end
end
