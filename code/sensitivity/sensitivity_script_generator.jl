using Printf

# Define the parameters
distances = [250, 500, 1000]
proportions = [0.1, 0.2, 0.5]
project_dir = "/project/umontana_climate_analogs/climate_analogs"
#project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
# Define a template for the script
script_template1 = """
using ProgressMeter
using RCall


 project_dir = "/project/umontana_climate_analogs/climate_analogs"
# project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
begin
  using Pkg; Pkg.activate(project_dir)
  Pkg.instantiate(); Pkg.precompile()
  end\n\n\n


begin  
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
#load in custom Julia functions
include(joinpath(project_dir, "code", "src", "climate_analogs_jl", "climate_analogs.jl"))
end\n\n



using AlertPushover
using Alert
"""


script_template2 = """
\n
annuals_futures_path = joinpath(project_dir, "data", "sensitivity", "annuals_futures_sens.Rds")
normal_futures_path = joinpath(project_dir, "data", "sensitivity", "normal_futures_sens.Rds")
analog_pool_path = joinpath(project_dir, "data", "sensitivity", "analog_pool_sens.Rds")
@rput annuals_futures_path
@rput normal_futures_path
@rput analog_pool_path\n\n

annuals_futures = rcopy(R"readRDS(annuals_futures_path)");\n\n

# grab historical focal cells
normal_futures = rcopy(R"readRDS(normal_futures_path)");\n\n


# analogs
analog_pool = rcopy(R"readRDS(analog_pool_path)");\n\n

try 

"""

script_template3 = """
catch e \n
    println("loaded data") \n
end \n

#covert data types to lightest possible
for i in 1:length(annuals_futures)
  ## convert x and y to float16
  annuals_futures[i][!, :x] = Float32.(annuals_futures[i][!, :x])
  annuals_futures[i][!, :y] = Float32.(annuals_futures[i][!, :y])
  ## convert all other columns to INT16
  annuals_futures[i][!, Not(:x, :y)] = Int16.(annuals_futures[i][!, Not(:x, :y)])
end;\n\n\n



#covert data types to lightest possible
normal_futures[!, :x] = Float32.(normal_futures[!, :x]);
normal_futures[!, :y] = Float32.(normal_futures[!, :y]);
## convert all other columns to INT16 after rounding
normal_futures[!, Not(:x, :y)] = Int16.(round.(normal_futures[!, Not(:x, :y)], digits = 0));\n\n\n


analog_pool[!, :x] = Float32.(analog_pool[!, :x]);
analog_pool[!, :y] = Float32.(analog_pool[!, :y]);
analog_pool[!, Not(:x, :y)] = Int16.(round.(analog_pool[!, Not(:x, :y)], digits = 0));\n\n

# This takes 4.5 hours to ~1.4m analogs across washington points on 25 cores
# calculate n_analog_use
"""
script_template4 = """
n_analog_pool::Int32 = round(((2 * max_dist) / 0.270)^2 * proportion_landscape, digits = 0) |> Integer \n
try \n
"""
script_template5 = """
catch e \n
    println("starting analog calculation") \n
end\n
# library(profvis)
#start = time()
# profvis({
  # test first 12
  
  #  for i in 1:size(annuals_futures,1)
  #    annuals_futures[i] = annuals_futures[i][1:12, :]
  #  end
  
  analog_results = ClimateAnalogs.find_analogs( annuals_futures,
  normal_futures,
  analog_pool,
  ["aet", "def", "tmax", "tmin"],
  n_analog_pool,
  100,
  0, # In KM! 
  max_dist, # In KM!
  """
  
  
  
  # Iterate over all combinations of distances and proportions
  for distance in distances
    for proportion in proportions
        fstring1 = @sprintf("max_dist = %d # in KM", distance)
        fstring2 = @sprintf("""alert("datasets for %dx0x%dp loaded")""", distance, Int(proportion * 100))
        fstring3 = @sprintf("proportion_landscape = %.2f", proportion)
        fstring4 = @sprintf("""alert("starting analogs %dx0x%dp")""", distance, Int(proportion * 100))
        fstring5 = @sprintf("""(project_dir * "/data/sensitivity/%dkm_0min_%dp"), ClimateAnalogs.calculate_analogs_threaded); #""", distance, Int(proportion * 100))
        fstring6 = @sprintf("""alert("analog results %dx0x%d complete")""", distance, Int(proportion * 100))

        script_content = 
          script_template1 * "\n" * 
          fstring1 * "\n" * 
          script_template2 * "\n" * 
          fstring2 * "\n\n" * 
          script_template3 * "\n\n" * 
          fstring3 * "\n" * 
          script_template4 * "\n" * 
          fstring4 * "\n\n\n" *
          script_template5 * 
          fstring5 * "\n\n" * 
          fstring6 * "\n"
        
        # Write the script to a file
        filename = project_dir * "/code/sensitivity/" *
         @sprintf("%d_0m_%dp.jl", distance, Int(proportion * 100))
        open(filename, "w") do f
            write(f, script_content)
        end
    end
end
