using ProgressMeter
using RCall


 project_dir = "/project/umontana_climate_analogs/climate_analogs"
# project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
begin
  using Pkg; Pkg.activate(project_dir)
  Pkg.instantiate(); Pkg.precompile()
  end





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
end





using AlertPushover
using Alert

max_dist = 500 # in KM


annuals_futures_path = joinpath(project_dir, "data", "sensitivity", "annuals_futures_sens.Rds")
normal_futures_path = joinpath(project_dir, "data", "sensitivity", "normal_futures_sens.Rds")
analog_pool_path = joinpath(project_dir, "data", "sensitivity", "analog_pool_sens.Rds")
@rput annuals_futures_path
@rput normal_futures_path
@rput analog_pool_path



annuals_futures = rcopy(R"readRDS(annuals_futures_path)");



# grab historical focal cells
normal_futures = rcopy(R"readRDS(normal_futures_path)");




# analogs
analog_pool = rcopy(R"readRDS(analog_pool_path)");



try 


alert("datasets for 500x0x10p loaded")

catch e 

    println("loaded data") 

end 


#covert data types to lightest possible
for i in 1:length(annuals_futures)
  ## convert x and y to float16
  annuals_futures[i][!, :x] = Float32.(annuals_futures[i][!, :x])
  annuals_futures[i][!, :y] = Float32.(annuals_futures[i][!, :y])
  ## convert all other columns to INT16
  annuals_futures[i][!, Not(:x, :y)] = Int16.(annuals_futures[i][!, Not(:x, :y)])
end;






#covert data types to lightest possible
normal_futures[!, :x] = Float32.(normal_futures[!, :x]);
normal_futures[!, :y] = Float32.(normal_futures[!, :y]);
## convert all other columns to INT16 after rounding
normal_futures[!, Not(:x, :y)] = Int16.(round.(normal_futures[!, Not(:x, :y)], digits = 0));





analog_pool[!, :x] = Float32.(analog_pool[!, :x]);
analog_pool[!, :y] = Float32.(analog_pool[!, :y]);
analog_pool[!, Not(:x, :y)] = Int16.(round.(analog_pool[!, Not(:x, :y)], digits = 0));



# This takes 4.5 hours to ~1.4m analogs across washington points on 25 cores
# calculate n_analog_use


proportion_landscape = 0.10
n_analog_pool::Int32 = round(((2 * max_dist) / 0.270)^2 * proportion_landscape, digits = 0) |> Integer 

try 


alert("starting analogs 500x0x10p")


catch e 

    println("starting analog calculation") 

end

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
  (project_dir * "/data/sensitivity/500km_0min_10p"), ClimateAnalogs.calculate_analogs_threaded); #

alert("analog results 500x0x10 complete")
