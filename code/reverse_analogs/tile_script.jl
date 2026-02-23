using ProgressMeter
using RCall
using Distributed
using CodecZlib
cpus_per_task = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
addprocs(cpus_per_task)

#####################
# PRIMARY SCRIPT FOR CONTEMPORARY AND FUTURE REVERSE ANALOGS. PATHS AND INPUTS CHANGED AS NEEDED
######################
println("Number of processes: " * string(nprocs()))

@everywhere project_dir = "/project/umontana_climate_analogs/climate_analogs"
# project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
@everywhere begin
  using Pkg
  Pkg.activate(project_dir)
  Pkg.instantiate()
  Pkg.precompile()
end


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
end



using AlertPushover
using Alert

max_dist = 500
if length(ARGS) == 0
  error("Stopping: Args must be supplied")
elseif length(ARGS) > 1
  error("Stopping: Only 1 arg can be supplied")
else
  tile_id = ARGS[1]
end



# path to rds files for the tile
annuals_future_path = joinpath(project_dir, "data", "reverse_analogs", "inputs", "annuals_future_$(tile_id).Rds")
normal_future_path = joinpath(project_dir, "data", "reverse_analogs", "inputs", "normal_future_$(tile_id).Rds")
analog_pool_path = joinpath(project_dir, "data", "reverse_analogs", "inputs", "analog_pool_$(tile_id).Rds")
@rput annuals_future_path
@rput normal_future_path
@rput analog_pool_path



annuals_future = rcopy(R"readRDS(annuals_future_path)");



# grab future focal normals
normal_future = rcopy(R"readRDS(normal_future_path)");




# analogs
analog_pool = rcopy(R"readRDS(analog_pool_path)");



try

  alert("datasets for tile $(tile_id) loaded")


catch e

  println("loaded data")

end


#covert data types to lightest reasonable
for i in 1:length(annuals_future)
  ## convert x and y to float16
  annuals_future[i][!, :x] = Float32.(annuals_future[i][!, :x])
  annuals_future[i][!, :y] = Float32.(annuals_future[i][!, :y])
  ## convert all other columns to INT16
  annuals_future[i][!, Not(:x, :y)] = Int16.(annuals_future[i][!, Not(:x, :y)])
end;






#covert data types to lightest reasonable
normal_future[!, :x] = Float32.(normal_future[!, :x]);
normal_future[!, :y] = Float32.(normal_future[!, :y]);
## convert all other columns to INT16 after rounding
normal_future[!, Not(:x, :y)] = Int16.(round.(normal_future[!, Not(:x, :y)], digits=0));


# covert data types to lightest reasonable
analog_pool[!, :x] = Float32.(analog_pool[!, :x]);
analog_pool[!, :y] = Float32.(analog_pool[!, :y]);
analog_pool[!, Not(:x, :y)] = Int16.(round.(analog_pool[!, Not(:x, :y)], digits=0));



# This takes 4.5 hours to ~1.4m analogs across reverse_analogs points on 25 cores
# calculate n_analog_use which is the sample size to draw from the analog pool for each focal cell.
proportion_landscape = 0.05
n_analog_pool::Int32 = round(((2 * max_dist) / 0.270)^2 * proportion_landscape, digits=0) |> Integer

try
  alert("starting analogs tile $(tile_id)")
catch e
  println("starting analog calculation")
end


# on a sample of 1.4 million analogs this takes 4.5 hours on 25 cores
analog_results = ClimateAnalogs.find_analogs(annuals_future,
  normal_future,
  analog_pool,
  ["aet", "def", "tmax", "tmin"],
  n_analog_pool,# number of analog pixels to sample
  100, # number of analogs to keep
  0, # In KM!
  max_dist, # In KM!
  (project_dir * "/data/reverse_analogs/outputs"),
  "reverse_analogs_$(tile)");

alert("analog results tile $(tile_id) complete")
