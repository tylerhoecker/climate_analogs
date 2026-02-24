using ProgressMeter
using RCall
using Distributed
using CodecZlib
addprocs(30)

# PRIMARY VALIDATION SCRIPT FOR EXCLUSION RADII
println("Number of processes: " * string(nprocs()))

# @everywhere project_dir = "/project/umontana_climate_analogs/climate_analogs"
@everywhere project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
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
tile_number = joinpath(project_dir, "data", "validation", "numbered_template_sample.Rds")
annuals_contemporary_path = joinpath(project_dir, "data", "validation", "annuals_contemporary_sample.Rds")
normal_contemporary_path = joinpath(project_dir, "data", "validation", "normal_contemporary_sample.Rds")
@rput tile_number
@rput annuals_contemporary_path
@rput normal_contemporary_path
annuals_contemporary = rcopy(R"readRDS(annuals_contemporary_path)");
normal_contemporary = rcopy(R"readRDS(normal_contemporary_path)");
tile_ids = rcopy(R"readRDS(tile_number)")
rename!(tile_ids, Dict(:aet => :tile_id, :x => :x, :y => :y))
for i in 1:length(annuals_contemporary)
  ## convert x and y to float16
  annuals_contemporary[i][!, :x] = Float32.(annuals_contemporary[i][!, :x])
  annuals_contemporary[i][!, :y] = Float32.(annuals_contemporary[i][!, :y])
  ## convert all other columns to INT16
  annuals_contemporary[i][!, Not(:x, :y)] = Int16.(round.(annuals_contemporary[i][!, Not(:x, :y)], digits=0))
end;


#covert data types to lightest possible
normal_contemporary[!, :x] = Float32.(normal_contemporary[!, :x]);
normal_contemporary[!, :y] = Float32.(normal_contemporary[!, :y]);
## convert all other columns to INT16 after rounding
normal_contemporary[!, Not(:x, :y)] = Int16.(round.(normal_contemporary[!, Not(:x, :y)], digits=0));


# calculate n_analog_use
proportion_landscape = 0.05
n_analog_pool::Int32 = round(((2 * max_dist) / 0.270)^2 * proportion_landscape, digits=0) |> Integer

for i in 1:20
  tile = i
  tile_inds = findall(row -> row[3] == tile, eachrow(tile_ids))
  if length(tile_inds) == 0
    println("no points for tile $(tile), skipping")
    continue
  end
  annuals_contemporaryᵢ = [annual[tile_inds, :] for annual in annuals_contemporary]
  normal_contemporaryᵢ = normal_contemporary[tile_inds, :]
  analog_pool_path = joinpath(project_dir, "data", "validation", "analog_pool_sample_$(tile).Rds")
  @rput analog_pool_path
  try
    analog_pool = rcopy(R"readRDS(analog_pool_path)")
  catch e
    println("could not read analog pool for tile $(tile), likely does not exist")
    continue
  end
  analog_pool[!, :x] = Float32.(analog_pool[!, :x])
  analog_pool[!, :y] = Float32.(analog_pool[!, :y])
  analog_pool[!, Not(:x, :y)] = Int16.(round.(analog_pool[!, Not(:x, :y)], digits=0))
  # This takes 4.5 hours to ~1.4m analogs across reverse_analogs points on 25 cores
  #   for i in 1:size(annuals_contemporaryᵢ,1)
  #     annuals_contemporaryᵢ[i] = annuals_contemporaryᵢ[i][1:12, :]
  #   end

  t1 = time()
  ClimateAnalogs.find_analogs(annuals_contemporaryᵢ,
    normal_contemporaryᵢ,
    analog_pool,
    ["aet", "def", "tmax", "tmin"],
    n_analog_pool,# number of analog pixels to sample
    100, # number of analogs to keep
    0, # In KM!
    max_dist, # In KM!
    (project_dir * "/data/validation/outputs"), "validation_0km_exclusion_$(tile)")
  t2 = time()
  println("Time taken for tile $(tile): $((t2 - t1)/60) minutes")


  t1 = time()
  ClimateAnalogs.find_analogs(annuals_contemporaryᵢ,
    normal_contemporaryᵢ,
    analog_pool,
    ["aet", "def", "tmax", "tmin"],
    n_analog_pool,# number of analog pixels to sample
    100, # number of analogs to keep
    25, # In KM!
    max_dist, # In KM!
    (project_dir * "/data/validation/outputs"), "validation_25km_exclusion_$(tile)")
  t2 = time()
  println("Time taken for tile $(tile): $((t2 - t1)/60) minutes")

end
