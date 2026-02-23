#####################
# PRIMARY SCRIPT FOR CONTEMPORARY AND FUTURE REVERSE ANALOGS. PATHS AND INPUTS CHANGED AS NEEDED
######################

library(tidyverse)
library(data.table)
library(doParallel)
library(future)

cpus_per_task <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
registerDoParallel(cpus_per_task)

paste0("Number of processes: ", as.character(cpus_per_task))

project_dir <- "INSERT/PATH"
# project_dir = "/home/jeff/Github/ClimateAnalogues/climate_analogs"
# Load in custom R functions
source(file.path(project_dir, "code", "src", "climate_analogs_R", "climate_analogs.R"))




max_dist <- 500
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
    stop("Stopping: Args must be supplied")
} else if (length(commandArgs(trailingOnly = TRUE)) > 1) {
    stop("Stopping: Only 1 arg can be supplied")
} else {
    tile_id <- commandArgs(trailingOnly = TRUE)[1]
}

# Path to rds files for the tile
annuals_future_path <- file.path(project_dir, "data", "reverse_analogs", "inputs", paste0("annuals_future_", tile_id, ".Rds"))
normal_future_path <- file.path(project_dir, "data", "reverse_analogs", "inputs", paste0("normal_future_", tile_id, ".Rds"))
analog_pool_path <- file.path(project_dir, "data", "reverse_analogs", "inputs", paste0("analog_pool_", tile_id, ".Rds"))
annuals_future <- readRDS(annuals_future_path)
# Grab historical focal cells
normal_future <- readRDS(normal_future_path)
# Analogs
analog_pool <- readRDS(analog_pool_path)


try(
    {
        message("datasets for tile ", tile_id, " loaded")
    },
    silent = TRUE
)

message("loaded data")



# Calculate n_analog_use which is the sample size to draw from the analog pool for each focal cell.
proportion_landscape <- 0.05
n_analog_pool <- round(((2 * max_dist) / 0.270)^2 * proportion_landscape, digits = 0) |> as.numeric()

tryCatch(
    {
        message("starting analogs tile ", tile_id)
    },
    error = function(e) {
        message("starting analog calculation")
    }
)

analog_results <- find_analogs(
    annuals_future,
    normal_future,
    analog_pool,
    c("aet", "def", "tmax", "tmin"),
    n_analog_pool, # number of analog pixels to sample
    100, # number of analogs to keep
    0, # In KM!
    max_dist, # In KM!
    paste0(project_dir, "/data/reverse_analogs/outputs/reverse_analogs_", tile_id)
)
