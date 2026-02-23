setwd("/project/umontana_climate_analogs/climate_analogs")
library(terra)
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(foreach)
source(file.path(getwd(), "code/src/climate_analogs_R/as.data.table.R"))

##############
# PRIMARY SCRIPT TO CREATE CONTEMPORARY AND FUTURE REVERSE ANALOG DATASETS FOR ANALYSIS
# PATHS ARE CHANGED AS NEEDED TO REFLECT REVERSE ANALOGS OR CONTEMPORARY ANALOGS
##############


climate_dir <- file.path(getwd(), "data/climate")
# create template
template <- rast(file.path(getwd(), "data/climate/topoterra_hist_1961-2022.tif"),
  lyrs = 1
)
focal_mask <- rast("data/lf_mask_focal.tif")
analog_mask <- rast("data/lf_mask_analog.tif")

template_focal <- template * focal_mask
template_analog <- template * analog_mask
# build extents
col_cells <- ceiling(dim(template)[2] / 4)
row_cells <- ceiling(dim(template)[1] / 5)

extents <- getTileExtents(template, c(row_cells, col_cells))





args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("tile number must be supplied as an args")
} else if (length(args) > 1) {
  stop("Only accepts one argument.")
} else {
  tile_id <- args[1]
}


local_ext <- ext(extents[as.numeric(tile_id), ]) %>%
  vect(crs = crs(template))
classify_mat <- matrix(c(NA, NA), ncol = 2)
template_binary <- classify(template_focal, classify_mat, others = 1) %>%
  mask(local_ext)
template_buffer <- classify(template_analog, classify_mat, others = 1) %>%
  mask(buffer(local_ext, 500 * 1000))

total_cells <- global(template_binary, "sum", na.rm = T) |> as.integer()
sample_size <- round(1 * total_cells)


template_sample <- template_binary


# check that sum of template_sample = sample_size
global(template_sample, "sum", na.rm = TRUE) |> as.integer() == sample_size
years <- 1985:2015
ann_fut <- paste0(climate_dir, "/topoterra_2C_", years, ".tif") |>
  map(\(x) {
    rast(x) * template_sample
  })
annuals_future <- vector("list", length(ann_fut))
for (i in seq_along(ann_fut)) {
  annuals_future[[i]] <- as.data.table(ann_fut[[i]], xy = TRUE)
  annuals_future[[i]] <- na.omit(annuals_future[[i]])
  ann_fut[[i]] <- T
  gc()
}
# # # This is slow, save for convenience for now


normal_future <- list.files(
  climate_dir, "topoterra_2C_1985-2015.tif",
  full.names = TRUE
) |>
  rast() %>%
  "*"(., template_sample) |>
  as.data.table(xy = TRUE) |>
  na.omit()

analog_pool <- list.files(
  climate_dir, "topoterra_hist_1985-2015.tif",
  full.names = TRUE
) |>
  rast() %>%
  "*"(., template_buffer) %>%
  as.data.table(xy = TRUE) |>
  na.omit()
saveRDS(analog_pool, file.path(getwd(), paste0("data/reverse_analogs/inputs/analog_pool_", tile_id, ".Rds")))
# verify that x and y coordinates of normal_future and annuals_future are the same
# map it to all annual contemporary
match_check <- map_lgl(annuals_future, \(x) {
  future_i <- x[, c("x", "y")]
  normal_future_i <- normal_future[, c("x", "y")]
  all.equal(future_i, normal_future_i)
}) |> sum() / length(annuals_future)
if (match_check == 1) {
  print("All x and y coordinates match")
  saveRDS(
    normal_future,
    file.path(getwd(), paste0("data/reverse_analogs/inputs/normal_future_", tile_id, ".Rds"))
  )
  saveRDS(
    annuals_future,
    file.path(getwd(), paste0("data/reverse_analogs/inputs/annuals_future_", tile_id, ".Rds"))
  )
} else {
  print("Some x and y coordinates do not match")
}
