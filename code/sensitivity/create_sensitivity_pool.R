setwd("/project/umontana_climate_analogs/climate_analogs")
library(terra)
library(data.table)
library(tidyverse)
source(file.path(getwd(), "code/src/climate_analogs_R/as.data.table.R"))

climate_dir <- file.path(getwd(), "data/climate")
#create template
template <- rast(file.path(getwd(),"data/climate/topoterra_hist_1961-2022.tif"),
                 lyrs = 1)
classify_mat <- matrix(c(NA, NA), ncol = 2)
template_binary <- classify(template, classify_mat, others = 1)

total_cells <- global(template_binary, "sum", na.rm = T) |> as.integer()
sample_size <- round(0.10 * total_cells)


template_sample <- template_binary
v <- cells(template_sample)
v_sample <- sample(v, size = sample_size, replace = FALSE)
v_antisample <- setdiff(v, v_sample)
# check that length(v_sample) = length(v)
length(v_sample) + length(v_antisample) == length(v)

template_sample[v_antisample] <- NA
plot(c(template_binary, template_sample))
# check that sum of template_sample = sample_size
global(template_sample, "sum", na.rm = TRUE) |> as.integer() == sample_size

ann_futures <- grep("topoterra_2C_[0-9]{4}.tif",
                list.files(climate_dir),
                value = TRUE) |>
  paste0(climate_dir, "/", ... = _) |>
  map(\(x) {
    rast(x) * template_sample
  }
  )

# Convert each element into a data.table using a for loop for memory efficiency
annuals_futures <- vector("list", length(ann_futures))
for (i in seq_along(ann_futures)) {
  annuals_futures[[i]] <- as.data.table(ann_futures[[i]], xy = TRUE)
  annuals_futures[[i]] <- na.omit(annuals_futures[[i]])
  ann_futures[[i]] <- NULL
  gc()
}
# # # This is slow, save for convenience for now


normal_futures <- list.files(
  climate_dir, "topoterra_2C_[0-9]{4}-[0-9]{4}.tif",
  full.names = TRUE
) |>
  rast() %>%
  "*"(., template_sample) |>
  as.data.table(xy = TRUE) |>
  na.omit()

analog_pool <- list.files(
                          climate_dir, "topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
  full.names = TRUE) |>
  rast() |>
  as.data.table(xy = TRUE) |>
  na.omit()
saveRDS(analog_pool, file.path(getwd(), "data/sensitivity/analog_pool_sens.Rds"))
  #verify that x and y coordinates of normal_futures and annuals_futures are the same
  #map it to all annual futures
match_check <- map_lgl(annuals_futures, \(x) {
    future_i <- x[, c("x", "y")]
    normal_futures_i <- normal_futures[, c("x", "y")]
    all.equal(future_i, normal_futures_i)

    }
    ) |> sum() / length(annuals_futures)
if(match_check == 1){
  print("All x and y coordinates match")
saveRDS(normal_futures,
        file.path(getwd(), "data/sensitivity/normal_futures_sens.Rds" ))
saveRDS(annuals_futures,
        file.path(getwd(), "data/sensitivity/annuals_futures_sens.Rds"))
} else {
  print("Some x and y coordinates do not match")
}



