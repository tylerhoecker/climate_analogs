# setwd("/project/umontana_climate_analogs/climate_analogs")
library(terra)
library(data.table)
library(dplyr)
library(purrr)
source(file.path(getwd(), "code/src/climate_analogs_R/as.data.table.R"))

#


climate_dir <- file.path(getwd(), "data/climate")
# create template
template <- rast(file.path(getwd(), "data/climate/topoterra_hist_1961-2022.tif"),
  lyrs = 1
)
lf_mask <- rast("data/lf_water_mask_4326.tiff") %>%
  crop(., template) %>%
  classify(matrix(c(-1111, NA), ncol = 2), others = 1)

template <- template * lf_mask




states <- vect("data/western_states/western_states.shp")
washington <- subset(states, NAME == "Washington", NSE = TRUE) %>%
  project(crs(template))
classify_mat <- matrix(c(NA, NA), ncol = 2)
template_binary <- classify(template, classify_mat, others = 1) %>%
  mask(washington)
template_buffer <- classify(template, classify_mat, others = 1) %>%
  mask(buffer(washington, 500 * 1000))

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
years <- 1961:1990
ann_contemp <- paste0(climate_dir, "/topoterra_hist_", years, ".tif") |>
  map(\(x) {
    rast(x) * template_sample
  })
annuals_contemporary <- vector("list", length(ann_contemp))
for (i in seq_along(ann_contemp)) {
  annuals_contemporary[[i]] <- as.data.table(ann_contemp[[i]], xy = TRUE)
  annuals_contemporary[[i]] <- na.omit(annuals_contemporary[[i]])
  ann_contemp[[i]] <- T
  gc()
}
# # # This is slow, save for convenience for now


normal_contemporary <- list.files(
  climate_dir, "topoterra_hist_1961-1990.tif",
  full.names = TRUE
) |>
  rast() %>%
  "*"(., template_sample) |>
  as.data.table(xy = TRUE) |>
  na.omit()

analog_pool <- list.files(
  climate_dir, "topoterra_hist_1961-1990.tif",
  full.names = TRUE
) |>
  rast() %>%
  "*"(., template_buffer) %>%
  as.data.table(xy = TRUE) |>
  na.omit()
saveRDS(analog_pool, file.path(getwd(), "data/washington/analog_pool_washington.Rds"))
# verify that x and y coordinates of normal_contemporary and annuals_contemporary are the same
# map it to all annual contemporary
match_check <- map_lgl(annuals_contemporary, \(x) {
  future_i <- x[, c("x", "y")]
  normal_contemporary_i <- normal_contemporary[, c("x", "y")]
  all.equal(future_i, normal_contemporary_i)
}) |> sum() / length(annuals_contemporary)
if (match_check == 1) {
  print("All x and y coordinates match")
  saveRDS(
    normal_contemporary,
    file.path(getwd(), "data/washington/normal_contemporary_washington.Rds")
  )
  saveRDS(
    annuals_contemporary,
    file.path(getwd(), "data/washington/annuals_contemporary_washington.Rds")
  )
} else {
  print("Some x and y coordinates do not match")
}
