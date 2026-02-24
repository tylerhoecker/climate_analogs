# setwd("/project/umontana_climate_analogs/climate_analogs")
library(terra)
library(data.table)
library(dplyr)
library(purrr)
source(file.path(getwd(), "code/src/climate_analogs_R/as.data.table.R"))

###### Create sample pixels and number them by tile for WNA validation ######


climate_dir <- file.path(getwd(), "data/climate")
# create template
template <- rast(file.path(getwd(), "data/climate/TopoTerra_1961-2022.tif"),
  lyrs = 1
)

# sample test

# sample pixels
focal_mask <- rast("data/lf_mask_focal.tif")
analog_mask <- rast("data/lf_mask_analog.tif")
sample_pixels <- function(template, sample_size) {
  validcells <- cells(template)
  sample_indices <- sample(validcells, size = sample_size, replace = FALSE)
  template[sample_indices] <- 1
  template[-sample_indices] <- NA
  return(template)
}

sample_size <- 100000
template_sample <- sample_pixels(
  template,
  sample_size
)

template_focal <- template * focal_mask
template_analog <- template * analog_mask
# build extents
col_cells <- ceiling(dim(template)[2] / 4)
row_cells <- ceiling(dim(template)[1] / 5)

extents <- getTileExtents(template, c(row_cells, col_cells))


# number the sampled pixels by which tile they belong to
numbered_template_sample <- template_sample %>%
  classify(matrix(c(NA, NA), ncol = 2), others = 0)
c <- 0
for (i in seq_len(nrow(extents))) {
  extentt <- extents[i, ] |> ext()
  # Crop the template_sample to the current extent
  cropped_sample <- crop(template_sample, extentt)
  # Get the indices of non-NA cells in the cropped sample
  valid_cells <- which(!is.na(values(cropped_sample)))
  l <- length(valid_cells)
  c <- c + l
  print(paste0(
    "Numbering tile ", i, " with ", l, " cells. ",
    "Cumulative is ", c, " cells."
  ))
  # Map the valid cells back to the original template_sample
  global_indices <- cellFromXY(template_sample, xyFromCell(cropped_sample, valid_cells))
  # Only update cells that are currently NA
  numbered_template_sample[global_indices] <- i
}
# now make it data table and save RDS
numbered_template_sample_dt <- as.data.table(numbered_template_sample, xy = TRUE) |>
  na.omit()
saveRDS(numbered_template_sample_dt, file.path(getwd(), "data/validation/numbered_template_sample.Rds"))


ann_contemp <- paste0(climate_dir, "/TopoTerra_", 1985:2015, ".tif") |>
  map(\(x) {
    rast(x) * template_sample
  })
annuals_contemporary <- vector("list", length(ann_contemp))
for (i in seq_along(ann_contemp)) {
  annuals_contemporary[[i]] <- as.data.table(ann_contemp[[i]], xy = TRUE)
  annuals_contemporary[[i]] <- na.omit(annuals_contemporary[[i]])
  ann_contemp[[i]] <- TRUE
  gc()
}
saveRDS(annuals_contemporary, file.path(getwd(), "data/validation/annuals_contemporary_sample.Rds"))


normal_contemporary <- list.files(
  climate_dir, "TopoTerra_1985-2015.tif",
  full.names = TRUE
) |>
  rast() %>%
  "*"(., template_sample) |>
  as.data.table(xy = TRUE) |>
  na.omit()
saveRDS(
  normal_contemporary,
  file.path(getwd(), "data/validation/normal_contemporary_sample.Rds")
)


# build 20 analog pools so they can be loaded as needed
analog_pool <- list.files(
  climate_dir, "TopoTerra_1985-2015.tif",
  full.names = TRUE
) |>
  rast()
for (i in 1:20) {
  extentt <- extents[i, ] |>
    ext() |>
    vect(crs = crs("EPSG:4326")) |>
    buffer(500 * 1000)
  analog_pool_t <- analog_pool %>%
    crop(., extentt) %>%
    as.data.table(xy = TRUE) |>
    na.omit()
  saveRDS(analog_pool_t, file.path(getwd(), paste0("data/validation/analog_pool_sample_", i, ".Rds")))
}
