# developing a simple sampling procedure that represents the climate 
# space of the western US well with TopoTerra.
library(terra)
library(tidyverse)
library(data.table)
source("code/as.data.table.R")
# Read in climate data
climate_normal <- list.files(
    "data/climate",
    "topoterra_hist_[0-9]{4}-[0-9]{4}.tif",
    full.names = TRUE
) |>
    rast()

# Convert to data.table
climate_dt <- climate_normal |>
    as.data.table(xy = TRUE, na.rm = TRUE)

# create and plot climate envelope with PCA
climate_pca <- climate_dt |>
    select(-c("x", "y")) |>
    as.matrix() |>
    prcomp(center = TRUE, scale. = TRUE)
# plot the first two PCs as contour using ggplot
climate_pca <- data.table(
  PC1 = climate_pca$x[, 1],
  PC2 = climate_pca$x[, 2]
) |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_density_2d() +
  theme_minimal() +
  labs(title = "Climate envelope of TopoTerra",
       x = "PC1",
       y = "PC2")
print(climate_pca)
