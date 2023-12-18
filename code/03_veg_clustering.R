library(motif)
library(stars)
library(tmap)
library(tidyverse)

# Read in climate analog info as a set of tiffs
current <- read_stars('../output/gee/wa_current_bps.tiff')
future <- read_stars('../output/gee/wa_bps1.tiff')
analog_data <- c(current, future)

# Define clusters
eco_signature = lsp_signature(analog_data,
                              type = "incove",
                              window = 50)

# Calculate distances among clusters
eco_dist = lsp_to_dist(eco_signature, dist_fun = "jensen-shannon")

# Heirarchical cluster
eco_hclust = hclust(eco_dist, method = "ward.D2")
plot(eco_hclust)

clusters = cutree(eco_hclust, k = 15)

eco_grid_sf = lsp_add_clusters(eco_signature,
                               clusters)


tm_clu = tm_shape(eco_grid_sf) +
  tm_polygons("clust", style = "cat", palette = "Set2", title = "Cluster:") +
  tm_layout(legend.position = c("LEFT", "BOTTOM"))
tm_clu

eco_grid_sf2 = eco_grid_sf %>%
  group_by(clust) %>%
  summarize()

tm_shape(analog_data) +
  tm_raster(style = "cat") +
  tm_facets(ncol = 2) +
  tm_shape(eco_grid_sf2) +
  tm_borders(col = "black") +
  tm_layout(legend.show = FALSE, 
            title.position = c("LEFT", "TOP"))

# Quality
eco_grid_sfq = lsp_add_quality(eco_grid_sf, eco_dist, type = "cluster")

eco_grid_sfq2 = eco_grid_sfq %>%
  group_by(clust) %>%
  summarise(inhomogeneity = mean(inhomogeneity),
            distinction = mean(distinction),
            quality = mean(quality))

st_write(eco_grid_sfq2, '../output/gee/clusters.shp')

tm_inh = tm_shape(eco_grid_sfq2) +
  tm_polygons("inhomogeneity", style = "cont", palette = "magma")

tm_iso = tm_shape(eco_grid_sfq2) +
  tm_polygons("distinction", style = "cont", palette = "-inferno")

tm_qua = tm_shape(eco_grid_sfq2) +
  tm_polygons("quality", style = "cont", palette = "Greens")

tm_cluster3 = tmap_arrange(tm_clu, tm_qua, tm_inh, tm_iso, ncol = 2)
tm_cluster3


