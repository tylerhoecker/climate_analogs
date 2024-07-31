# make a mean sigma raster

# libraries
library(terra)
library(tidyverse)
library(viridis)

# data (using "incorrect" outputs while I wait for the real deal to come in)
path<-"/home/svetlanayegorova/Documents/OR_Future_Forests_Proj/OR_vegetation_shifts/"
state_bndry<-vect(paste0(path, "data/state_boundaries/cb_2018_us_state_5m.shp"))


# select Oregon's boudary data only:
OR_bndry<-subset(state_bndry, state_bndry$NAME=="Oregon")
# reproject OR_bndry shapefile to match csr of the climate data: 
OR_new<-terra::project(OR_bndry, "EPSG:9674")
rm(OR_bndry)
# raw_data<-readRDS("./scinet_output/compiled_outputs_with_veg/analog_output_with_sz3.RDS")
# raw_data2<-readRDS("./scinet_output/compiled_outputs_with_veg/analog_output_with_sz2.RDS")
# 
# # do a few checks: 
# check<-raw_data%>%
#   group_by(focal_x, focal_y)%>%
#   summarize(count=n())
# 
# range(check$count)
# # great, 100 per focal plot. 
# 
# check<-raw_data2%>%
#   filter(Sigma<=2)%>%
#   group_by(focal_x, focal_y)%>%
#   summarize(count=n())
# 
# range(check$count)
# 
# # calculate mean sigma and make a raster: 
# mean_sig1<-raw_data%>%
#   filter(Sigma<=2)%>%
#   group_by(focal_x, focal_y)%>%
#   summarize(mean_sigma=mean(Sigma), count=n())
# 
# str(mean_sig)
# mean_sig_tib<-tibble(mean_sig)
# 
# str(mean_sig_tib)
# # make a raster: 
# sig_raster<-rast(mean_sig, type="xyz", crs="EPSG:9674")
# plot(sig_raster)


# function to combine output for raster creation: 

get_sigma<-function(i){
  # regular data: 
  raw_dat<-readRDS(paste0(path, "scinet_output/compiled_otpt_glss/glb_srch_ss_", i, ".RDS"))
  
  # validation: 
  # raw_dat<-readRDS(paste0("./scinet_output/compiled_outputs_with_veg/validation_with_sz", i, ".RDS"))
  mean_sig<-raw_dat%>%
    filter(Sigma<=2)%>%
    mutate(dist=sqrt((focal_x-analog_x)^2+(focal_y-analog_y)^2)/1000)%>%
    group_by(focal_x, focal_y)%>%
    summarize(mean_sigma=mean(Sigma), 
              mean_km=mean(dist), 
              count=n())
  # saveRDS(mean_sig, paste0("./scinet_output/sigma_chunks/sigma_ch_", i, ".RDS"))
  saveRDS(mean_sig, paste0(path, "/scinet_output/sigma_chunks/gl_ss_sigma_ch_", i, ".RDS"))
  
  rm(raw_dat)
  gc()
  print(i)
}

lapply(1:18, get_sigma)

sig_out<-data.frame()

for (i in 1:18){
  sig<-readRDS(paste0(path, "scinet_output/sigma_chunks/gl_ss_sigma_ch_", i, ".RDS"))
  sig_out<-rbind(sig_out, sig)
}
head(sig_out)

# how many/where are the pixels without analogs? 
nasigma<-subset(sig_out, count<100)
nrow(nasigma)
nasig_v<-vect(nasigma, geom=c("focal_x", "focal_y"), crs="EPSG:9674")

sig_raster<-rast(sig_out, type="xyz", crs="EPSG:9674")
plot(sig_raster, "mean_sigma")
# sig_rast<-mask(sig_raster, OR_new)

par(mfrow=c(1, 1))
plot(sig_raster, "mean_sigma", col=viridis(10), main="mean sigma")
points(nasig_v, cex=0.5, col="yellow")
# writeVector(nasig_v, "./scinet_output/non_analog_pixels.shp")

plot(sig_raster, "mean_km", col=viridis(20), main="mean distance to analog")

# historgram of mean dist to analog
hist(sig_out$mean_km, breaks=50, main="mean distance to analog (glbl srch)")
mean(sig_out$mean_km) # 198 km
median(sig_out$mean_km) #182 km

hist(nasigma$count)
# add OR boundary for reference: 
plot(OR_new, add=T)

# writeRaster(sig_raster, paste0(path, "/scinet_output/sigma_raster_global_ss.tiff"), overwrite=TRUE)
# sigma<-rast("./scinet_output/sigma_raster.tiff")
# 
# sigma_df<-as.data.frame(sigma)
# saveRDS(sigma_df, "./scinet_output/sigma_df.RDS")
