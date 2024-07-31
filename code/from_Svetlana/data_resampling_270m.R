# this file resamples climate and vegetation data to 270 m resolution. 

# libraries: 
library(terra)



# start with vegetation data, use msk layer to mask climate data as well:
tmin_270<-rast("./data/270m/tmin_ref_270.tiff")
veg<-rast("./data/R_06_PNV/R06_PNV_20211.tif")

# complete mask raster, which excludes NAs, rock, water, developed areas: 
complete_mask<-rast("./outputs/270m_climate/complete_mask4climatedata.tiff")

#  state boundaries: 
state_bndry<-vect("./data/state_boundaries/cb_2018_us_state_5m.shp")

################################################################################

# select Oregon's boudary data only:
OR_bndry<-subset(state_bndry, state_bndry$NAME=="Oregon")

# reproject OR_bndry shapefile to match csr of the climate data: 
OR_new<-terra::project(OR_bndry, "EPSG:9674")
rm(OR_bndry)

veg_or<-crop(veg, tmin_270)

veg270<-aggregate(veg_or, fact=9, fun="modal" ) # assign the most common vegetation class among the 27 tiles considered

veg270_resamp<-resample(veg270, tmin_270, method="mode")
veg270<-veg270_resamp
plot(veg270)

# mask out pixels that are developed, rock or water
veg270_update<-mask(veg270, complete_mask, maskvalues=0)

writeRaster(veg270_update, "./outputs/270m_veg_masked.tiff", overwrite=TRUE)
# writeRaster(msk1, "./outputs/270m_developed_rock_water_maks.tiff", overwrite=TRUE)
################################################################################


# annual future climate data: 

tmin.list<-list.files(path=".",
                      recursive = TRUE, 
                      pattern="tmmin_")


tmax.list<-list.files(path = ".", 
                      recursive = TRUE, 
                      pattern="tmmax_")

aet.list<-list.files(path = ".", 
                     recursive = TRUE, 
                     pattern="aet_")

def.list<-list.files(path = ".", 
                     recursive = TRUE, 
                     pattern="deficit_")

tmin12.list<-list.files(path=".", 
                        recursive = TRUE, 
                        pattern="tmmin-12")
tmax7.list<-list.files(path=".", 
                       recursive = TRUE,
                       pattern = "tmmax-07")

# load future rasters (2041-2070) aggregate and save them

for(f in 31:60){
  
  year<-2010+f
  
  tmin<-rast(tmin.list[f])
  aggregate(tmin, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmin_", year, "_270.tiff"), overwrite=TRUE)

  tmax<-rast(tmax.list[f])
  aggregate(tmax, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmax_", year, "_270.tiff"), overwrite=TRUE)


  aet<-rast(aet.list[f])
  aggregate(aet, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "aet_", year, "_270.tiff"), overwrite=TRUE)


  def<-rast(def.list[f])
  aggregate(def, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "def_", year, "_270.tiff"), overwrite=TRUE)


  tmin12<-rast(tmin12.list[f])
  aggregate(tmin12, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmin12_", year, "_270.tiff"), overwrite=TRUE)

  tmax7<-rast(tmax7.list[f])
  aggregate(tmax7, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmax7_", year, "_270.tiff"))

  }
 


# load REFERENCE rasters (1981-2010) aggregate and save them

for(f in 1:30){
  
  year<-1980+f
  
  tmin<-rast(tmin.list[f])
  aggregate(tmin, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmin_", year, "_270.tiff"), overwrite=TRUE)

  tmax<-rast(tmax.list[f])
  aggregate(tmax, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmax_", year, "_270.tiff"))


  aet<-rast(aet.list[f])
  aggregate(aet, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "aet_", year, "_270.tiff"))


  def<-rast(def.list[f])
  aggregate(def, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "def_", year, "_270.tiff"))

  tmin12<-rast(tmin12.list[f])
  aggregate(tmin12, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmin12_", year, "_270.tiff"), overwrite=TRUE)  
  
  tmax7<-rast(tmax7.list[f])
  aggregate(tmax7, fact=3, FUN = "mean", cores=6, filename=paste0("./data/270m/", "tmax7_", year, "_270.tiff"))  
  
}


# now reference period normals: 
aet.r<-rast("./data/normals/1981-2010/aet_1981-2010.tif")
aggregate(aet.r, fact=3, FUN = "mean", cores=6, filename="./data/270m/aet_ref_270.tiff")  

def.r<-rast("./data/normals/1981-2010/deficit_1981-2010.tif")
aggregate(def.r, fact=3, FUN = "mean", cores=6, filename="./data/270m/def_ref_270.tiff")  

tmax.r<-rast("./data/normals/1981-2010/tmmax_1981-2010.tif")
aggregate(tmax.r, fact=3, FUN = "mean", cores=6, filename="./data/270m/tmax_ref_270.tiff")  

tmin.r<-rast("./data/normals/1981-2010/tmmin_1981-2010.tif")
aggregate(tmin.r, fact=3, FUN = "mean", cores=6, filename="./data/270m/tmin_ref_270.tiff")  

tmin12.r<-rast("./data/normals/1981-2010/tmmin-12.tif")
aggregate(tmin12.r, fact=3, FUN = "mean", cores=6, filename="./data/270m/tmin12_ref_270.tiff")  


tmax7.r<-rast("./data/normals/1981-2010/tmmax-07.tif")
aggregate(tmax7.r, fact=3, FUN = "mean", cores=6, filename="./data/270m/tmax7_ref_270.tiff")  



# grab the climate files that have been aggregated to 270, mask out rock, water, developed
# reference and future annual and reference means, 

# Then clip focal data to Oregon boundaries to save time. 


# aet_270.list<-list.files(path="./data/270m/",
#                          recursive = FALSE, 
#                          pattern=glob2rx("aet_*270*"))

# libraries, again: 
library(terra)

# data: 
msk<-rast("./outputs/270m_climate/complete_mask4climatedata.tiff")

# load state boundaries: 
state_bndry<-vect("./data/state_boundaries/cb_2018_us_state_5m.shp")

################################################################################

# select Oregon's boudary data only:
OR_bndry<-subset(state_bndry, state_bndry$NAME=="Oregon")
# reproject OR_bndry shapefile to match csr of the climate data: 
OR_new<-terra::project(OR_bndry, "EPSG:9674")
rm(OR_bndry)
################################################################################

# AET
aet.list1<-list.files(path="./data/270m/",
                  recursive=TRUE,
                  pattern="aet_")

list2<-list.files(path="./data/270m/",
                  recursive=TRUE,
                  pattern="msk")

aet_270.list<-setdiff(aet.list1, list2)

# Def: 
def.list1<-aet.list1<-list.files(path="./data/270m/",
                                 recursive=TRUE,
                                 pattern="def_")



def_270.list<-setdiff(def.list1, list2)


#tmin12
tmin.list1<-list.files(path="./data/270m/",
                                  recursive=TRUE,
                                  pattern="tmin12_")

tmin12_270.list<-setdiff(tmin.list1, list2)

# tmax7: 
tmax.list1<-list.files(path="./data/270m/",
                           recursive = TRUE, 
                           pattern="tmax7")

tmax7_270.list<-setdiff(tmax.list1, list2)


for(a in 1:30){
  
  year<-1980+a
  
  aet<-rast(paste0("./data/270m/", aet_270.list[a]))
  aet1<-mask(aet, msk,  maskvalues=0)
  mask(aet1, OR_new, filename=paste0("./data/270m/masked/", "aet_msk", year, "_270.tiff" ), overwrite=TRUE)

  
  def<-rast(paste0("./data/270m/",def_270.list[a]))
  def1<-mask(def, msk, maskvalues=0)
  mask(def1, OR_new, filename=paste0("./data/270m/masked/", "def_msk", year, "_270.tiff"), overwrite=TRUE)
  
  tmin<-rast(paste0("./data/270m/",tmin12_270.list[a]))
  tmin1<-mask(tmin, msk, maskvalues=0)
  mask(tmin1, OR_new, filename=paste0("./data/270m/masked/", "tmin12_msk", year ,"_270.tiff"), overwrite=TRUE)
  
  tmax<-rast(paste0("./data/270m/",tmax7_270.list[a]))
  tmax1<-mask(tmax, msk, maskvalues=0)
  mask(tmax1, OR_new, filename=paste0("./data/270m/masked/", "tmax7_msk", year,"_270.tiff" ), overwrite=TRUE)
  
}


# future annual data: 



for(a in 31:60){
  
  year<-2010+a
  
  aet<-rast(paste0("./data/270m/",aet_270.list[a]))
  aet1<-mask(aet, msk, maskvalues=0)
  mask(aet1, OR_new, filename=paste0("./data/270m/masked/", "aet_msk", year, "_270.tiff" ), overwrite=TRUE)
  
  
  def<-rast(paste0("./data/270m/",def_270.list[a]))
  def1<-mask(def, msk,  maskvalues=0)
  mask(def1, OR_new, filename=paste0("./data/270m/masked/", "def_msk", year, "_270.tiff"), overwrite=TRUE)
  
  tmin<-rast(paste0("./data/270m/",tmin12_270.list[a]))
  tmin1<-mask(tmin, msk,  maskvalues=0)
  mask(tmin1, OR_new, filename=paste0("./data/270m/masked/", "tmin12_msk", year ,"_270.tiff"), overwrite=TRUE)
  
  tmax<-rast(paste0("./data/270m/", tmax7_270.list[a]))
  tmax1<-mask(tmax, msk,  maskvalues=0)
  mask(tmax1, OR_new, filename=paste0("./data/270m/masked/", "tmax7_msk", year,"_270.tiff" ), overwrite=TRUE)
  
}

# reference normals: 


aet_ref<-rast(paste0("./data/270m/", aet_270.list[61]))
mask(aet_ref, msk,  maskvalues=0, filename=paste0("./data/270m/masked/", "aet_msk_ref_270.tiff" ), overwrite=TRUE)

def<-rast(paste0("./data/270m/", def_270.list[61]))
mask(def, msk,  maskvalues=0, filename=paste0("./data/270m/masked/", "def_msk_ref_270.tiff"), overwrite=TRUE)

tmin<-rast(paste0("./data/270m/", tmin12_270.list[61]))
mask(tmin, msk,  maskvalues=0, filename=paste0("./data/270m/masked/", "tmin12_msk_ref_270.tiff"), overwrite=TRUE)

tmax<-rast(paste0("./data/270m/",tmax7_270.list[61]))
mask(tmax, msk,  maskvalues=0, filename=paste0("./data/270m/masked/", "tmax7_msk_ref_270.tiff" ), overwrite=TRUE)

