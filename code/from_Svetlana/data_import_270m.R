# 5/11/2023
# Plan: 
# Import future period annual data
# Stack annual raster data for each variable type
# Import reference period normals
# Make a coordinate list from the rasters
# Extract future data matrices from the raster stacks


##### Libraries #############
# library(rgdal)
library(terra)
library(tidyverse)
library(data.table)
##### Data ##################

# create a new sqlite database: 
#b_db<-dbConnect(RSQLite::SQLite(), "./4_the_big_analysis/data/big_db.sqlite")


# read in future climate data with terra and list.files() function: 

# make a list of files to read in
# only reading tmmin files
tmin.list<-list.files(path="./data/270m/masked/",
                      recursive = TRUE, 
                      pattern="tmin12_msk")


tmax.list<-list.files(path = "./data/270m/masked/", 
                      recursive = TRUE, 
                      pattern="tmax7_msk")

aet.list<-list.files(path = "./data/270m/masked/", 
                     recursive = TRUE, 
                     pattern="aet_msk")

def.list<-list.files(path = "./data/270m/masked/", 
                     recursive = TRUE, 
                     pattern="def_msk")



# tmin12.list<-list.files(path="./data/270m/",
#                       recursive = TRUE, 
#                       pattern="tmin12_msk")


# tmax7.list<-list.files(path = "./data/270m/", 
#                       recursive = TRUE, 
#                       pattern="tmax7_")

# read in the reference period rasters (this info will be used for validation),
# name them: 
for(f in 2:31){
  year<-1979+f
  
  tmin12<-rast(paste0("./data/270m/masked/", tmin.list[f]))
  names(tmin12)<- tmin.list[f]
  assign(paste0("tmin12.270.", year), tmin12)

  tmax7<-rast(paste0("./data/270m/masked/",tmax.list[f]))
  names(tmax7)<-tmax.list[f]
  assign(paste0("tmax7.270.", year), tmax7)

  aet<-rast(paste0("./data/270m/masked/", aet.list[f]))
  names(aet)<- aet.list[f]
  assign(paste0("aet.270.", year), aet)

  def<-rast(paste0("./data/270m/masked/",def.list[f]))
  names(def)<- tmin.list[f]
  assign(paste0("def.270.", year), def)
}

##########read in future climate rasters: 
for(f in 32:61){
  
  year<-2009+f
  
  tmin12<-rast(paste0("./data/270m/masked/", tmin.list[f]))
  names(tmin12)<-tmin.list[f]
  assign(paste0("tmin12.270.", year), tmin12)
  
  tmax7<-rast(paste0("./data/270m/masked/",tmax.list[f]))
  names(tmax7)<-tmax.list[f]
  assign(paste0("tmax7.270.", year), tmax7)
  
  aet<-rast(paste0("./data/270m/masked/",aet.list[f]))
  names(aet)<-aet.list[f]
  assign(paste0("aet.270.", year), aet)

  def<-rast(paste0("./data/270m/masked/",def.list[f]))
  names(def)<-def.list[f]
  assign(paste0("def.270.", year), def)

}


gc()


# stack future period rasters:
tmin12.list <- list()
tmax7.list<-list()
aet.list<-list()
def.list<-list()

# tmin12.list <- list()
# tmax7.list<-list()

years<-seq(2041,2070)

for(year in years){
  tmin12.list<-c(tmin12.list, get(paste0("tmin12.270.", year)))
  tmax7.list<-c(tmax7.list, get(paste0("tmax7.270.", year)))
  aet.list<-c(aet.list, get(paste0("aet.270.", year)))
  def.list<-c(def.list, get(paste0("def.270.", year)))
}

# stack: 
# tmin.stack<-rast(tmin.list)
# tmax.stack<-rast(tmax.list)
aet.stack<-rast(aet.list)
def.stack<-rast(def.list)
tmin12.stack<-rast(tmin12.list)
tmax7.stack<-rast(tmax7.list)

# extract data for future annual data: 
aet.future.df<-as.data.table(aet.stack, xy=TRUE, cells=TRUE)

# put into long format and save:

# AET
colnames(aet.future.df)[4:33]<-c(2041:2070)
aet.long<-pivot_longer(aet.future.df, c(4:33),names_to=c("year"), values_to = c("aet"))
# fwrite(aet.long, file="./4_the_big_analysis/data/aet_annual_future.csv")
# dbWriteTable(b_db, "aet_future_annual", aet.long)
saveRDS(aet.long, "./4_the_big_analysis/data/remasked/aet_annual_future.RDS")
rm(aet.future.df, aet.long)
rm(aet.stack)
gc()

# deficit
def.future.df<-as.data.table(def.stack, xy=TRUE, cells=TRUE)
colnames(def.future.df)[4:33]<-c(2041:2070)
def.long<-pivot_longer(def.future.df, c(4:33),names_to=c("year"), values_to = c("def"))
# fwrite(def.long, file="./4_the_big_analysis/data/def_annual_future.csv")
saveRDS(def.long, "./4_the_big_analysis/data/remasked/def_annual_future.RDS")

rm(def.future.df)
# dbWriteTable(b_db, "def_future_annual", def.long, overwrite=TRUE)
# dbListTables(b_db)
rm(def.long)
rm(def.stack)
gc()

# tmin12
tmin12.future.df<-as.data.table(tmin12.stack, xy=TRUE, cells=TRUE)
colnames(tmin12.future.df)[4:33]<-c(2041:2070)

tmin12.long<-pivot_longer(tmin12.future.df, c(4:33),names_to=c("year"), values_to = c("tmin12"))
rm(tmin12.future.df)
# fwrite(tmin12.long, file="./4_the_big_analysis/data/tmin12_annual_future.csv")
saveRDS(tmin12.long,  "./4_the_big_analysis/data/remasked/tmin12_annual_future.RDS")

# dbWriteTable(b_db, "tmin12_future_annual", tmin12.long, overwrite=TRUE)
# dbListTables(b_db)
rm(tmin12.long)
gc()

# tmax7
tmax7.future.df<-as.data.table(tmax7.stack, xy=TRUE, cells=TRUE)
colnames(tmax7.future.df)[4:33]<-c(2041:2070)
tmax7.long<-pivot_longer(tmax7.future.df, c(4:33),names_to=c("year"), values_to = c("tmax7"))
rm(tmax7.future.df)
# fwrite(tmax7.long, file="./4_the_big_analysis/data/tmax7_annual_future.csv")
saveRDS(tmax7.long, "./4_the_big_analysis/data/remasked/tmax7_annual_future.RDS")
# 
# dbWriteTable(b_db, "tmax7_future_annual", tmax7.long)
# dbListTables(b_db)
rm(tmax7.long)
gc()

#################################################################
############### read in reference period rasters ################


aet.270.ref<-rast("./data/270m/masked/aet_msk_ref_270.tiff")
def.270.ref<-rast("./data/270m/masked/def_msk_ref_270.tiff")
tmin.270.ref<-rast("./data/270m/masked/tmin12_msk_ref_270.tiff")
tmax.270.ref<-rast("./data/270m/masked/tmax7_msk_ref_270.tiff")


# pull out the data: 
aet_rnorms<-as.data.table(aet.270.ref, xy=TRUE, cells=TRUE)
def_rnorms<-as.data.table(def.270.ref, xy=TRUE, cells=TRUE)
tmin12_rnorms<-as.data.table(tmin.270.ref, xy=TRUE, cells=TRUE)
tmax7_rnorms<-as.data.table(tmax.270.ref, xy=TRUE, cells=TRUE)

# combine: 
ref_climate_norms<-cbind(aet_rnorms, 
                         def_rnorms[, 4], 
                         tmin12_rnorms[, 4], 
                         tmax7_rnorms[, 4])
# save
saveRDS(ref_climate_norms, "./4_the_big_analysis/data/remasked/ref_climate_normals_remasted.RDS")


############ Extract and write reference period annual data ##################

# # stack reference period rasters:
# tmin.list <- list()
# tmax.list<-list()

aet.list<-list()
def.list<-list()
tmin12.list <- list()
tmax7.list<-list()

years<-seq(1981,2010)

for(year in years){
  # tmin.list<-c(tmin.list, get(paste0("tmin.270.", year)))
  # tmax.list<-c(tmax.list, get(paste0("tmax.270.", year)))
  
  aet.list<-c(aet.list, get(paste0("aet.270.", year)))
  def.list<-c(def.list, get(paste0("def.270.", year)))
  tmin12.list<-c(tmin12.list, get(paste0("tmin12.270.", year)))
  tmax7.list<-c(tmax7.list, get(paste0("tmax7.270.", year)))
}

# stack: 
tmin12.stack<-rast(tmin12.list)
tmax7.stack<-rast(tmax7.list)
aet.stack<-rast(aet.list)
def.stack<-rast(def.list)


# Extract reference period annual data for ALL pixels: 
# put into long format and save:

# AET
aet.ref.df<-as.data.table(aet.stack, xy=TRUE, cells=TRUE)
colnames(aet.ref.df)[4:33]<-c(1981:2010)
aet.long<-pivot_longer(aet.ref.df, c(4:33),names_to=c("year"), values_to = c("aet"))
# fwrite(aet.long, file="./4_the_big_analysis/data/aet_annual_ref.csv")
saveRDS(aet.long, "./4_the_big_analysis/data/aet_annual_ref.RDS")

rm(aet.ref.df, aet.long)
rm(aet.stack)
gc()

# DEF
def.ref.df<-as.data.table(def.stack, xy=TRUE, cells=TRUE)
colnames(def.ref.df)[4:33]<-c(1981:2010)
def.long<-pivot_longer(def.ref.df, c(4:33),names_to=c("year"), values_to = c("def"))
# fwrite(def.long, file="./4_the_big_analysis/data/def_annual_ref.csv")
saveRDS(def.long, "./4_the_big_analysis/data/def_annual_ref.RDS")

rm(def.ref.df, def.long)
rm(def.stack)
gc()

# Tmin12
tmin12.ref.df<-as.data.table(tmin12.stack, xy=TRUE, cells=TRUE)
colnames(tmin12.ref.df)[4:33]<-c(1981:2010)
tmin12.long<-pivot_longer(tmin12.ref.df, c(4:33),names_to=c("year"), values_to = c("tmin12"))
# fwrite(tmin12.long, file="./4_the_big_analysis/data/tmin12_annual_ref.csv")
saveRDS(tmin12.long, "./4_the_big_analysis/data/tmin12_annual_ref.RDS")
rm(tmin12.ref.df, tmin12.long)
rm(tmin12.stack)
gc()

# Tmax7
tmax7.ref.df<-as.data.table(tmax7.stack, xy=TRUE, cells=TRUE)
colnames(tmax7.ref.df)[4:33]<-c(1981:2010)
tmax7.long<-pivot_longer(tmax7.ref.df, c(4:33),names_to=c("year"), values_to = c("tmax7"))
# fwrite(tmax7.long, file="./4_the_big_analysis/data/tmax7_annual_ref.csv")
saveRDS(tmax7.long, "./4_the_big_analysis/data/tmax7_annual_ref.RDS")
rm(tmax7.ref.df, tmax7.long)
rm(tmax7.stack)
gc()



#### combine reference annual data and add to the sqlite database: 
# 
# aet_r<-fread("./4_the_big_analysis/data/aet_annual_ref.csv", select=c("x", "y", "aet"))
# def_r<-fread("./4_the_big_analysis/data/def_annual_ref.csv", select="def")
# tmin_r<-fread("./4_the_big_analysis/data/tmin12_annual_ref.csv", select="tmin12")
# tmax_r<-fread("./4_the_big_analysis/data/tmax7_annual_ref.csv", select="tmax7")
# 
# 
# ref_annual_climate<-cbind(aet_r, def_r, tmin_r, tmax_r)
# rm(aet_r, def_r, tmin_r, tmax_r)
# gc()
# 
# b_db<-dbConnect(RSQLite::SQLite(), "./4_the_big_analysis/data/big_db.sqlite")
# dbWriteTable(b_db, "ref_annual_climate", ref_annual_climate, overwrite=TRUE)
# dbListTables(b_db)
# dbDisconnect(b_db)
# rm(ref_annual_climate)
# gc()

#################################### 
# Reference Period Normals: ####

# aet_ref_norm<-rast(paste0("./data/270m/", aet.list[1]))
# def_ref_norm<-rast(paste0("./data/270m/", def.list[1]))
# tmin12_ref_norm<-rast(paste0("./data/270m/", tmin.list[1]))
# tmax7_ref_norm<-rast(paste0("./data/270m/", tmax.list[1]))
# 
# 
# b_db<-dbConnect(RSQLite::SQLite(), "./4_the_big_analysis/data/big_db.sqlite")
# dbWriteTable(b_db, "ref_clim_270", ref_climate_norms)
# dbDisconnect(b_db)
# fwrite(ref_climate_norms, file="./4_the_big_analysis/data/normals_ref_climate.csv")

                            
######### previous code for 270 m data, unmasked #############3

# # put it into long form and save: 
# colnames(tmin.matrix.r)[2:31]<-c(1981:2010)
# tmin.long<-pivot_longer(tmin.matrix.r, c(2:31),names_to=c("year"), values_to = c("tmin"))
# 
# colnames(tmax.matrix.r)[2:31]<-c(1981:2010)
# tmax.long<-pivot_longer(tmax.matrix.r, c(2:31),names_to=c("year"), values_to = c("tmax"))
# 
# colnames(aet.matrix.r)[2:31]<-c(1981:2010)
# aet.long<-pivot_longer(aet.matrix.r, c(2:31),names_to=c("year"), values_to = c("aet"))
# 
# colnames(def.matrix.r)[2:31]<-c(1981:2010)
# def.long<-pivot_longer(def.matrix.r, c(2:31),names_to=c("year"), values_to = c("def"))
# 
# colnames(tmin12.matrix.r)[2:31]<-c(1981:2010)
# tmin12.long<-pivot_longer(tmin12.matrix.r, c(2:31),names_to=c("year"), values_to = c("tmin12"))
# 
# colnames(tmax7.matrix.r)[2:31]<-c(1981:2010)
# tmax7.long<-pivot_longer(tmax7.matrix.r, c(2:31),names_to=c("year"), values_to = c("tmax7"))
# 
# 
# ### combine and save: 
# # ref_annual_270<-cbind(aet.long, def.long[, "def"], tmin.long[, "tmin"], tmax.long[, "tmax"], tmin12.long[, "tmin12"], tmax7.long[, "tmax7"])
# # write.csv(ref_annual_270, "./outputs/ref_annual_data_270m.csv")
# 
# ref_annual_270<-read.csv("./outputs/270m_climate/ref_annual_data_270m.csv")
# ref_annual_270<-cbind(ref_annual_270, tmin12.long[############ Reference period normals data ####################################
# ## extract reference period data at 270 m resolution: 
# # extract climate data from AET raster: 
# 
# aet.r<-as.data.frame(aet.270.ref, xy=TRUE, cells=TRUE)
# def.r<-as.data.frame(def.270.ref, xy=TRUE, cells=TRUE)
# tmin.r<-as.data.frame(tmin.270.ref, xy=TRUE, cells=TRUE)
# tmax.r<-as.data.frame(tmax.270.ref, xy=TRUE, cells=TRUE)
# 
# tmin12.r<-as.data.frame(tmin12.270.ref, xy=TRUE, cells=TRUE)
# tmax7.r<-as.data.frame(tmax7.270.ref, xy=TRUE, cells=TRUE)
# 
# # combine it: 
# ref_clim270<-cbind(aet.r, def.r[, 4], tmin.r[, 4], tmax.r[, 4], tmin12.r[, 4], tmax7.r[, 4])
# colnames(ref_clim270)<-c("cell", "x", "y", "aet_270", "def_270", "tmin_270", "tmax_270", "tmin12_270", "tmax7_270")
# 
# write.csv(ref_clim270, "./outputs/reference_climate_270m.csv")

# # ref_270<-read.csv("./outputs/270m_climate/reference_climate_270m.csv")
# # ref_clim270<-cbind(ref_270, tmin12.r[, 4], tmax7.r[, 4])
# colnames(ref_clim270)[9:10]<-c("tmin12_270", "tmax7_270")
# ref_clim270<-ref_clim270[1:10]
# # put ref_clim270 into the sqlite database: 
# mydb<-dbConnect(RSQLite::SQLite(), "./outputs/new_db.sqlite")
# dbListTables(mydb)
# 
# dbWriteTable(mydb, "ref_clim270", ref_clim270[, c(2:10)], overwrite=TRUE)
# dbDisconnect(mydb)
# 
# # save tmin12 and tmax7 reference normals for transfer to scapegoat. 
# tminmax127<-cbind(tmin12.r, tmax7.r[, 4])
# colnames(tminmax127)[5]<-"tmin7_270"
# write.csv(tminmax127, "./outputs/tmin12_tmax7_ref_normals_270.csv")
