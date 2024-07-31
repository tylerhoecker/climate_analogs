period<-"ref"

for(i in 1:18){

aet<-readRDS(paste0("./remasked/ref_annual/aet_", period, "_annual_ch_", i, ".RDS"))
def<-readRDS(paste0("./remasked/ref_annual/def_", period, "_annual_ch_", i, ".RDS"))
tmin<-readRDS(paste0("./remasked/ref_annual/tmin12_", period, "_annual_ch_", i, ".RDS"))
tmax<-readRDS(paste0("./remasked/ref_annual/tmax7_", period, "_annual_ch_", i, ".RDS"))

chunk<-cbind(aet[, c("x", "y", "aet")], def[, "def"], tmin[, "tmin12"], tmax[, c("tmax7", "ch_n")])

saveRDS(chunk, paste0("./remasked/ref_annual/", period, "_chunk_", i, ".RDS"))

}
