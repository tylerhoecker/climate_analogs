aet<-readRDS("./remasked/aet_annual_future.RDS")
aet.short<-aet[, c("x", "y","aet")]
rm(aet)
saveRDS(aet.short, "./remasked/aet_annual_future_short.RDS")
rm(aet.short)
gc()
print("aet done")

def<-readRDS("./remasked/def_annual_future.RDS")
def.short<-def[, "def"]
rm(def)
saveRDS(def.short, "./remasked/def_annual_future_short.RDS")
rm(def.short)
gc()

print("def")

tmin<-readRDS("./remasked/tmin12_annual_future.RDS")
tmin.short<-tmin[, "tmin12"]
rm(tmin)
saveRDS(tmin.short, "./remasked/tmin_annual_future_short.RDS")
rm(tmin.short)
gc()

print("tmin")


tmax<-readRDS("./remasked/tmax7_annual_future.RDS")
tmax.short<-tmax[, "tmax7"]
rm(tmax)
saveRDS(tmax.short, "./remasked/tmax_annual_future_short.RDS")
rm(tmax.short)
gc()

print("tmax")

#fut_dat<-cbind(aet.short, def.short, tmin.short, tmax.short)
#fut_dat$n<-1:nrow(fut_dat)

#rm(aet.short, def.short, tmin.short, tmax.short)
#gc()

#saveRDS(fut_dat, "./remasked/fut_annual_data.RDS")
#print("saved")

# divide data into 20  chunks 


#fut_dat$ch_n<-ceiling(fut_dat$n/20)

#save_funct()<-function(n){
#	  chunk<-subset(fut_dat, ch_n=n)
#  saveRDS(chunk, paste0(".remasked/fut_annual_ch_", n, ".RDS"))
#    rm(chunk)
#    print(n)
#}

#lapply(X=1:20, FUN=save_func)

