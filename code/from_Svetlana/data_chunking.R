 







save_func<-function(a, datain){ 
  chunk2write<-subset(datain, ch_n==a)
  saveRDS(chunk2write, paste0("./remasked/ref_annual/", clim_var, "_ref_annual_ch_", a, ".RDS"))
  rm(chunk2write)
  gc()
   print(a)
}

chunk<-173395*30


clim<-c("aet", "def", "tmin12", "tmax7")
period<-"ref"

# aet: 
datain<-readRDS(paste0("./remasked/ref_annual/", clim[1], "_annual_", period, ".RDS"))
clim_var<-clim[1]
datain<-datain[, c(2:3,5)]
datain$n<-1:nrow(datain)
datain$ch_n<-ceiling(datain$n/chunk)

lapply(X=1:31, FUN=save_func, datain)
rm(datain)
gc()

for(i in 2:4){
# the rest of them
datain<-readRDS(paste0("./remasked/ref_annual/", clim[i], "_annual_", period, ".RDS"))
clim_var<-clim[i]
datain$n<-1:nrow(datain)
datain<-datain[, c(5:6)]

datain$ch_n<-ceiling(datain$n/chunk)

lapply(X=1:31, FUN=save_func, datain)
rm(datain)
gc()
}


# tmin: 
#datain<-readRDS(paste0("./remasked/", clim[3], "_annual_future.RDS"))
#clim_var<-clim[3]
#datain$n<-1:nrow(datain)
#datain<-datain[, c(5:6)]

#datain$ch_n<-ceiling(datain$n/chunk)

#lapply(X=1:31, FUN=save_func, datain)
#rm(datain)
#gc()


# tmax
#datain<-readRDS(paste0("./remasked/", clim[4], "_annual_future.RDS"))
#clim_var<-clim[4]
#datain$n<-1:nrow(datain)
#datain<-datain[, c(5:6)]

#datain$ch_n<-ceiling(datain$n/chunk)

#lapply(X=1:31, FUN=save_func, datain)
#rm(datain)
#gc()





 

