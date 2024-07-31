# 8/25/23 the function to process final results. 

# libraries
library(terra)
library(dplyr)

## data common to all steps of the function: 

path<-'/home/svetlanayegorova/Documents/OR_Future_Forests_Proj/OR_vegetation_shifts/'

#270m vegetation: 
# veg<-rast(paste0(path, "outputs/270m_veg_masked.tiff"))

# simplified veg: 
veg<-rast(paste0(path, "/scinet_output/simple_veg_full_raster.tiff"))
veg<-subset(veg, "SV_N")



get_veg<-function(i) {
  # get data
  test_data<-readRDS(paste0(path,"scinet_output/compiled_otpt_glss/glb_srch_ss_", i, ".RDS"))
  
  # validation data: 
  # test_data<-readRDS(paste0(path,"scinet_output/validation/valid_compiled_", i, ".RDS"))
  
  # organize to intersect with veg
  out_data1<-test_data%>%
    group_by(focal_x, focal_y)
  
  n<-nrow(out_data1)
  
  # process in two halves to avoid overwhelming RAM
  # first half: 
  focal_veg<-terra::extract(veg, out_data1[1:round(n/2), c("focal_x", "focal_y")])
  analog_veg<-terra::extract(veg, out_data1[1:round(n/2), c("analog_x", "analog_y")])
  
  # second half: 
  focal_veg1<-terra::extract(veg, out_data1[(round(n/2)+1):n, c("focal_x", "focal_y")])
  analog_veg1<-terra::extract(veg, out_data1[(round(n/2)+1):n, c("analog_x", "analog_y")])
  
  # combine
  focal_veg_all<-rbind(focal_veg, focal_veg1)
  analog_veg_all<-rbind(analog_veg, analog_veg1)
  
  # regular pnv data lines: 
  # out_data_cmplt<-cbind(out_data1, focal_veg_all[, "Subzone"], analog_veg_all[, "Subzone"])
  # colnames(out_data_cmplt)[7:8]<-c("focal_sz", "analog_sz")

  # simplified veg lines: 
  out_data_cmplt<-cbind(out_data1, focal_veg_all[, "SV_N"], analog_veg_all[, "SV_N"])
  colnames(out_data_cmplt)[7:8]<-c("focal_sv", "analog_sv")

  
  # remove focal points with NA veg: 
  # out_data_cmplt<-out_data_cmplt%>%
  #   filter(!is.na(focal_sz))
  
  # save the "full" result with 100 analogs: 
  # saveRDS(out_data_cmplt, paste0(path, "scinet_output/compiled_outputs_with_veg/validation_with_sz", i, ".RDS"))
  saveRDS(out_data_cmplt, paste0(path, "scinet_output/compiled_outputs_with_veg/gl_ss_with_simple_veg_", i, ".RDS"))
  
  print(i)
  rm(out_data_cmplt, focal_veg_all, analog_veg_all, focal_veg, focal_veg1, analog_veg, analog_veg1, test_data, out_data1)
  gc()
}

lapply(X=1:18, FUN=get_veg)
