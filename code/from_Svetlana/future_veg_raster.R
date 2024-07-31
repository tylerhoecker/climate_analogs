#libraries
library(terra)
library(tidyverse)

# 
path<-"/home/svetlanayegorova/Documents/OR_Future_Forests_Proj/OR_vegetation_shifts"
plurality_veg<-function(i){
  
  # validation path:
  # fut_veg<-readRDS(paste0(path, "/scinet_output/vegzone_chunks/validation_chunk_", i, ".RDS"))
  
  # regular data path: 
  fut_veg<-readRDS(paste0(path, "/scinet_output/vegzone_chunks/glss_sv_chunk_", i, ".RDS"))
  # are there instances where I have ties, i.e., more than 5 zones projected?
  tie_test<-fut_veg%>%
    group_by(focal_x, focal_y)%>%
    slice_max(analog_vote, n=1)

  # pull out tied projections, keep the projection that matches the focal vegzone, otherwise,
  # decide randomly:
  dupes<-tie_test %>%
    group_by(focal_x, focal_y) %>%
    summarize(count=n())%>%
    filter(count>1) %>%
    ungroup()
  
  # get focal plots with tied projections:
  fut_veg_dupes<-semi_join(tie_test, dupes, join_by(focal_x==focal_x, focal_y==focal_y))

  # if projection matches existing veg, keep that.
  stable_veg<-fut_veg_dupes%>%
    #filter(focal_vegzone==analog_vegzone)
  # %>%
   filter(focal_sv==analog_sv)
  # 

  
  # for the rest of them pick non-na projection or, pick a random veg:
  random_veg<-anti_join(fut_veg_dupes, stable_veg, join_by(focal_x==focal_x, focal_y==focal_y))
  nrow(random_veg)

  random_veg<-random_veg%>%
    # filter(!is.na(analog_vegzone))%>% don't remove NA's 
    slice_max(order_by=analog_vote, with_ties = FALSE)
  
  # # get the future veg that does not have ties: 
  non_dupes<-tie_test %>%
    group_by(focal_x, focal_y) %>%
    summarize(count=n())%>%
    filter(count==1) %>%
    ungroup()
  
  non_dupes1<-right_join(tie_test, non_dupes, join_by(focal_x==focal_x, focal_y==focal_y))
  
  non_dupe_veg<-rbind(non_dupes1, random_veg, stable_veg)
  
  # combine unduplicated data:
  # non_dupe_veg<-anti_join(fut_veg, dupes, join_by(focal_x==focal_x, focal_y==focal_y))
 # saveRDS(non_dupe_veg, paste0(path, "/scinet_output/winning_veg_chunks/validation_w_dist_", i, ".RDS"))
  saveRDS(non_dupe_veg, paste0(path, "/scinet_output/winning_veg_chunks/glss_sv_chunk_", i, ".RDS"))
  
  }

# test<-rast(non_dupe_veg, type="xyz", crs="EPSG:9674")
# plot(test, "analog_vegzone")
# hist(non_dupe_veg$analog_vegzone)
# plurality_veg(1)
# test<-readRDS(paste0(path, "/scinet_output/winning_veg_chunks/chunk_w_dist_1.RDS"))
# test_r<-rast(test, type="xyz", crs="EPSG:9674")
# plot(test_r)

lapply(X=1:18, plurality_veg)
# test_pred<-lapply(X=5:6, plurality_veg)

# get the saved files and make a raster: 

# valid.list1<-list.files(path=paste0(path, "/scinet_output/winning_veg_chunks/"),
#                       recursive = FALSE, 
#                       pattern="validation_chunk")

file.list1<-list.files(path=paste0(path, "/scinet_output/winning_veg_chunks/"),
                        recursive = FALSE, 
                        pattern="glss_sv_chunk_")

collect<-data.frame()

for(i in 1:18){
  inrds<-readRDS(paste0(path, "/scinet_output/winning_veg_chunks/", file.list1[i]))
  collect<-rbind(collect, inrds)
  # inrds_r<-rast(inrds, type="xyz", crs="EPSG:9674" )
  # assign(paste0("chunk", i), inrds_r)
  print(i)
}

# collect$change<-ifelse(collect$focal_vegzone==collect$analog_vegzone, 0, 1)
collect$change<-ifelse(collect$focal_sv==collect$analog_sv, 0, 1)
fut_veg<-rast(collect, type="xyz", crs="EPSG:9674")

plot(fut_veg, "analog_sv")
plot(fut_veg, "focal_sv")
writeRaster(fut_veg, paste0(path, "/scinet_output/future_simple_veg_glss.tiff") )

# ch.list<-list()
# 
# # #25 does not have regular y cell sizes
# # inrds<-readRDS(paste0(path, file.list[25]))
# # inrds_r<-rast(inrds, type="xyz", crs="EPSG:9674" )
# # assign(paste0("chunk", i), inrds_r)
# # print(i)
# # 
# 
# 
# # for(i in 1:18){
# #   ch.list<-c(ch.list, paste0("chunk", i))
# # }
# # test<-mosaic(chunk1, chunk2, chunk3, chunk4, chunk5, chunk6, chunk7, 
# #              chunk8, chunk9, chunk10, chunk11, chunk12, chunk13, chunk14, 
# #              chunk15, chunk16, chunk17, chunk18)
# 
# test<-mosaic(chunk1, chunk2, chunk3, chunk4, chunk5, chunk6, chunk7, 
#              chunk8, chunk9, chunk10, chunk11, chunk12, chunk13, chunk14, 
#              chunk15, chunk16, chunk17, chunk18)
# 
# plot(test)
# 
# test$change<-ifel(test$focal_vegzone==test$analog_vegzone, 0, 1)
# # test$change<-ifel(test$focal_sv==test$analog_sv, 0, 1)
# # writeRaster(test, "./scinet_output/simple_veg_future_veg.tiff",  datatype="INT1U", 
# #             overwrite=TRUE)
# test_valid$valid<-ifel(test_valid$focal_vegzone==test_valid$analog_vegzone, 1, 0)
# # 
# # par(mfrow=c(1, 2))
# # plot(test_valid, "valid", main="validation success")
# # plot(test, "change", main="future turnover")
# # 
# # test_valid.df<-as.data.frame(test_valid)
# # 
# # val_summary<-test_valid.df%>%
# #   group_by(change)%>%
# #   summarize(count=n())
# # 
# # val_pct_right<-1886892/sum(1886892, 1003239)
# # 
# # # save the validation raster: 
# # 
# # # writeRaster(test_valid, paste0(path, "/scinet_output/validation_raster.tiff"))
# # # read the future veg raster to avoid confusion: 
# # fut_veg<-rast(paste0(path, "/scinet_output/complete_winning_veg_raster_final.tiff"))
# # plot(fut_veg, "change")
# # # save the raster: 
# # # writeRaster(test, paste0(path, "/scinet_output/complete_winning_veg_raster_final_w_dist.tiff"), overwrite=TRUE)
# # # 
# # # # write the input file to make raster: 
# # # saveRDS(pred_bind, "./scinet_output/winning_veg_table.RDS")
# # ################################################################################
# # # par(mfrow=c(2, 2))
# # # plot(test, "analog_vote", main="analog vote", col=rev(purples(10)))
# # # plot(test, "mean_km", main="distance to used analogs")
