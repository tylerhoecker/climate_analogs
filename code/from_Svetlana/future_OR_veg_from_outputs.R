# get five most supported future veg projection for each pixel 

# libraries: 
library(dplyr)
library(tidyverse)

# data
path<-"/home/svetlanayegorova/Documents/OR_Future_Forests_Proj/OR_vegetation_shifts"
zonetable<-read.csv(paste0(path, "/data/R6_vegzone_code_table.txt"))



get_future_vegzone<-function(i){
  # regular data:
  raw_data2<-readRDS(paste0(path, "/scinet_output/compiled_outputs_with_veg/gl_ss_with_veg_",i,".RDS"))
  
  # validation data: 
  # raw_data2<-readRDS(paste0(path, "/scinet_output/compiled_outputs_with_veg/validation_with_sz",i,".RDS"))
  
  test<-raw_data2%>%
    filter(!is.na(focal_sz)) # remove data that falls outside of the OR border
  
  # add zone info to the table, summarize subzone and zone support: 
  dat_vegz_focal<-left_join(test, zonetable[, c("Subzone", "Vegzone")], join_by("focal_sz"=="Subzone"))
  colnames(dat_vegz_focal)[9]<-"focal_vegzone"
  dat_vegz<-left_join(dat_vegz_focal, zonetable[, c("Subzone", "Vegzone")], join_by("analog_sz"=="Subzone"))
  colnames(dat_vegz)[10]<-"analog_vegzone"
  
  # summarize future vegzone, pull out 5 top-most projected vegzones: 
  fut_veg<-dat_vegz%>%
  group_by(focal_x, focal_y, focal_vegzone, analog_vegzone)%>%
  summarize(analog_vote=n())%>%
  slice_max(order_by = analog_vote, n=5, with_ties = FALSE)
  
  saveRDS(fut_veg, paste0(path, "/scinet_output/vegzone_chunks/gl_ss_chunk_", i, ".RDS"))
  # saveRDS(fut_veg, paste0(path, "/scinet_output/vegzone_chunks/validation_chunk_", i, ".RDS"))
  print(i)
}


lapply(X=1:2, get_future_vegzone)

# are there instances where I have ties, i.e., more than 5 zones projected? 
# tie_test<-fut_veg%>%
#   group_by(focal_x, focal_y)%>%
#   summarize(n=n())%>%
#   filter(n>5)

# # pull out tied projections, keep the projection that matches the focal vegzone, otherwise, 
# # decide randomly:
# dupes<-fut_veg %>%
#   group_by(focal_x, focal_y) %>%
#   summarize(count=n())%>%
#   filter(count>1) %>%
#   ungroup()
# 
# # get focal plots with tied projections:
# fut_veg_dupes<-semi_join(fut_veg, dupes, join_by(focal_x==focal_x, focal_y==focal_y))
# 
# # if projection matches existing veg, keep that. 
# stable_veg<-fut_veg_dupes%>%
#   filter(focal_vegzone==analog_vegzone)
# 
# # for the rest of them, pick a random veg: 
# random_veg<-anti_join(fut_veg_dupes, stable_veg, join_by(focal_x==focal_x, focal_y==focal_y))
# nrow(random_veg)
# 
# random_veg<-random_veg%>%
#   slice_max(order_by=analog_vote, with_ties = FALSE)
# 
# # combine unduplicated data: 
# non_dupe_veg<-anti_join(fut_veg, dupes, join_by(focal_x==focal_x, focal_y==focal_y))
# 
 