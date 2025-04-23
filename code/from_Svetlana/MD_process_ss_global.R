# 7/1/23 



# libraries
library(tidyverse)
library(MASS)
library(stats)
library(future.apply)
library(chi)
library(data.table)


# define data paths: 
#############################################################################

# define the md function: takes reference climate of all potential analogs(ref.matix) 
# and mean (mu) and inverted covariance matrix (sx) of focal future cliamte as arguments. 

the.md.function4matrix <- function(ref.matrix, mu, sx) {
  
  D <- sqrt(mahalanobis(ref.matrix, mu, sx, inverted=T))
  
  P <- round(pchi(D, df=4), 3)  #distribution function, df is #of climate variables
  
  sigma<-ifelse(P==1, 25, round(qchi(P, df=1), 2)) 
  return(list(Sigma=sigma, D=D))
  
}



# set up the "metaMDfunction", which takes row number from future_means as an argument:  
metaMDfunction1<-function(f) {
  # select a focal plot: 
  focal_f<-future_means[f, c("x", "y")]

  
  # select future annual climate data: 
  focal_future_data<-future_climr%>%
    filter(x==focal_f$x, y==focal_f$y)
  
  
# take a 10 percent random sample of the reference climate data: 
 ref_circle<-ref_clim%>%
	 slice_sample(prop=0.1)

  # organize reference climate matrix: 
  ref.matrix<-as.matrix(ref_circle[, c(4:7)])
  

  # organize focal variable means and covariance: 
  mu<-as.vector(colMeans(focal_future_data[, 3:6]))
  
  sx<-solve(cov(focal_future_data[, 3:6])) # invert the covariance matrix here to save time in MD calculation
  
  out<-the.md.function4matrix(ref.matrix, mu, sx)
  
  
  # transform output from a list to a table format: 
  out_t<-as_tibble((out))
  # colnames(out_t)<-c("P", "Sigma", "MD")
  
  # add focal and analog coordinates to the table: 
  out_t$focal_x<-pull(focal_f[1, 1])
  out_t$focal_y<-pull(focal_f[1, 2])
  
  out_t$analog_x<-ref_circle$x
  out_t$analog_y<-ref_circle$y
  
  # only select 100-most climate similar analogs: 
  out_t<-out_t%>%
    arrange(Sigma)%>%
    slice_head(n=100)
  
  # out_t<-out_t[order(Sigma)][1:200]
  rm(ref_circle, ref.matrix)
  # add new data to the output file
  
  # analog_output<-rbind(analog_output, out_t) 
  return(out_t)
  
} 


############# chunking the input and output to metaMDfunction###############

plan(callr, workers = 96)

options(future.globals.maxSize= 1000000000)


# reference climate: 
ref_clim<-readRDS("./remasked/ref_climate_normals_remasted.RDS")

# go through the chunks of future annual data: 
for(n in 1:18) {

future_climr<-readRDS(paste0("./remasked/input_chunks/remasked_future_chunk_", n, ".RDS"))
colnames(future_climr)[4:5]<-c("def", "tmin12")

# calculate future climate means: 
future_means<-future_climr%>%
  group_by(x, y)%>%
  summarise(aet_mean=mean(aet), 
            def_mean=mean(def), 
            tmin_mean=mean(tmin12), 
            tmax_mean=mean(tmax7))%>%
  ungroup()

## input data: 
## chunking up the input
  n_chunks<-ceiling(nrow(future_means)/2000)

  
# # evaluation and writing loop: 

 for(a in 1:n_chunks){
   start<-a*2000-1999
   finish<-ifelse(a<n_chunks, a*2000, nrow(future_means))
   meta_out<-future_lapply(X=start:finish, FUN=metaMDfunction1, future.seed=TRUE)
   meta_out1<-do.call("rbind", meta_out)
   print(a)
   saveRDS(meta_out1, paste0("./project_outputs/glbl_srch_otpt/gs_ss_output_",n, "_", start, "_", finish, ".RDS"))
 }

}



