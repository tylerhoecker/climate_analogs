# compile supercomputer results: 
# 8/23/23
#compile<-function(i) {

file_list<-list.files(
		     # path="./project_outputs/remasked_opt/", 
		     # path="./validation/", 
		      path="./project_outputs/glbl_srch_otpt", 
		     # pattern=paste0("valid_output_", i, "_"), 
		      pattern="gs_ss_output",
		      full.names=TRUE)

out_data<-do.call("rbind", lapply(file_list, readRDS))

saveRDS(out_data, "./compiled_outputs/glb_srch_ss_test.RDS")
#saveRDS(out_data, paste0("./compiled_outputs/validation/valid_compiled_", i, ".RDS"))
#print(i)

#}

#compile(1)
#lapply(X=1:31, FUN=compile)
#lapply(X=c(1:18), FUN=compile)
