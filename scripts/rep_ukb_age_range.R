#! /usr/bin/env Rscript
writeLines("Loading required packages...")
library(data.table)
library(yaml) 

args<-commandArgs(T)
scripts_dir<-args[1]
source(paste0(scripts_dir,"/","pipeline_functions.R"))


config<-get_paramaters_from_config()
heading("Showing parameters from the config file...")
print(config)

heading("Read in Predicted Outcomes for Testing Data")
data<-fread("../p03_make_clock/st03_pred_obs_resid_testing_1.tsv",data.table=FALSE)
head(data)
dim(data)

#need to make new subset file with just those over 40
heading("Extracting Restricted Age Range from Testing Sample")
new_data<-data[data$observed_outcome>=40 & data$observed_outcome<=75,]
head(new_data)
dim(new_data)

heading("Writing Restricted Set to New File")
write.table(new_data,"st03_pred_obs_resid_testing_1_res.tsv")

#want to write output stats for kk plots
for_kk_plot<-get_stats("st03_pred_obs_resid_testing_1_res.tsv",config$cohort,config$omic_assay)
head(for_kk_plot)
write.table(for_kk_plot,paste0("st03_obs_pred_outcome_correlation_stats_",config$cohort,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
#writeLines()
#

slope_stats<-slope_extract("st03_pred_obs_resid_testing_1_res.tsv",config$cohort,config$omic_assay)
head(slope_stats)
write.table(slope_stats,paste0("st03_obs_pred_outcome_fit_stats_",config$cohort,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")


cat("\n\n\nDone!\n\n\n")

