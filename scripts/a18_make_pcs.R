writeLines("Loading required packages...")
library(data.table)
library(yaml) 
library(ggplot2)
library(RColorBrewer)
library(Hmisc)


args<-commandArgs(T)
scripts_dir<-args[1]
source(paste0(scripts_dir,"/","pipeline_functions.R"))

config<-get_paramaters_from_config()
heading("Showing parameters from the config file...")
print(config)

#need to correct al the omics measures 
#read in omics data 
heading("Reading in Omics Data...")
omics_data<-read.table("../p02_correct_for_covariates/st02_corrected_omics_data.tsv",header=T,stringsAsFactors=F)
writeLines("Showing preview of omics data...")
head(omics_data[,1:8])
dim(omics_data)

if(!is.na(config$validation_sample)){
	heading("Reading in Validation Omics Data...")
	val_omics_data<-read.table(paste0("../p02_correct_for_covariates/st02_corrected_omics_data_",config$validation_sample,".tsv"),header=T,stringsAsFactors=F)
	writeLines("Showing preview of Validation omics data...")
	head(val_omics_data[,1:8])
	dim(val_omics_data)
}

for(n in args[-1]){
	extract_pcs(omics_data,n)
}

if(!is.na(config$validation_sample)){
	for(n in args[-1]){
		extract_pcs(val_omics_data,n,config$validation_sample)
	}	
}


cat("\n\n\nDone!\n\n\n")
