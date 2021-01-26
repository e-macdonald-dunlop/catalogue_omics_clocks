#! /usr/bin/env Rscript
writeLines("Loading required packages...")
library(data.table)
library(yaml) 
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(gplots)
library(corrplot)
library(e1071)
library(nFactors) 


args<-commandArgs(T)
scripts_dir<-args[1]
source(paste0(scripts_dir,"/","pipeline_functions.R"))


config<-get_paramaters_from_config()
heading("Showing parameters from the config file...")
print(config)

#now need to read in the omics phenotypes

#need to read in covariatesheading("Reading in Omics Data...")
omics_data<-get_omics_data(config$omics_file,config$to_drop)
writeLines("Showing preview of omics data...")
head(omics_data)
dim(omics_data)

if(!is.na(config$validation_sample)){
	heading("Reading in Validation Omics Data...")
	val_omics_data<-get_omics_data(config$validation_omics_file,config$to_drop)
	writeLines("Showing preview of Validation omics data...")
	head(val_omics_data)
	dim(val_omics_data)

}

heading("Reading in Covariate Data...")
covariate_data<-get_covariate_data(config$covariate_file)
writeLines("Showing preview of covariate data...")
head(covariate_data)
dim(covariate_data)

if(!is.na(config$validation_sample)){
	heading("Reading in Validation Covariate Data...")
	val_covariate_data<-get_covariate_data(config$validation_covariate_file)
	writeLines("Showing preview of Validation covariate data...")
	head(val_covariate_data)
	dim(val_covariate_data)
}


#need to sort bespoke before do anything else

omics_data<-do_bespoke(omics_data,config$bespoke_omics_rscript)
if(!is.na(config$validation_sample)){
	val_omics_data<-do_bespoke(val_omics_data,config$bespoke_omics_rscript,TRUE)
}

covariate_data<-do_bespoke(covariate_data,config$bespoke_pheno_rscript)
heading("Dropping column in Covariate Data not needed for Analysis...")
covariate_data<-covariate_data[,c("iid",config$covariate_list)]
writeLines("Showing preview of trimmed covariate data...")
head(covariate_data)
dim(covariate_data)

if(!is.na(config$validation_sample)){
	val_covariate_data<-do_bespoke(val_covariate_data,config$bespoke_pheno_rscript,TRUE)
	heading("Dropping column in Validation Covariate Data not needed for Analysis...")
	val_covariate_data<-val_covariate_data[,c("iid",config$covariate_list)]
	writeLines("Showing preview of trimmed Validation covariate data...")
	head(val_covariate_data)
	dim(val_covariate_data)
}

#this is the one to deal with missingness 
#want: missingness heatmaps + raw vs outlier threshold density plots + to write out files with no missing values 

#raw vs outlier threshold density plots
heading("Plotting Distributions of Omics Measures...")
omics_matrix<-make_distribution_plots(omics_data,config$zscore_raw)

if(!is.na(config$validation_sample)){
	heading("Plotting Distributions of Validation Omics Measures...")
	val_omics_matrix<-make_distribution_plots(val_omics_data,config$zscore_raw,config$validation_sample)
}

#plot missingness
heading("Missingness...")
writeLines("Plotting Missingness...")
missingness_by_predictor(omics_data)
missingness_by_sample(omics_data)

if(!is.na(config$validation_sample)){
	writeLines("Plotting Validation Missingness...")
  	missingness_by_predictor(val_omics_data)
  	missingness_by_sample(val_omics_data)
}
# apply zscore threshold

writeLines("Applying Z-score Filter...")
omics_data[,-1]<-apply_zscore_filter(omics_matrix,config$zscore_raw)

if(!is.na(config$validation_sample)){
	writeLines("Applying Z-score Filter to validation omics data...")
	val_omics_data[,-1]<-apply_zscore_filter(val_omics_matrix,config$zscore_raw)
}

#make NA heatmaps
writeLines("Plotting missingness heatmap pre outlier removal...")
png("NA_heatmap.png") #,width=20,height=20
plot_obj<-plot_NA_heatmap(omics_data)
print(plot_obj)
dev.off()

if(!is.na(config$validation_sample)){
	writeLines("Plotting missingness heatmap pre outlier removal for validation...")
	png(paste0("NA_heatmap_",config$validation_sample,".png")) #,width=20,height=20
	plot_obj<-plot_NA_heatmap(val_omics_data)
	print(plot_obj)
	dev.off()
}

#remove missing values
writeLines("Removing outliers...")
omics_data<-remove_single_value_predictors(omics_data)
clean_data<-omics_data

if(!is.na(config$validation_sample)){
	writeLines("Removing outliers from validation...")
	
	val_clean_data<-val_omics_data
}

heading("Creating Omics Data Desciptive Plots...")
writeLines("Plotting correlation heatmap...")
make_correlation_heatmap(omics_data)
writeLines("Plotting correlation of omics and outcome phenotype...")
make_outcome_correlation_heatmap(omics_data,covariate_data,config$outcome)
make_outcome_predictor_plots(omics_data,covariate_data,config$outcome)
writeLines("Plotting local regression plots...")
make_batch_plot(omics_data)

#here for tab one
if(!is.na(config$validation_sample)){
	writeLines("Plotting validation correlation heatmap...")
	make_correlation_heatmap(val_omics_data,config$validation_sample)
	writeLines("Plotting validation Screeplot...")
	make_screeplot(val_omics_data,config$validation_sample)
	writeLines("Plotting correlation of validation omics and outcome phenotype...")
	make_outcome_correlation_heatmap(val_omics_data,val_covariate_data,config$outcome,config$validation_sample)
	make_outcome_predictor_plots(val_omics_data,val_covariate_data,config$outcome,config$validation_sample)
	make_batch_plot(val_omics_data)
}

#write out clean data
heading("Writing Output Files... ")
write.table(clean_data,"st01_qc_omics_data.tsv",col.names=T,row.names=F,quote=F,sep="\t")
write.table(covariate_data,"st01_qc_pheno_data.tsv",col.names=T,row.names=F,quote=F,sep="\t")
writeLines("Writing QC'd omics data to p01_basic_qc/st01_qc_omics_data.tsv")
writeLines("Writing QC'd covariate data to p01_basic_qc/st01_qc_pheno_data.tsv")

if(!is.na(config$validation_sample)){
	write.table(val_clean_data,paste0("st01_qc_omics_data_",config$validation_sample,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
	write.table(val_covariate_data,paste0("st01_qc_pheno_data_",config$validation_sample,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
	writeLines(paste0("Writing QC'd validation omics data to st01_qc_omics_data_",config$validation_sample,".tsv"))
	writeLines(paste0("Writing QC'd validation covariate data to st01_qc_pheno_data_",config$validation_sample,".tsv"))
}


cat("\n\n\nDone!\n\n\n")

