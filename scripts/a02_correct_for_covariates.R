#! /usr/bin/env Rscript
writeLines("Loading required packages...")
library(data.table)
library(yaml) 
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(gplots)
library(corrplot)
library(e1071)
library(nFactors) 

args<-commandArgs(T)
scripts_dir<-args[1]
source(paste0(scripts_dir,"pipeline_functions.R"))

config<-get_paramaters_from_config()
heading("Showing parameters from the config file...")
print(config)

#need to correct all the omics measures 
#read in omics data 
heading("Reading in Omics Data...")
omics_data<-read.table("../p01_basic_qc/st01_qc_omics_data.tsv",header=T,stringsAsFactors=F)
writeLines("Showing preview of omics data...")
head(omics_data)
dim(omics_data)

if(!is.na(config$validation_sample)){
	heading("Reading in Validation Omics Data...")
	val_omics_data<-fread(paste0("../p01_basic_qc/st01_qc_omics_data_",config$validation_sample,".tsv"),data.table=FALSE)
	writeLines("Showing preview of Validation omics data...")
	head(val_omics_data)
	dim(val_omics_data)
}

heading("Reading in Covariate Data...")
covariate_data<-fread("../p01_basic_qc/st01_qc_pheno_data.tsv",data.table=FALSE)
writeLines("Showing preview of covariate data...")
head(covariate_data)
dim(covariate_data)

if(!is.na(config$validation_sample)){
	heading("Reading in Validation Covariate Data...")
	val_covariate_data<-read.table(paste0("../p01_basic_qc/st01_qc_pheno_data_",config$validation_sample,".tsv"),header=T,stringsAsFactors=F)
	writeLines("Showing preview of Validation covariate data...")
	head(val_covariate_data)
	dim(val_covariate_data)
}


##
# want to test the association of covariates with omic
#want to output the summary
# + drop covariates if not signif in majority of omics


if(length(config$pre_correction_covariates)>0){
	covs<-config$pre_correction_covariates
	#assoc_covariates(omics_data,covariate_data,covs)
	#covariate_data<-do_bespoke(covariate_data,config$bespoke_pheno_rscript)
	output<-correct_data(omics_data,covariate_data,covs)
	corrected_data<-output$corr_data
	#########
	heading("Plotting Distributions of Corrected Omics Measures...")
	corrected_matrix<-make_distribution_plots(corrected_data,config$zscore_resid)

	# apply zscore threshold
	writeLines("Applying Z-score Filter...")
	if(!is.na(config$zscore_resid)){
		corrected_data[,-1]<-apply_zscore_filter(corrected_matrix,config$zscore_resid)
	}
	#make NA heatmaps
	#writeLines("Plotting missingness heatmap pre outlier removal...")
	#png("NA_heatmap.png") #,width=20,height=20
	#plot_obj<-plot_NA_heatmap(corrected_data)
	#print(plot_obj)
	#dev.off()

	write.table(corrected_data,"st02_corrected_omics_data_pre_scale.tsv",col.names=T,row.names=F,quote=F,sep="\t")

	#remove missing values
	writeLines("Removing outliers...")
	corrected_data<-do_bespoke_qc(corrected_data,config$bespoke_omics_qc_rscript)
	#writeLines("Plotting missingness heatmap post outlier removal...")
	#png("NA_heatmap_post_zscore.png") #,width=20,height=20
	#plot_obj<-plot_NA_heatmap(corrected_data)
	#print(plot_obj)
	#dev.off()

	heading("Creating Omics Data Desciptive Plots...")
	writeLines("Plotting correlation heatmap...")
	make_correlation_heatmap(corrected_data)
	writeLines("Plotting Screeplot...")
	make_screeplot(corrected_data)


	#########
	corrected_data[-1]<-scale_centre(corrected_data[-1])
	make_outcome_predictor_plots(corrected_data,covariate_data,config$outcome)
	corrected_data<-corrected_data[colSums(!is.na(corrected_data)) > 0]
	write.table(corrected_data,"st02_corrected_omics_data.tsv",col.names=T,row.names=F,quote=F,sep="\t")
	heading("Writing out files..")
	writeLines("Writing corrected omics data to st02_corrected_omics_data_.tsv")
	#do_clustering(corrected_data)
	if(!is.na(config$validation_sample)){
		#covs<-config$pre_correction_covariates
		#assoc_covariates(val_omics_data,val_covariate_data,covs)
		val_output<-correct_data(val_omics_data,val_covariate_data,output$covariates_chosen,sample="validation")
		val_corrected_data<-val_output$corr_data
		#########
		heading("Plotting Distributions of Corrected Omics Measures...")
		val_corrected_matrix<-make_distribution_plots(val_corrected_data,config$zscore_resid,config$validation_sample)

		# apply zscore threshold
		writeLines("Applying Z-score Filter...")
		if(!is.na(config$zscore_resid)){
			val_corrected_data[,-1]<-apply_zscore_filter(val_corrected_matrix,config$zscore_resid)
		}
		#make NA heatmaps
		writeLines("Plotting missingness heatmap pre outlier removal...")
		png(paste0("NA_heatmap_",config$validation_sample,".png")) #,width=20,height=20
		plot_obj<-plot_NA_heatmap(val_corrected_data)
		print(plot_obj)
		dev.off()

		#remove missing values
		writeLines("Removing outliers...")
		val_corrected_data<-remove_NAs_validation(val_corrected_data,config$discovery_predictors)
		writeLines("Plotting missingness heatmap post outlier removal...")
		png("NA_heatmap_post_zscore.png") #,width=20,height=20
		plot_obj<-plot_NA_heatmap(val_corrected_data)
		print(plot_obj)
		dev.off()

		#########
	
		val_corrected_data[-1]<-scale_centre(val_corrected_data[-1])
		val_corrected_data<-val_corrected_data[colSums(!is.na(val_corrected_data)) > 0]
		write.table(val_corrected_data,paste0("st02_corrected_omics_data_",config$validation_sample,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
		heading("Writing out files..")
		writeLines(paste0("Writing corrected omics data to st02_corrected_omics_data_",config$validation_sample,".tsv"))
		#do_clustering(val_corrected_data)
	}
}else{
	if(!file.exists("../p01_basic_qc/st01_qc_omics_data.tsv"))stop("Error: No qc omics file found in p01_basic_qc!")
	######
	heading("Plotting Distributions of Corrected Omics Measures...")
	omics_matrix<-make_distribution_plots(omics_data,config$zscore_resid)

	# apply zscore threshold
	writeLines("Applying Z-score Filter...")
	if(!is.na(config$zscore_resid)){
		omics_data[,-1]<-apply_zscore_filter(omics_matrix,config$zscore_resid)
	}
	#make NA heatmaps
	writeLines("Plotting missingness heatmap pre outlier removal...")
	png("NA_heatmap.png") #,width=20,height=20
	plot_obj<-plot_NA_heatmap(omics_data)
	print(plot_obj)
	dev.off()

	#remove missing values
	writeLines("Removing outliers...")
	omics_data<-do_bespoke_qc(omics_data,config$bespoke_omics_qc_rscript)
	writeLines("Plotting missingness heatmap post outlier removal...")
	png("NA_heatmap_post_zscore.png") #,width=20,height=20
	plot_obj<-plot_NA_heatmap(omics_data)
	print(plot_obj)
	dev.off()

	heading("Creating Omics Data Desciptive Plots...")
	writeLines("Plotting correlation heatmap...")
	make_correlation_heatmap(omics_data)
	writeLines("Plotting Screeplot...")
	make_screeplot(omics_data)

	######
	
	omics_data[-1]<-scale_centre(omics_data[-1])
	omics_data<-omics_data[colSums(!is.na(omics_data)) > 0]
	write.table(omics_data,"st02_corrected_omics_data.tsv",col.names=T,row.names=F,quote=F,sep="\t")
	heading("Writing out files..")
	writeLines("No covariates specified, copying ../p01_basic_qc/st01_qc_omics_data.tsv to p02_correct_for_covariates/st02_corrected_omics_data.tsv...")
	#do_clustering(omics_data)
	if(!is.na(config$validation_sample)){
		if(!file.exists(paste0("../p01_basic_qc/st01_qc_omics_data_",config$validation_sample,".tsv")))stop("Error: No qc omics file found in p01_basic_qc!")
		######
		heading("Plotting Distributions of Corrected Omics Measures...")
		val_omics_matrix<-make_distribution_plots(val_omics_data,config$zscore_resid,config$validation_sample)

		# apply zscore threshold
		writeLines("Applying Z-score Filter...")
		if(!is.na(config$zscore_resid)){
			val_omics_data[,-1]<-apply_zscore_filter(val_omics_matrix,config$zscore_resid)
		}
		#make NA heatmaps
		writeLines("Plotting missingness heatmap pre outlier removal...")
		png(paste0("NA_heatmap_",config$validation_sample,".png")) #,width=20,height=20
		plot_obj<-plot_NA_heatmap(val_omics_data)
		print(plot_obj)
		dev.off()

		#remove missing values
		writeLines("Removing outliers...")
		val_omics_data<-remove_NAs_validation(val_omics_data,config$discovery_predictors)
		writeLines("Plotting missingness heatmap post outlier removal...")
		png("NA_heatmap_post_zscore.png") #,width=20,height=20
		plot_obj<-plot_NA_heatmap(val_omics_data)
		print(plot_obj)
		dev.off()

		######
	
		val_omics_data[-1]<-scale_centre(val_omics_data[-1])
		val_omics_data<-val_omics_data[colSums(!is.na(val_omics_data)) > 0]
		write.table(val_omics_data,paste0("st02_corrected_omics_data_",config$validation_sample,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
		heading("Writing out files..")
		writeLines(paste0("No covariates specified, copying ../p01_basic_qc/st01_qc_omics_data_",config$validation_sample,".tsv to p02_correct_for_covariates/st02_corrected_omics_data_",config$validation_sample,".tsv..."))
		#do_clustering(val_omics_data)
	}
}


cat("\n\n\nDone!\n\n\n")

