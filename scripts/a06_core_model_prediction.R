writeLines("Loading required packages...")
library(data.table)
library(yaml) 
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(glmnet)
library(caret)

args<-commandArgs(T)
scripts_dir<-args[1]
source(paste0(scripts_dir,"pipeline_functions.R"))

config<-get_paramaters_from_config()
heading("Showing parameters from the config file...")
print(config)

heading("Reading in Predictor Inclusion Frequency Data...")
freq_data<-read.table("../p05_predictor_inclusion/st05_frequency_table.tsv",header=T)
writeLines("Showing preview of predictor inclusion frequency data...")
head(freq_data)
dim(freq_data)

#take only those in over 99% of iterations
freq_data<-freq_data[which(freq_data[,"count"]>=495),]
freq_data<-freq_data[which(freq_data[,"predictor"]!="(Intercept)"),]
head(freq_data)
dim(freq_data)

#need to make a list
core_predictors<-freq_data[,"predictor"]

#need to read in omics data
heading("Reading in Omics Data...")
omics_data<-read.table("../p02_correct_for_covariates/st02_corrected_omics_data.tsv",header=T,stringsAsFactors=F)
writeLines("Showing preview of omics data...")
head(omics_data)
dim(omics_data)

heading("Getting Core Predictor Omics Data...")
core_omics_data<-omics_data[,c("iid",as.character(core_predictors))]
writeLines("Showing preview of core omics data...")
head(core_omics_data)
dim(core_omics_data)

omics_data<-core_omics_data

#effectively a03_make_clock script

if(!is.na(config$validation_sample)){
	heading("Reading in Validation Omics Data...")
	val_omics_data<-read.table(paste0("../p02_correct_for_covariates/st02_corrected_omics_data_",config$validation_sample,".tsv"),header=T,stringsAsFactors=F)
	val_omics_data<-val_omics_data[,c("iid",as.character(core_predictors))]
	writeLines("Showing preview of Validation omics data...")
	head(val_omics_data)
	dim(val_omics_data)
}

heading("Reading in Covariate Data...")
covariate_data<-read.table("../p01_basic_qc/st01_qc_pheno_data.tsv",header=T,stringsAsFactors=F)
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


#this has all the covariates - now need just the outcome + standard + additional is not the pre correction 
heading("Getting Outcome Phenotype & Standard Covariates...")
covariate_data<-covariate_data[,c("iid",config$outcome,config$standard_covariates,config$additional_covariates)]
writeLines("Showing preview of required covariates...")
head(covariate_data)
dim(covariate_data)


if(!is.na(config$validation_sample)){
	writeLines("Showing preview of validation required covariates...")
	val_covariate_data<-val_covariate_data[,c("iid",config$outcome,config$standard_covariates,config$additional_covariates)]
	head(val_covariate_data)
	dim(val_covariate_data)
}

heading("Merging Omics & Outcome Data...")
rownames(omics_data)<-omics_data[,"iid"]
rownames(covariate_data)<-covariate_data[,"iid"]
data<-merge(omics_data,covariate_data,by="iid")
data<-data[complete.cases(data[,-1]),]
omics_data<-omics_data[,-1]
covariate_data<-covariate_data[,-1]
writeLines("Showing preview of merged data..")
head(data)
dim(data)

if(!is.na(config$validation_sample)){
	rownames(val_omics_data)<-val_omics_data[,"iid"]
	rownames(val_covariate_data)<-val_covariate_data[,"iid"]
	val_data<-merge(val_omics_data,val_covariate_data,by="iid")
	val_data<-val_data[complete.cases(val_data[,-1]),]
	val_omics_data<-val_omics_data[,-1]
	val_covariate_data<-val_covariate_data[,-1]
	writeLines("Showing preview of merged validation data..")
	head(val_data)
	dim(val_data)
}


heading("Split Into Training & Testing")
split_data<-split_training_testing(data,"core")
data_training<-split_data$training_data
data_testing<-split_data$testing_data
writeLines("Showing preview of training data...")
head(data_training)
dim(data_training)
writeLines("Showing preview of testing data...")
head(data_testing)
dim(data_testing)


heading("Split Predictors & Outcome Matrices")
training_matrix<-split_predictor_outcome_matrices(data_training,config$outcome,c(colnames(omics_data[,colnames(omics_data)!="iid"]),config$standard_covariates,config$additional_covariates))
testing_matrix<-split_predictor_outcome_matrices(data_testing,config$outcome,c(colnames(omics_data[,colnames(omics_data)!="iid"]),config$standard_covariates,config$additional_covariates))
writeLines("Showing preview of training predictor matrix...")
head(training_matrix$predictor_matrix)
dim(training_matrix$predictor_matrix)
writeLines("Showing preview of training outcome matrix...")
head(training_matrix$outcome_matrix)
dim(training_matrix$outcome_matrix)
writeLines("Showing preview of testing predictor matrix...")
head(testing_matrix$predictor_matrix)
dim(testing_matrix$predictor_matrix)
writeLines("Showing preview of testing outcome matrix...")
head(testing_matrix$outcome_matrix)
dim(testing_matrix$outcome_matrix)

if(!is.na(config$validation_sample)){
	validation_matrix<-split_predictor_outcome_matrices(val_data,config$outcome,c(colnames(val_omics_data[,colnames(val_omics_data)!="iid"]),config$standard_covariates,config$additional_covariates))
	writeLines("Showing preview of validation predictor matrix...")
	head(validation_matrix$predictor_matrix)
	dim(validation_matrix$predictor_matrix)
	writeLines("Showing preview of validation outcome matrix...")
	head(validation_matrix$outcome_matrix)
	dim(validation_matrix$outcome_matrix)
}


heading("Penalised Regression...")

penalised_regression(training_matrix$predictor_matrix,training_matrix$outcome_matrix,testing_matrix$predictor_matrix,testing_matrix$outcome_matrix,config$outcome,config$method,config$if_elastic_net)

heading("Number of Variables Selected...")
x<-n_variables_selected(1)
writeLines(x)


#want to write output stats for kk plots
for_kk_plot<-get_stats("st03_pred_obs_resid_testing_1.tsv",config$cohort,config$omic_assay)
head(for_kk_plot)
write.table(for_kk_plot,paste0("st03_obs_pred_outcome_correlation_stats_",config$cohort,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
#writeLines()
#

slope_stats<-slope_extract("st03_pred_obs_resid_testing_1.tsv",config$cohort,config$omic_assay)
head(slope_stats)
write.table(slope_stats,paste0("st03_obs_pred_outcome_fit_stats_",config$cohort,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")



cat("\n\n\nDone!\n\n\n")
