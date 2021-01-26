writeLines("Loading required packages...")
library(glmnet)
library(data.table)
library(yaml) 
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)

args<-commandArgs(T)
scripts_dir<-args[1]
source(paste0(scripts_dir,"pipeline_functions.R"))

config<-get_paramaters_from_config()
heading("Showing parameters from the config file...")
print(config)
#need to correct al the omics measures 
#read in omics data 


heading("Need to search through the files in the model_coefficients directory...")

writeLines("Reading in coefficients from model 1...")
model_data<-read.table(paste0("../p04_clock_500_iterations/st03_model_coefficients_1.tsv"))
head(model_data)
dim(model_data)

model_data<-cbind(rownames(model_data), model_data)
for_plot<-data.frame(predictor=row.names(model_data),count=0)
writeLines("Showing preview of model 1...")
head(for_plot)

heading("Now reading in coefficients from models 2-500")
	
for(i in 2:500){
  	#read in the files
  	model_data<-read.table(paste0("../p04_clock_500_iterations/st03_model_coefficients_",i,".tsv"))
  	writeLines(paste0("Successfuly read in coefficients from model ",i))
  	head(model_data)
  	dim(model_data)
  	model_data <- cbind(rownames(model_data),model_data)
  	# count the non-zero coefficients
  	for_plot$count<-ifelse(model_data$X1!=0,for_plot$count+1,for_plot$count)
  	head(for_plot)
}

writeLines("Preview of inclusion data from plot...")
head(for_plot)

heading("write frequency table out to a file...")

write.table(for_plot,paste0("st05_frequency_table.tsv"),col.names = T,row.names = F,quote = F,sep="\t")
writeLines(paste0("Inclusion frequencies written out to p05_predoctor_inclusion/st05_frequency_table.tsv"))

heading("Plot bar chart..")

# remove intercept
heading("Removing intercept from data frame for plot...")
head(for_plot)
for_plot<-for_plot[2:nrow(for_plot),]

#all those not zeros
heading("Removing all omics wiith an inclusion frquency of zero from data frmae for plot...")
for_plot<-for_plot[which(for_plot$count!=0),]
writeLines("Preview data frame for plot...")
head(for_plot)
dim(for_plot)

#need to order them largest to smallest
#population[order(population$age),]
heading("Sorting omics based on inclusion frequency...")
for_plot<-for_plot[order(for_plot$count),]
exp_stages<-for_plot$predictor
writeLines("Preview data frame for plot...")
head(for_plot)

#plot
png(paste0("st05_predictor_inclusion.png"),width=400,height=700)
plot_obj<-ggplot(data=for_plot,aes(x=predictor,y=count)) +
	geom_bar(stat="identity", fill="steelblue") +
  	ggtitle("Predictor Inclusion") +
  	xlab("Predictor") +
  	ylab("Model Inclusion Frequency") +
  	scale_x_discrete (limits =exp_stages) + 
  	coord_flip() +
  	#theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) + 
  	theme_minimal()
  print(plot_obj)
dev.off()
heading(paste0("Saved bar chart to p05_predictor_inclusion/st05_predictor_inclusion.png"))


