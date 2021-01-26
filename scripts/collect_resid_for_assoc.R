library(data.table)
library(yaml) 
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(gdata)
library(Hmisc)

args<-commandArgs(T)
#args<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/clock_pheno_assoc/final/data/"
if(!file.exists(args[1])){dir.create(args[1])}
setwd(args[1])
source("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/pipeline_functions.R")

args<-c("dexa/new","horvath_cpgs","lipidomics","metabolon_metabolomics_new_new","pheno/fewer","hannum_cpgs","igg_glycomics","metabolon_complex_lipids_new","nmr","protein_new","combined_new")

heading("Get Clocks")
#read in iids

base<-fread("/exports/igmm/eddie/wilson-lab/data/base_data/orcades/phenotypes/orcades_base_phenotypes.tsv",data.table=F)
base<-base[,"iid",drop=F]

base<-base[!grepl("NIMS",base$iid),,drop=FALSE]

full<-base


template_path<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/XXX/YYY/QQQ/st03_pred_obs_resid_ZZZ_1.tsv"
template_path_2<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/XXX/YYY/QQQ/st03_ZZZ_data_1.tsv"
methods<-c("fixed_alpha") #,"lasso","cv_alpha"
data_set<-c("testing","training")
clock_types<-c("p03_make_clock","p06_core_model_prediction","a_minus_b/p03_make_clock","b_only/p03_make_clock","a_only/p03_make_clock","3_pcs/p03_make_clock","5_pcs/p03_make_clock","10_pcs/p03_make_clock","20_pcs/p03_make_clock") #

 
for(panel in args){
	writeLines(paste0(panel,"\n"))
	for(clock_type in clock_types){
		writeLines(paste0(clock_type,"\n"))
		for(method in methods){
			if(clock_type=="p06_core_model_prediction"){
				#testing
				if(file.exists(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","testing",template_path))))))){
					panel_data<-read.table(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","testing",template_path))))),header=T,stringsAsFactors=F)
					iids<-fread(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","testing",gsub("1","core",template_path_2)))))),header=T,select="iid",stringsAsFactors=F,data.table=F)
				  	panel_data<-cbind(iids,panel_data)
				  	for_full<-panel_data[,c("iid","resid")]
				  	#training
				  	panel_data<-read.table(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",template_path))))),header=T,stringsAsFactors=F)
				  	iids<-fread(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",gsub("1","core",template_path_2)))))),header=T,select="iid",stringsAsFactors=F,data.table=F)
				  	panel_data<-cbind(iids,panel_data)
				  	panel_data<-panel_data[,c("iid","resid")]
				  	for_full<-rbind(for_full,panel_data)
				  	names(for_full)[2]<-paste0(panel,"_",method,"_",clock_type)
				  	full<-merge(full,for_full,by="iid",all.x=T)
				}	
			}else{
				#testing
				if(file.exists(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","testing",template_path))))))){
					panel_data<-read.table(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","testing",template_path))))),header=T,stringsAsFactors=F)
					iids<-fread(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","testing",template_path_2))))),header=T,select="iid",stringsAsFactors=F,data.table=F)
				  	panel_data<-cbind(iids,panel_data)
				  	for_full<-panel_data[,c("iid","resid")]
				  	#training
				  	panel_data<-read.table(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",template_path))))),header=T,stringsAsFactors=F)
				  	iids<-fread(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",template_path_2))))),header=T,select="iid",stringsAsFactors=F,data.table=F)
				  	panel_data<-cbind(iids,panel_data)
				  	panel_data<-panel_data[,c("iid","resid")]
				  	for_full<-rbind(for_full,panel_data)
				  	names(for_full)[2]<-paste0(panel,"_",method,"_",clock_type)
				  	full<-merge(full,for_full,by="iid",all.x=T)
				}
			}
		}
	}
}



head(full)
dim(full)

write.table(full,"all_clock_types_resid_training_testing_13_07_2020.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#read in age_at_vene

data<-fread("/exports/igmm/eddie/wilson-lab/data/processing/orcades/phenotypes/p03_age_at_vene/age_month_of_vene.tsv",data.table=F)
head(data)
dim(data)

age<-data[,c("iid","age_at_vene")]
head(age)
dim(age)

final_data<-merge(full,age,by="iid",all.x=T)

#need to change names
colnames(final_data)<-gsub("dexa/new","dexa",colnames(final_data))
colnames(final_data)<-gsub("metabolon_metabolomics_new_new","metabolon_metabolomics_new",colnames(final_data))
colnames(final_data)<-gsub("pheno/fewer","pheno",colnames(final_data))

#need to remove slashes fromthe names
colnames(final_data)<-gsub("\\/","_",colnames(final_data))

##############################
# Need to read in age_at_dexa
##############################
#read in dob from base
#base<-fread("/exports/igmm/eddie/wilson-lab/data/base_data/orcades/phenotypes/orcades_base_phenotypes.tsv",select=c("iid","date_of_birth"),data.table=FALSE)
#date<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/dexa_date.txt",data.table=FALSE)
#df<-merge(base,date,by="iid",all=TRUE)
#df$age_at_dexa<-as.numeric(as.Date(df$date,format="%d/%m/%Y")-as.Date(df$date_of_birth,format="%d/%m/%Y"))/365
#df<-df[,c("iid","age_at_dexa")]
#final_data<-merge(final_data,df,by="iid",all=TRUE)

write.table(final_data,"all_clocks_resid_training_testing_age_13_07_2020.tsv",col.names=T,row.names=F,quote=F,sep="\t")


###########
#file for peter

#df<-fread("all_clocks_resid_training_testing_age_12_06_2020.tsv",data.table=FALSE)
#head(df)
#dim(df)

#df<-df[,grepl("fixed_alpha_p03_make_clock|age_at_vene",colnames(df))]
#head(df)
#dim(df)

#colnames(df)<-gsub("_fixed_alpha_p03_make_clock","",colnames(df))

#write.table(df,"all_clocks_resid_training_testing_age_standard_12_06_2020.tsv",col.names=T,row.names=F,quote=F,sep="\t")

