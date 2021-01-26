library(data.table)
library(psych) # for descriptive statistics
library(ppcor) # this pacakge computes partial and semipartial correlations.
library(ggplot2)
library(corrplot)

args<-commandArgs(T)

#args<-c("/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/overlap/",FALSE)

wd<-args[1]
st<-args[2]

if(!file.exists(wd)){dir.create(wd)}

setwd(wd)
source("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/pipeline_functions.R")
args<-c("dexa/new","horvath_cpgs","lipidomics","metabolon_metabolomics_new_new","pheno/fewer","hannum_cpgs","igg_glycomics","metabolon_complex_lipids_new","nmr","protein_new") #,"combined"

heading("Gather Clock Residuals")
if(!file.exists("./a01_gather_resid")){dir.create("./a01_gather_resid")}
#read in iids

base<-fread("/exports/igmm/eddie/wilson-lab/data/base_data/orcades/phenotypes/orcades_base_phenotypes.tsv",data.table=F)
base<-base[,"iid",drop=F]

base<-base[!grepl("NIMS",base$iid),,drop=FALSE]

full<-base


template_path<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/XXX/YYY/QQQ/st03_pred_obs_resid_ZZZ_1.tsv"
template_path_2<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/XXX/YYY/QQQ/st03_ZZZ_data_1.tsv"
methods<-c("fixed_alpha") #,"lasso","cv_alpha"
data_set<-c("testing","training") #
clock_types<-c("p03_make_clock") #,"p06_core_model_prediction","a_minus_b/p03_make_clock","b_only/p03_make_clock","a_only/p03_make_clock","3_pcs/p03_make_clock","5_pcs/p03_make_clock","10_pcs/p03_make_clock","20_pcs/p03_make_clock"

 
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
				  	for_full<-panel_data[,c("iid","pred_outcome")]
				  	#training
				  	panel_data<-read.table(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",template_path))))),header=T,stringsAsFactors=F)
				  	iids<-fread(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",gsub("1","core",template_path_2)))))),header=T,select="iid",stringsAsFactors=F,data.table=F)
				  	panel_data<-cbind(iids,panel_data)
				  	panel_data<-panel_data[,c("iid","pred_outcome")]
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
				  	for_full<-panel_data[,c("iid","pred_outcome")]
				  	#training
				  	panel_data<-read.table(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",template_path))))),header=T,stringsAsFactors=F)
				  	iids<-fread(paste0(gsub("XXX",method,gsub("YYY",panel,gsub("QQQ",clock_type,gsub("ZZZ","training",template_path_2))))),header=T,select="iid",stringsAsFactors=F,data.table=F)
				  	panel_data<-cbind(iids,panel_data)
				  	panel_data<-panel_data[,c("iid","pred_outcome")]
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

write.table(full,"./a01_gather_resid/all_clock_types_resid_training_testing_09_07_2020.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#read in age_at_vene

data<-fread("/exports/igmm/eddie/wilson-lab/data/processing/orcades/phenotypes/p03_age_at_vene/age_month_of_vene.tsv",data.table=F)
head(data)
dim(data)

age<-data[,c("iid","age_at_vene")]
head(age)
dim(age)

final_data<-merge(full,age,by="iid",all.x=T)

#need to change names
#colnames(final_data)<-gsub("dexa/dexa_new","dexa_new",colnames(final_data))
colnames(final_data)<-gsub("dexa/new","dexa",colnames(final_data))
colnames(final_data)<-gsub("pheno/fewer","pheno",colnames(final_data))
colnames(final_data)<-gsub("_fixed_alpha_p03_make_clock","",colnames(final_data))

#need to remove slashes fromthe names
colnames(final_data)<-gsub("\\/","_",colnames(final_data))

write.table(final_data,"./a01_gather_resid/all_clocks_resid_training_testing_age_09_07_2020.tsv",col.names=T,row.names=F,quote=F,sep="\t")

#now have prepped all of the files
#need an option to standardise if want to standardise

traits<-colnames(final_data)[!colnames(final_data)%in%c("iid")]


mapper<-function(old_name){
	df<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/clock_names.txt",data.table=FALSE)
	new_name<-df[df$short==old_name,"long"]
	return(new_name)
}


if(st==TRUE){
	for(trait in traits){
		new_trait_name<-paste0(trait,"_s")
		sd<-sd(final_data[,trait],na.rm=TRUE)
		print(sd)
		final_data[,new_trait_name]<-final_data[,trait]/sd
	}
	df<-final_data[,c("iid",paste0(traits,"_s"))]
	names(df)<-c("iid",traits)
}else{
	df<-final_data
}

#make nice names
#old_names<-colnames(df)[!colnames(df)%in%c("iid","age_at_vene")]
#new_names<-unlist(lapply(old_names,mapper))

#colnames(df)[!colnames(df)%in%c("iid","age_at_vene")]<-new_names

#first do it all unstandardised

heading("Fitting Full Model")
if(!file.exists("./a02_full_model")){dir.create("./a02_full_model")}
#first want to fit the full model + save the output
full_model<-lm(age_at_vene~.,data=df[,-1])
summary(full_model)

coeffs<-data.frame(summary(full_model)$coefficients)
#coeffs<-as.matrix(coeffs)
coeffs<-cbind(rownames(coeffs), coeffs)
names(coeffs)<-c("predictor","beta","se","t","p")
head(coeffs)
#save the coefficients
write.table(coeffs,"./a02_full_model/full_model_coefficients.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
summary(full_model)$r.squared
text<-paste0(summary(full_model)$r.squared)
write.table(text,"./a02_full_model/full_model_r2.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
#plot the betas
ggplot(coeffs[-1,],aes(x=predictor,y=beta))+
	geom_bar(stat="identity",fill="#f8766d") +
	geom_errorbar(aes(ymin=beta-se, ymax=beta+se),width=.2,position=position_dodge(.9)) +
	coord_flip() +
	theme_classic()
ggsave("./a02_full_model/full_model_coefficients_bar_chart.pdf",height=6,width=7)

#want an outcome matrix
#want a predictor matrix
outcome_mat<-as.matrix(df[,"age_at_vene"])
predictor_mat<-as.matrix(df[,!colnames(df)%in%c("iid","age_at_vene")])

heading("Full Correlations")
if(!file.exists("./a03_full_correlations")){dir.create("./a03_full_correlations")}
#correlation matrix
cor_mat<-cor(predictor_mat,use="pairwise.complete.obs")
#write out correlation matrix
write.table(cor_mat,"./a03_full_correlations/full_clock_correlation_matrix.tsv",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

pdf("./a03_full_correlations/full_correlation_pairwise_plot.pdf")
pairs(cor_mat)
dev.off()


pdf("./a03_full_correlations/full_clock_correlation_heatmap.pdf")
corrplot(cor_mat,tl.col="black",order="hclust",hclust.method="average",tl.cex=1,addrect=3) # addrect = 4,
dev.off()

df_new<-df[complete.cases(df),]
predictor_mat<-as.matrix(df_new[,!colnames(df_new)%in%c("iid","age_at_vene")])
cor_mat<-cor(predictor_mat,use="pairwise.complete.obs")
#write out correlation matrix
write.table(cor_mat,"./a03_full_correlations/full_clock_correlation_matrix_comp_cases.tsv",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

pdf("./a03_full_correlations/full_correlation_pairwise_plot_comp_cases.pdf")
pairs(cor_mat)
dev.off()


pdf("./a03_full_correlations/full_clock_correlation_heatmap_comp_cases.pdf")
corrplot(cor_mat,tl.col="black",order="hclust",hclust.method="average",tl.cex=1,addrect=3) # addrect = 4,
dev.off()


heading("Semi Partial Correaltions")
if(!file.exists("./a04_semi_partial_correlations")){dir.create("./a04_semi_partial_correlations")}
#probelm is need complete cases for this


clocks<-colnames(df_new)[!colnames(df_new)%in%c("iid","age_at_vene")]


get_part_corrs<-function(predictor){
	outcome_mat<-as.matrix(df_new[,colnames(df_new)=="age_at_vene"])
	predictor_mat<-as.matrix(df_new[,colnames(df_new)==predictor])
	other_pred_mat<-as.matrix(df_new[,!colnames(df_new)%in%c("iid","age_at_vene",predictor)])
	res<-spcor.test(outcome_mat,predictor_mat,other_pred_mat)
	res$predictor<-predictor
	return(res)
}

res<-lapply(clocks,get_part_corrs)
res<-do.call(rbind,res)

head(res)
dim(res)
res$sr2<-res$estimate^2

write.table(res,"./a04_semi_partial_correlations/sp_correlations.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


heading("Mergeing Full Model and Semi Partial Correaltions")
x<-fread("./a02_full_model/full_model_coefficients.tsv",data.table=FALSE)
x<-x[-1,]
data<-merge(res,x,by="predictor",all=TRUE)
head(data)
dim(data)

write.table(data,"./a04_semi_partial_correlations/full_model_sp_correlations_merged.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#make plot
#make nice names for plot
res$predictor<-unlist(lapply(res$predictor,mapper))

ggplot(res,aes(x=predictor,y=sr2,fill=predictor)) +
	geom_bar(stat="identity",fill="#f8766d") +
	coord_flip() +
	theme_classic()
ggsave("./a04_semi_partial_correlations/sp_correlations_bar_chart.pdf",height=6,width=6)
