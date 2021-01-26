library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

args<-commandArgs(T)
#args<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/figures/final/"

if(!file.exists(args[1])){dir.create(args[1])}
setwd(args[1])
source("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/pipeline_functions.R")


#compare full + core models for each 
heading("Reading in Data")
omics<-c("dexa/new","hannum_cpgs","horvath_cpgs","igg_glycomics","lipidomics","metabolon_complex_lipids_new","metabolon_metabolomics_new_new","nmr","pheno/fewer","protein_new","combined_new") #
methods<-c("fixed_alpha","lasso","cv_alpha")
clock_types<-c("p03_make_clock","p06_core_model_prediction","a_minus_b/p03_make_clock","b_only/p03_make_clock","a_only/p03_make_clock","3_pcs/p03_make_clock","5_pcs/p03_make_clock","10_pcs/p03_make_clock","20_pcs/p03_make_clock")
template<-"../../XXX/YYY/ZZZ/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"

results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric(),method=character(),clock_type=character())


for(method in methods){
	print(method)
	for(clock_type in clock_types){
		print(clock_type)
		for(omic in omics){
			print(omic)
			if(file.exists(gsub("XXX",method,gsub("ZZZ",clock_type,gsub("YYY",omic,template))))){
				print(gsub("XXX",method,gsub("ZZZ",clock_type,gsub("YYY",omic,template))))
				print("Found file!")
				df<-fread(gsub("XXX",method,gsub("ZZZ",clock_type,gsub("YYY",omic,template))),data.table=FALSE)
				head(df)
				dim(df)
				#if(omic%in%c("dexa/dexa_test","pheno/fewer")){
					#df$omic<-gsub("dexa/dexa_test","DEXA",gsub("pheno/fewer","Clinomics",df$omic))
				#}
				df$method<-method
				df$clock_type<-clock_type
				results<-rbind(results,df)
			}
		}
	}	
}

head(results)
dim(results)

heading("Plotting Method KK Plot")
if(!file.exists("./a02_method_kk_plot")){dir.create("./a02_method_kk_plot")}

## just the orcades - different colurs for different methods

df<-results[results$study=="ORCADES" & results$clock_type=="p03_make_clock",]
head(df)
dim(df)
df$omic<-gsub("_Raw","",df$omic)
df$omic<-gsub("_Updated","",df$omic)
df$omic<-gsub("Pheno","Clinomics",df$omic)
df$omic<-gsub("_"," ",df$omic)
df$omic<-gsub("protein new","Proteins",df$omic)
df$omic<-gsub("Megaomics","Mega Omics",df$omic)
df$omic<-gsub("DNA Methylation","DNAme",df$omic) 

df$method<-gsub("fixed_alpha","Fixed Alpha",df$method)
df$method<-gsub("cv_alpha","CV Alpha",df$method)
df$method<-gsub("lasso","LASSO",df$method)
names(df)<-gsub("method","Method",names(df))


get_nice_names<-function(old_name){
	df<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/clock_names.txt",data.table=FALSE)
	new_name<-df[df$long_old==old_name,"long"]
	return(new_name)
}

for(i in 1:nrow(df)){
	df[i,"omic"]<-get_nice_names(df[i,"omic"])
}

heading("Writing out of method comp data to file")
out<-df[,c("omic","study","r","lower","upper","p","Method")]
write.table(out,"./a02_method_kk_plot/method_comp_table.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=Method)) +
   		geom_pointrange(position=position_dodge(.2), stat = "identity",fatten=0.8) + #
    	#scale_x_discrete (limits =exp_stages) +
    	geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 12))  # use a white background
	return(fp)
}

#want to order by r
order_want<-df %>% filter(Method=="Fixed Alpha") %>% arrange(desc(r)) %>% select(omic)

df$Method<-factor(df$Method,levels=c("Fixed Alpha","LASSO","CV Alpha"))
df$omic<-factor(df$omic,levels=c("MS Fatty Acids Lipidomics","DEXA","MS Complex Lipidomics","NMR Metabolomics","UPLC IgG Glycomics","Clinomics","MS Metabolomics","DNAme Horvath CpGs","PEA Proteomics","DNAme Hannum CpGs","Mega Omics"))
y<-kk_plot(df)
y
ggsave("./a02_method_kk_plot/method_comp_kk.pdf",width=6,height=4)
ggsave("./a02_method_kk_plot/method_comp_kk.png",width=6,height=4) #

#want a plot for standard via core
heading("Plotting Standard vs Core")
if(!file.exists("./a03_standard_core_kk_plot")){dir.create("./a03_standard_core_kk_plot")}

df<-results[results$clock_type%in%c("p03_make_clock","p06_core_model_prediction"),]
head(df)
dim(df)


kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=clock_type)) +
   		geom_pointrange(position=position_dodge(.2), stat = "identity") + 
    	#scale_x_discrete (limits =exp_stages) +
    	scale_colour_manual(values=c("black","#f8766d")) +
    	geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 12))  # use a white background
	return(fp)
}

#panel per method

df_f<-df[df$method=="fixed_alpha",]
df_f$omic<-gsub("_Raw","",df_f$omic)
df_f$omic<-gsub("_Updated","",df_f$omic)
df_f$omic<-gsub("_"," ",df_f$omic)
df_f$omic<-gsub("Pheno","Clinomics",df_f$omic)
df_f$omic<-gsub("protein new","Proteins",df_f$omic)
df_f$omic<-gsub("Megaomics","Mega Omics",df_f$omic)
df_f$omic<-gsub("DNA Methylation","DNAme",df_f$omic)

get_nice_names<-function(old_name){
	df<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/clock_names.txt",data.table=FALSE)
	new_name<-df[df$long_old==old_name,"long"]
	return(new_name)
}

for(i in 1:nrow(df_f)){
	df_f[i,"omic"]<-get_nice_names(df_f[i,"omic"])
}


df_f$clock_type<-gsub("p03_make_clock","Standard",df_f$clock_type)
df_f$clock_type<-gsub("p06_core_model_prediction","Core",df_f$clock_type)
df_f$clock_type<-factor(df_f$clock_type,levels=c("Standard","Core"))
df_f$omic<-factor(df_f$omic,levels=c("MS Fatty Acids Lipidomics","DEXA","MS Complex Lipidomics","NMR Metabolomics","UPLC IgG Glycomics","Clinomics","MS Metabolomics","DNAme Horvath CpGs","PEA Proteomics","DNAme Hannum CpGs","Mega Omics")) #"DEXA",
#df_f<-df_f[!grepl("DEXA",df_f$omic),]


heading("Writing out data to file")
out<-df_f[,c("omic","study","r","lower","upper","p","clock_type")]
write.table(out,"./a03_standard_core_kk_plot/standard_vs_core_fixed_alpha_table.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

kk_plot(df_f)
ggsave("./a03_standard_core_kk_plot/standard_vs_core_fixed_alpha.pdf",width=6,height=4)
ggsave("./a03_standard_core_kk_plot/standard_vs_core_fixed_alpha.png",width=6,height=4) #

#########################
#Standard vs core with numbers of variabels
#########################
heading("Read in N variabeles selected data")
v<-fread("../n_variables_selected_table.tsv",data.table=FALSE)
head(v)
dim(v)

#get nice names

get_names<-function(old_name){
	df<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/clock_names.txt",data.table=FALSE)
	new_name<-df[df$short==old_name,"long"]
	return(new_name)
}

for(i in 1:nrow(v)){
	v[i,"omic"]<-get_names(v[i,"clock"])
}

get_labels<-function(omic){
	row<-v[v$omic==omic,]
	n_s<-row[,"n_sel"]
	n_c<-row[,"n_sel_core"]
	label<-paste0(omic," ",n_s,"|",n_c)
	return(label)
}

for(i in 1:nrow(df_f)){
	df_f[i,"label"]<-get_labels(df_f[i,"omic"])
}

head(df_f)

kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=label, y=r, ymin=lower, ymax=upper,colour=clock_type)) +
   		geom_pointrange(position=position_dodge(.2), stat = "identity") + 
    	#scale_x_discrete (limits =exp_stages) +
    	scale_colour_manual(values=c("black","#f8766d")) +
    	geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 12))  # use a white background
	return(fp)
}

df_f$label<-factor(df_f$label,levels=c("MS Fatty Acids Lipidomics 27|22","DEXA 28|22","MS Complex Lipidomics 130|25","NMR Metabolomics 81|43","UPLC IgG Glycomics 50|25","Clinomics 12|11","MS Metabolomics 181|40","DNAme Horvath CpGs 155|67","PEA Proteomics 203|62","DNAme Hannum CpGs 50|32","Mega Omics 214|53")) #"DEXA",
kk_plot(df_f)

ggsave("./a03_standard_core_kk_plot/standard_vs_core_fixed_alpha_with_n.pdf",width=6,height=4)
ggsave("./a03_standard_core_kk_plot/standard_vs_core_fixed_alpha_with_n.png",width=6,height=4) #

#to write out table for supp want a separate column for N predictors included

get_n<-function(vec){
	row<-v[v$omic==vec$omic,]
	n<-ifelse(vec$clock_type=="Standard",row$n_sel,row$n_sel_core)
	return(n)
}

for(i in 1:nrow(df_f)){
	df_f[i,"n_sel"]<-get_n(df_f[i,])
}

heading("Writing out data to file")
out<-df_f[,c("omic","study","r","lower","upper","p","clock_type","n_sel")]
write.table(out,"a03_standard_core_kk_plot/standard_vs_core_fixed_alpha_table_n_predictors.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#df_l<-df[df$method=="lasso",]
#kk_plot(df_l)
#ggsave("standard_vs_core_lasso.png",width=12,height=8)

#df_c<-df[df$method=="cv_alpha",]
#kk_plot(df_c)
#ggsave("standard_vs_core_cv_alpha.png",width=12,height=8)

#######################################
#comp all clock types

heading("Creating Clock Type KK Plots")
df<-results
head(df)
dim(df)

kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=clock_type)) +
   		geom_pointrange(position=position_dodge(.1), stat = "identity",fatten=0.8) + 
    	#scale_x_discrete (limits =exp_stages) +
    	geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 12))  # use a white background
	return(fp)
}

df_f<-df[df$method=="fixed_alpha",]
#df_f<-df_f[!grepl("DEXA",df_f$omic),]
df_f$clock_type<-gsub("a_minus_b/p03_make_clock","A minus B",df_f$clock_type)
df_f$clock_type<-gsub("b_only/p03_make_clock","B Only",df_f$clock_type)
df_f$clock_type<-gsub("a_only/p03_make_clock","A Only",df_f$clock_type)
df_f$clock_type<-gsub("3_pcs/p03_make_clock","3 PCs",df_f$clock_type)
df_f$clock_type<-gsub("5_pcs/p03_make_clock","5 PCs",df_f$clock_type)
df_f$clock_type<-gsub("10_pcs/p03_make_clock","10 PCs",df_f$clock_type)
df_f$clock_type<-gsub("20_pcs/p03_make_clock","20 PCs",df_f$clock_type)
df_f$clock_type<-gsub("p03_make_clock","Standard",df_f$clock_type)
df_f$clock_type<-gsub("p06_core_model_prediction","Core",df_f$clock_type)
df_f$omic<-gsub("_Raw","",df_f$omic)
df_f$omic<-gsub("Megaomics","Mega Omics",df_f$omic)
df_f$omic<-gsub("protein_new","Proteins",df_f$omic)
df_f$omic<-gsub("_Updated","",df_f$omic)
df_f$omic<-gsub("_"," ",df_f$omic)
df_f$omic<-gsub("Pheno","Clinomics",df_f$omic)
df_f$omic<-gsub("DNA Methylation","DNAme",df_f$omic)


df_f$clock_type<-gsub("p03_make_clock","Standard",df_f$clock_type)
df_f$clock_type<-gsub("p06_core_model_prediction","Core",df_f$clock_type)
#df_f$clock_type<-factor(df_f$clock_type,levels=c("Standard","Core","A minus B","B Only","A Only","3 PCs","5 PCs","10 PCs","20 PCs"))
#df_f$omic<-factor(df_f$omic,levels=c("MS Lipidomics","NMR Metabolomics","IgG Glycomics","Metabolon Complex Lipids","Clinomics","Metabolon Metabolomics","DNA Methylation Horvath CpGs","Proteins","DNA Methylation Hannum CpGs","DEXA","Mega Omics"))

#kk_plot(df_f)
#ggsave("clock_type_fixed_alpha.pdf",width=6,height=4)

#df_l<-df[df$method=="lasso",]
#kk_plot(df_l)
#ggsave("clock_type_lasso.png",width=12,height=8)

#df_c<-df[df$method=="cv_alpha",]
#kk_plot(df_c)
#ggsave("clock_type_cv_alpha.png",width=12,height=8)

heading("Creating A & B KK Plot")
if(!file.exists("./a04_AB_kk_plot")){dir.create("./a04_AB_kk_plot")}

#want to just have standard vs A/B
df_a<-df_f[df_f$clock_type%in%c("Standard","A minus B","B Only","A Only"),]

df_a$clock_type<-factor(df_a$clock_type,levels=c("Standard","A minus B","B Only","A Only"))
df_a$omic<-factor(df_a$omic,levels=c("MS Lipidomics","DEXA","NMR Metabolomics","IgG Glycomics","Metabolon Complex Lipids","Clinomics","Metabolon Metabolomics","DNAme Horvath CpGs","Proteins","DNAme Hannum CpGs","Mega Omics")) #"DEXA",

kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=clock_type)) +
   		geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	#scale_x_discrete (limits =exp_stages) +
    	scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff")) + #,"#cb181d","#7a0177"
    	#geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 11)) + # use a white background
    	theme(legend.position="bottom")
	return(fp)
}

kk_plot(df_a)
ggsave("./a04_AB_kk_plot/clock_type_fixed_alpha_A_B.pdf",width=7,height=5)
ggsave("./a04_AB_kk_plot/clock_type_fixed_alpha_A_B.png",width=7,height=5) #

heading("Creating PC KK Plot")
if(!file.exists("./a05_PC_kk_plot")){dir.create("./a05_PC_kk_plot")}

#want to just have standard vs PCs
df_p<-df_f[df_f$clock_type%in%c("Standard","3 PCs","5 PCs","10 PCs","20 PCs"),]


get_nice_names<-function(old_name){
	df<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/clock_names.txt",data.table=FALSE)
	new_name<-df[df$long_old==old_name,"long"]
	return(new_name)
}

for(i in 1:nrow(df_p)){
	df_p[i,"omic"]<-get_nice_names(df_p[i,"omic"])
}


df_p$clock_type<-factor(df_p$clock_type,levels=c("Standard","3 PCs","5 PCs","10 PCs","20 PCs"))
df_p$omic<-factor(df_p$omic,levels=c("MS Fatty Acids Lipidomics","DEXA","MS Complex Lipidomics","NMR Metabolomics","UPLC IgG Glycomics","Clinomics","MS Metabolomics","DNAme Horvath CpGs","PEA Proteomics","DNAme Hannum CpGs","Mega Omics")) #"DEXA",

kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=clock_type)) +
   		geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	#scale_x_discrete (limits =exp_stages) +
    	scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff")) + #,"#cb181d","#7a0177"
    	#geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 11)) +  # use a white background
    	theme(legend.position="bottom")
	return(fp)
}

heading("Write out data to file")
out<-df_p[,c("omic","study","r","lower","upper","p","clock_type")]
write.table(out,"./a05_PC_kk_plot/clock_type_fixed_alpha_PCs_table.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

kk_plot(df_p)
ggsave("./a05_PC_kk_plot/clock_type_fixed_alpha_PCs.pdf",width=7,height=5)
ggsave("./a05_PC_kk_plot/clock_type_fixed_alpha_PCs.png",width=7,height=5) #

#################################################

heading("Creating Validation KK Plot")
if(!file.exists("./a06_validation_kk_plot")){dir.create("./a06_validation_kk_plot")}

#validation kk plot
results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric(),method=character())

#read in fixed alpha results
input_files<-c(
	"../../fixed_alpha/dexa/new/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/dexa/new/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/dexa/new/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	"../../fixed_alpha/hannum_cpgs/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/hannum_cpgs/validation/gs/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_GS.tsv",
	"../../fixed_alpha/hannum_cpgs/validation/egcut/p03_make_clock/st03_obs_pred_outcome_correlation_stats_EGC_methylation_hannum.tsv",
	"../../fixed_alpha/horvath_cpgs/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/horvath_cpgs/validation/gs/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_GS.tsv",
	"../../fixed_alpha/horvath_cpgs/validation/egcut/p03_make_clock/st03_obs_pred_outcome_correlation_stats_EGC_methylation_horvath.tsv",
	"../../fixed_alpha/igg_glycomics/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/igg_glycomics/validation/korcula/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_Korcula.tsv",
	"../../fixed_alpha/igg_glycomics/validation/vis/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_Vis.tsv",
	"../../fixed_alpha/nmr/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/nmr/validation/korcula/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_Korcula.tsv",
	"../../fixed_alpha/nmr/validation/egcut/p03_make_clock/st03_obs_pred_outcome_correlation_stats_EGC_NMR.tsv",
	"../../fixed_alpha/CVD_INF/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/CVD_INF/validation/vis/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_Vis.tsv",
	"../../fixed_alpha/EGCUT_panels/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/EGCUT_panels/validation/egcut/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_EGC_prote.tsv",
	"../../fixed_alpha/pheno/fewer/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/pheno/fewer/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	"../../fixed_alpha/pheno/fewer/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv") #,"../../fixed_alpha/dexa/new/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv","../../fixed_alpha/dexa/new/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv", "../../fixed_alpha/dexa/new/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv", "../../fixed_alpha/dexa/dexa_new/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv","../../fixed_alpha/dexa/dexa_new/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv","../../fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB.tsv","../../fixed_alpha/dexa/ukbb/validation/orcades/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"

for(file in input_files){
	df<-fread(file,data.table=F)
	df$method<-"fixed_alpha"
	results<-rbind(results,df)
}

df<-results

df$omic<-gsub("_Raw","",df$omic)
df$omic<-gsub("Pheno","Clinomics",df$omic)
df$study<-gsub("ORCADES_res","ORCADES Age Restricted",df$study)
df$omic<-gsub("CVD_INF","PEA Proteomics Subset 1",df$omic)
df$omic<-gsub("EGCUT_panels","PEA Proteomics Subset 2",df$omic)
df$omic<-gsub("OLINK_proteins","PEA Proteomics Subset 2",df$omic)
df$omic<-gsub("IgG_Glycomics","UPLC IgG Glycomics",df$omic)
df$omic<-gsub("_"," ",df$omic)
df$omic<-ifelse(df$study=="EGC_methylation_hannum","DNA Methylation Hannum CpGs",df$omic)
df$omic<-ifelse(df$study=="EGC_methylation_horvath","DNA Methylation Horvath CpGs",df$omic)

df$study<-gsub("EGC_prote","EBB",df$study)
df$study<-gsub("EGC_methylation_hannum","EBB",df$study)
df$study<-gsub("EGC_methylation_horvath","EBB",df$study)
df$study<-gsub("EGC_NMR","EBB",df$study)
#df$study<-gsub("EGC_methylation_horvath","EGCUT",df$study)
df$study<-gsub("GS","GS:SFHS",df$study)
df$study<-factor(df$study,levels=c("ORCADES","Vis","Korcula","GS:SFHS","EBB","UKBB","ORCADES Age Restricted"))
df$omic<-gsub("DNA Methylation","DNAme",df$omic)
#want to order to make ORCADES first 

kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=study)) +
   		geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	#scale_x_discrete (limits =exp_stages) +
    	scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff","#ff61c9","grey")) +
    	geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 11)) + # use a white background
    	theme(legend.position="bottom")
	return(fp)
}


heading("Write out data to file")
out<-df[,c("omic","study","r","lower","upper","p")]
write.table(out,"./a06_validation_kk_plot/validation_kk_plot_table.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


kk_plot(df)
ggsave("./a06_validation_kk_plot/validation_kk_plot.pdf",width=7,height=3)
ggsave("./a06_validation_kk_plot/validation_kk_plot.png",width=7,height=3) #

#need to split into panels of how good the replication is
#need ot make facets
df[df$omic%in%c("DEXA","NMR Metabolomics"),"fac"]<-"c"
df[df$omic%in%c("Clinomics","UPLC IgG Glycomics"),"fac"]<-"b"
df[df$omic%in%c("DNAme Hannum CpGs","DNAme Horvath CpGs","PEA Proteomics Subset 1","PEA Proteomics Subset 2"),"fac"]<-"a"

kk_plot(df) + facet_grid(fac~.,scales="free_y",space="free_y") + theme(strip.background = element_blank(),strip.text.y=element_blank())
ggsave("./a06_validation_kk_plot/validation_kk_plot_fac.pdf",width=7,height=4)
ggsave("./a06_validation_kk_plot/validation_kk_plot_fac.png",width=7,height=4) #

heading("Creating Training Testing KK Plot")
if(!file.exists("./a13_training_testing_kk_plot")){dir.create("./a13_training_testing_kk_plot")}
args<-c("dexa/new","hannum_cpgs","horvath_cpgs","igg_glycomics","lipidomics","metabolon_complex_lipids_new","metabolon_metabolomics_new_new","nmr","pheno/fewer","protein_new","combined_new") #

long_data<-data.frame(observed_outcome=numeric(),pred_outcome=numeric(),resid=numeric(),clock=character(),sample=character())

anat<-data.frame(omic=character(),sample=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric())


cor_extract<-function(data){
	return(cor.test(as.numeric(data$pred_outcome),as.numeric(data$observed_outcome)))
}


get_annotation_x<-function(data_testing,data_training,clock){
	cor_testing<-cor_extract(data_testing)
	cor_training<-cor_extract(data_training)
	anat<-data.frame(matrix(NA,nrow=2,ncol=6))
	names(anat)<-c("omic","sample","r","lower","upper","p")
	anat$sample<-c("Testing","Training")
	anat$r<-c(cor_testing$estimate,cor_training$estimate)
	anat$lower<-c(cor_testing$conf.int[1],cor_training$conf.int[1])
	anat$upper<-c(cor_testing$conf.int[2],cor_training$conf.int[2])
	anat$p<-c(cor_testing$p.value,cor_training$p.value)
	anat$omic<-clock
	
	return(anat)
}

for(panel in args){
	data_testing<-read.table(paste0("../../fixed_alpha/",panel,"/p03_make_clock/st03_pred_obs_resid_testing_1.tsv"),header=T,stringsAsFactors = F)
	data_testing[,'clock']<-panel
	data_testing[,'sample']<-"Testing"
	data_training<-read.table(paste0("../../fixed_alpha/",panel,"/p03_make_clock/st03_pred_obs_resid_training_1.tsv"),header=T,stringsAsFactors = F)
	data_training[,'clock']<-panel
	data_training[,'sample']<-"Training"
	long_data<-rbind(long_data,data_testing,data_training)
	anat_panel<-get_annotation_x(data_testing,data_training,panel)
	anat<-rbind(anat,anat_panel)
	print(head(anat))
}


get_nice_names<-function(old_name){
	df<-fread("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/clock_names.txt",data.table=FALSE)
	new_name<-df[df$short==old_name,"long"]
	return(new_name)
}

for(i in 1:nrow(anat)){
	anat[i,"omic"]<-get_nice_names(anat[i,"omic"])
}


kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=sample)) +
   		geom_pointrange(position=position_dodge(.2), stat = "identity") + #,fatten=0.8
   		scale_colour_manual(values=c("black","#bdbdbd")) + #,"#cb181d","#7a0177"
    	#scale_x_discrete (limits =exp_stages) +
    	geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	coord_flip() +  # flip coordinates (puts labels on y axis)
    	xlab(xlab) + ylab(ylab) +
    	#theme(axis.text.x = element_text(colour=a) ) + 
    	theme_bw((base_size = 12))  # use a white background
	return(fp)
}

anat$sample<-factor(anat$sample,levels=c("Testing","Training"))
anat$omic<-factor(anat$omic,levels=c("MS Fatty Acids Lipidomics","DEXA","NMR Metabolomics","UPLC IgG Glycomics","MS Complex Lipidomics","Clinomics","MS Metabolomics","DNAme Horvath CpGs","PEA Proteomics","DNAme Hannum CpGs","Mega Omics")) #"DEXA",

kk_plot(anat)
ggsave("./a13_training_testing_kk_plot/training_testing_kk.pdf",width=6,height=4)
ggsave("./a13_training_testing_kk_plot/training_testing_kk.png",width=6,height=4) #

heading("Write table out to file")
out<-anat
write.table(out,"./a13_training_testing_kk_plot/training_testing_table.tsv",col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
#need to do a separate dexa plot for all of the different variations


#heading("Creating Validation KK Plot")
#if(!file.exists("./a07_dexa_kk_plot")){dir.create("./a07_dexa_kk_plot")}

#get the testing stats

#results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric(),method=character())

#paths<-c("/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_pred_obs_resid_training.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/dexa_new/p03_make_clock/st03_pred_obs_resid_training_1.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/new/p03_make_clock/st03_pred_obs_resid_training_1.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/new/random/p03_make_clock/st03_pred_obs_resid_training.tsv")
#omics<-c("DEXA_ukbb","DEXA","DEXA_u","DEXA_ran")
#studies<-c("UKBB Training","ORCADES Training","ORCADES Training","ORCADES Training")

#get_stats<-function(resid_file,study,omic){
	#df<-fread(resid_file,data.table=F)
	#x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	#r<-x$estimate
	#lower<-x$conf.int[1]
	#upper<-x$conf.int[2]
	#p<-x$p.value
	#return(data.frame(omic=omic,study=study,r=r,lower=lower,upper=upper,p=p))
#}

#for(i in 1:length(paths)){
	#res<-get_stats(paths[i],studies[i],omics[i])
	#res$method<-"fixed_alpha"
	#results<-rbind(results,res)
#}

#read in fixed alpha results
#input_files<-c(
	#"../../fixed_alpha/dexa/new/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	#"../../fixed_alpha/dexa/new/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	#"../../fixed_alpha/dexa/new/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	#"../../fixed_alpha/dexa/dexa_new/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	#"../../fixed_alpha/dexa/dexa_new/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	#"../../fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	#"../../fixed_alpha/dexa/ukbb/validation/orcades/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	#"../../fixed_alpha/dexa/new/random/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv",
	#"../../fixed_alpha/dexa/new/random/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	#"../../fixed_alpha/dexa/new/random/p19_testing_ukb_range/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"
	#)

#for(file in input_files){
	#df<-fread(file,data.table=F)
	#df$method<-"fixed_alpha"
	#results<-rbind(results,df)
#}

#df<-results


#kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	#fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=study)) +
   		#geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	##scale_x_discrete (limits =exp_stages) +
    	##scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff","")) +
    	##geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	#coord_flip() +  # flip coordinates (puts labels on y axis)
    	#xlab(xlab) + ylab(ylab) +
    	##theme(axis.text.x = element_text(colour=a) ) + 
    	#theme_bw((base_size = 11)) + # use a white background
    	#theme(legend.position="bottom")
	#return(fp)
#}

#kk_plot(df)
#ggsave("./a07_dexa_kk_plot/dexa_kkplot.pdf",width=7,height=4)



#heading("Creating Validation KK Plot")
#if(!file.exists("./a08_dexa_kk_plot")){dir.create("./a08_dexa_kk_plot")}

#args<-c(
#"age_res_ran",
#"age_res_ran_date",
#"age_res_ran_date_occ",
#"age_res_ran_occ",
#"age_res_ran_out",
#"age_res_ran_out_date",
#"age_res_ran_out_date_occ",
#"age_res_ran_out_occ",
#"age_res_s",
#"age_res_s_date",
#"age_res_s_date_occ",
#"age_res_s_occ",
#"age_res_s_out",
#"age_res_s_out_date",
#"age_res_s_out_date_occ",
#"age_res_s_out_occ",
#"ran",
#"ran_date",
#"ran_date_occ",
#"ran_occ",
#"ran_out",
#"ran_out_date",
#"ran_out_date_occ",
#"ran_out_occ",
#"s",
#"s_date",
#"s_date_occ",
#"s_occ",
#"s_out",
#"s_out_date",
#"s_out_date_occ",
#"s_out_occ",
#"new"
#)

#template_testing<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"
#template_training<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES_training.tsv"
#template_val<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv"

#get the testing stats
#need to read in the testing correlations

#results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric())

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_testing),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"ORCADES Testing"
	#results<-rbind(results,df)
#}

#for(arg in args){
	#if(!file.exists(gsub("XXX",arg,template_training))){
		#next
	#}else{
		#df<-fread(gsub("XXX",arg,template_training),data.table=FALSE)
		#df$omic<-arg
		#df$study<-"ORCADES Training"
		#results<-rbind(results,df)
	#}
#}

#need to get the validation

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_val),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"UKBB"
	#results<-rbind(results,df)
#}

#need to get the ukbb training testing validation
#paths<-c(
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB_training.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/validation/orcades/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv")
#studies<-c("UKBB Testing","UKBB Training","ORCADES")


#for(i in 1:length(paths)){
	#df<-fread(paths[i],data.table=FALSE)
	#df$omic<-"ukbb"
	#df$study<-studies[i]
	#results<-rbind(results,df)
#}

#need to get the training for the s ones

#ss<-args[!grepl("ran",args)]

#template_x<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_pred_obs_resid_training_1.tsv"

#get_stats<-function(omic){
	#df<-fread(gsub("XXX",omic,template_x),data.table=F)
	#x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	#r<-x$estimate
	#lower<-x$conf.int[1]
	#upper<-x$conf.int[2]
	#p<-x$p.value
	#return(data.frame(omic=omic,study="ORCADES Training",r=r,lower=lower,upper=upper,p=p))
#}

#for(i in 1:length(ss)){
	#res<-get_stats(ss[i])
	#results<-rbind(results,res)
#}

#df<-results

#want to add facets
#date
#occ
#date_occ
#new
#ukbb

#df$fac<-NA
#df[grepl("date",df$omic),"fac"]<-"date"
#df[grepl("occ",df$omic),"fac"]<-"occ"
#df[grepl("date_occ",df$omic),"fac"]<-"date_occ"
#df[grepl("new",df$omic),"fac"]<-"new"
#df[grepl("ukbb",df$omic),"fac"]<-"ukbb"


#kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	#fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=study)) +
   		#geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	##scale_x_discrete (limits =exp_stages) +
    	##scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff","")) +
    	##geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	#coord_flip() +  # flip coordinates (puts labels on y axis)
    	#xlab(xlab) + ylab(ylab) +
    	##theme(axis.text.x = element_text(colour=a) ) + 
    	#theme_bw((base_size = 11)) + # use a white background
    	#theme(legend.position="bottom")
	#return(fp)
#}

#kk_plot(df) + facet_wrap(~fac,nrow=1,scales="free")
#ggsave("./a08_dexa_kk_plot/dexa_kkplot.pdf",width=15,height=5)


##########################################

#heading("Creating Validation KK Plot")
#if(!file.exists("./a09_dexa_date_kk_plot")){dir.create("./a09_dexa_date_kk_plot")}

#args<-c(
#"s",
#"s_date",
#"s_date_occ",
#"s_occ",
#"s_out",
#"s_out_date",
#"s_out_date_occ",
#"s_out_occ",
#"new",
#"new_date",
#"date_month",
#"date_season",
#"date_season_year",
#"no_tot_lean"
#)

#template_testing<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"#
#template_training<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES_training.tsv"
#template_val<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv"

#get the testing stats
#need to read in the testing correlations

#results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric())

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_testing),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"ORCADES Testing"
	#results<-rbind(results,df)
#}

#for(arg in args){
	#if(!file.exists(gsub("XXX",arg,template_training))){
		#next
	#}else{
		#df<-fread(gsub("XXX",arg,template_training),data.table=FALSE)
		#df$omic<-arg
		#df$study<-"ORCADES Training"
		#results<-rbind(results,df)
	#}
#}

#need to get the validation

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_val),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"UKBB"
	#results<-rbind(results,df)
#}

#need to get the ukbb training testing validation
#paths<-c(
#	"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",#
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB_training.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/validation/orcades/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv")
#studies<-c("UKBB Testing","UKBB Training","ORCADES")


#for(i in 1:length(paths)){
	#df<-fread(paths[i],data.table=FALSE)
	#df$omic<-"ukbb"
	#df$study<-studies[i]
	#results<-rbind(results,df)
#}

#need to get the training for the s ones

#ss<-args[!grepl("ran",args)]

#template_x<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_pred_obs_resid_training_1.tsv"

#get_stats<-function(omic){
	#df<-fread(gsub("XXX",omic,template_x),data.table=F)
	#x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	#r<-x$estimate
	#lower<-x$conf.int[1]
	#upper<-x$conf.int[2]
	#p<-x$p.value
	#return(data.frame(omic=omic,study="ORCADES Training",r=r,lower=lower,upper=upper,p=p))
#}

#for(i in 1:length(ss)){
	#res<-get_stats(ss[i])
	#results<-rbind(results,res)
#}

#df<-results

#want to add facets
#date
#occ
#date_occ
#new
#ukbb

#df$fac<-NA
#df[grepl("date",df$omic),"fac"]<-"date"
#df[grepl("occ",df$omic),"fac"]<-"occ"
#df[grepl("date_occ",df$omic),"fac"]<-"date_occ"
#df[grepl("new",df$omic),"fac"]<-"new"
#df[grepl("ukbb",df$omic),"fac"]<-"ukbb"


#kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	#fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=study)) +
   		#geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	##scale_x_discrete (limits =exp_stages) +
    	##scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff","")) +
    	##geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	#coord_flip() +  # flip coordinates (puts labels on y axis)
    	#xlab(xlab) + ylab(ylab) +
    	##theme(axis.text.x = element_text(colour=a) ) + 
    	#theme_bw((base_size = 11)) + # use a white background
    	#theme(legend.position="bottom")
	#return(fp)
#}

#kk_plot(df) + facet_wrap(~fac,nrow=1,scales="free")
#ggsave("././a09_dexa_date_kk_plot/dexa_kk_plot.pdf",width=15,height=5)



######################################

##########################################

#heading("Creating Validation KK Plot")
#if(!file.exists("./a10_dexa_date_kk_plot")){dir.create("./a10_dexa_date_kk_plot")}

#args<-c(
#"age_res_ran",
#"age_res_ran_date",
#"age_res_ran_date_occ",
#"age_res_ran_occ",
#"age_res_ran_out",
#"age_res_ran_out_date",
#"age_res_ran_out_date_occ",
#"age_res_ran_out_occ",
#"age_res_s",
#"age_res_s_date",
#"age_res_s_date_occ",
#"age_res_s_occ",
#"age_res_s_out",
#"age_res_s_out_date",
#"age_res_s_out_date_occ",
#"age_res_s_out_occ",
#"ran",
#"ran_date",
#"ran_date_occ",
#"ran_occ",
#"ran_out",
#"ran_out_date",
#"ran_out_date_occ",
#"ran_out_occ",
#"s",
#"s_date",
#"s_date_occ",
#"s_occ",
#"s_out",
#"s_out_date",
#"s_out_date_occ",
#"s_out_occ",
#"new",
#"new_date",
#"date_month",
#"date_season",
#"date_season_year",
#"no_tot_lean",
#"age_res_date_year",
#"age_res_date_year_both"
#)

#template_testing<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"
#template_training<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES_training.tsv"
#template_val<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv"

#get the testing stats
#need to read in the testing correlations

#results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric())

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_testing),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"ORCADES Testing"
	#results<-rbind(results,df)
#}

#for(arg in args){
	#if(!file.exists(gsub("XXX",arg,template_training))){
		#next
	#}else{
		#df<-fread(gsub("XXX",arg,template_training),data.table=FALSE)
		#df$omic<-arg
		#df$study<-"ORCADES Training"
		#results<-rbind(results,df)
	#}
#}

#need to get the validation

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_val),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"UKBB"
	#results<-rbind(results,df)
#}

#need to get the ukbb training testing validation
#paths<-c(
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB_training.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/validation/orcades/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv")
#studies<-c("UKBB Testing","UKBB Training","ORCADES")


#for(i in 1:length(paths)){
	#df<-fread(paths[i],data.table=FALSE)
	#df$omic<-"ukbb"
	#df$study<-studies[i]
	#results<-rbind(results,df)
#}

#need to get the training for the s ones

#ss<-args[!grepl("ran",args)]

#template_x<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_pred_obs_resid_training_1.tsv"

#get_stats<-function(omic){
	#df<-fread(gsub("XXX",omic,template_x),data.table=F)
	#x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	#r<-x$estimate
	#lower<-x$conf.int[1]
	#upper<-x$conf.int[2]
	#p<-x$p.value
	#return(data.frame(omic=omic,study="ORCADES Training",r=r,lower=lower,upper=upper,p=p))
#}

#for(i in 1:length(ss)){
	#res<-get_stats(ss[i])
	#results<-rbind(results,res)
#}

#df<-results

#want to add facets
#date
#occ
#date_occ
#new
#ukbb

#df$fac<-NA
#for_pj<-df[!grepl("ran",df$omic) & !grepl("occ",df$omic) & !grepl("both",df$omic),]
#for_pj$omic<-gsub("new_date","date_year",for_pj$omic)
#for_pj$omic<-gsub("new","zscore_3",for_pj$omic)
#for_pj[grepl("date",for_pj$omic),"fac"]<-"date"
#for_pj[grepl("zscore",for_pj$omic),"fac"]<-"zscore"
#for_pj[grepl("ukbb",for_pj$omic),"fac"]<-"ukbb"


#kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	#fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=study)) +
   		#geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	##scale_x_discrete (limits =exp_stages) +
    	##scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff","")) +
    	##geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	#coord_flip() +  # flip coordinates (puts labels on y axis)
    	#xlab(xlab) + ylab(ylab) +
    	##theme(axis.text.x = element_text(colour=a) ) + 
    	#theme_bw((base_size = 11)) + # use a white background
    	#theme(legend.position="bottom")
	#return(fp)
#}

#kk_plot(for_pj) + facet_wrap(~fac,nrow=1,scales="free")
#ggsave("./a10_dexa_date_kk_plot/dexa_kk_plot_pj.pdf",width=10,height=7)



######################################

##########################################

#heading("Creating Validation KK Plot")
#if(!file.exists("./a11_dexa_date_kk_plot")){dir.create("./a11_dexa_date_kk_plot")}

#args<-c(
#"s",
#"s_occ",
#"s_out",
#"s_out_occ",
#"new",
#"s_out",
#"new_date",
#"date_month",
#"date_season",
#"date_season_year"
#)

#template_testing<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"
#template_training<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_obs_pred_outcome_correlation_stats_ORCADES_training.tsv"
#template_val<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/validation/ukbb/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv"
#template_core_testing<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p06_core_model_prediction/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv"
#template_core_val<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/validation/ukbb_core/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_UKBB.tsv"
#get the testing stats
#need to read in the testing correlations

#results<-data.frame(omic=character(),study=character(),r=numeric(),lower=numeric(),upper=numeric(),p=numeric())

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_testing),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"ORCADES Testing"
	#results<-rbind(results,df)
#}

#for(arg in args){
	#if(!file.exists(gsub("XXX",arg,template_training))){
		#next
	#}else{
		#df<-fread(gsub("XXX",arg,template_training),data.table=FALSE)
		#df$omic<-arg
		#df$study<-"ORCADES Training"
		#results<-rbind(results,df)
	#}
#}

#need to get the validation

#for(arg in args){
	#df<-fread(gsub("XXX",arg,template_val),data.table=FALSE)
	#df$omic<-arg
	#df$study<-"UKBB"
	#results<-rbind(results,df)
#}

#need to get the ukbb training testing validation
#paths<-c(
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/p03_make_clock/st03_obs_pred_outcome_correlation_stats_UKBB_training.tsv",
	#"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/ukbb/validation/orcades/p03_predict_in_validation/st03_obs_pred_outcome_correlation_stats_ORCADES.tsv")
#studies<-c("UKBB Testing","UKBB Training","ORCADES")


#for(i in 1:length(paths)){
	#df<-fread(paths[i],data.table=FALSE)
	#df$omic<-"ukbb"
	#df$study<-studies[i]
	#results<-rbind(results,df)
#}

# need core

#for(arg in args){
	#if(!file.exists(gsub("XXX",arg,template_core_testing))){
		#next
	#}else{
		#df<-fread(gsub("XXX",arg,template_core_testing),data.table=FALSE)
		#df$omic<-paste0(arg,"_core")
		#df$study<-"ORCADES Testing"
		#results<-rbind(results,df)
	#}
#}

#need core validation
#for(arg in args){
	#if(!file.exists(gsub("XXX",arg,template_core_val))){
		#next
	#}else{
		#df<-fread(gsub("XXX",arg,template_core_val),data.table=FALSE)
		#df$omic<-paste0(arg,"_core")
		#df$study<-"UKBB"
		#results<-rbind(results,df)
	#}
#}


#need to get the training for the s ones

#ss<-args[!grepl("ran",args)]

#template_x<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p03_make_clock/st03_pred_obs_resid_training_1.tsv"

#get_stats<-function(omic){
	#df<-fread(gsub("XXX",omic,template_x),data.table=F)
	#x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	#r<-x$estimate
	#lower<-x$conf.int[1]
	#upper<-x$conf.int[2]
	#p<-x$p.value
	#return(data.frame(omic=omic,study="ORCADES Training",r=r,lower=lower,upper=upper,p=p))
#}

#for(i in 1:length(ss)){
	#res<-get_stats(ss[i])
	#results<-rbind(results,res)
#}

#get core training
#template_y<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/fixed_alpha/dexa/XXX/p06_core_model_prediction/st03_pred_obs_resid_training_1.tsv"
#s<-c("new","s_out")

#get_stats<-function(omic){
	#df<-fread(gsub("XXX",omic,template_y),data.table=F)
	#x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	#r<-x$estimate
	#lower<-x$conf.int[1]
	#upper<-x$conf.int[2]
	#p<-x$p.value
	#return(data.frame(omic=paste0(omic,"_core"),study="ORCADES Training",r=r,lower=lower,upper=upper,p=p))
#}

#for(i in 1:length(s)){
	#res<-get_stats(s[i])
	#results<-rbind(results,res)
#}

#df<-results

#want to add facets
#date
#occ
#date_occ
#new
#ukbb

#df$fac<-NA
#for_pj<-df[!grepl("ran",df$omic) & !grepl("occ",df$omic) & !grepl("both",df$omic),]
#for_pj$omic<-gsub("new_date","date_year",for_pj$omic)
#for_pj$omic<-gsub("new","zscore_3",for_pj$omic)
#for_pj[grepl("date",for_pj$omic),"fac"]<-"date"
#for_pj[grepl("zscore",for_pj$omic),"fac"]<-"zscore"
#for_pj[grepl("ukbb",for_pj$omic),"fac"]<-"ukbb"


#kk_plot<-function(df=NULL,xlab="Clock",ylab="r (95% CI)",vline=0){
	#fp <- ggplot(data=df, aes(x=omic, y=r, ymin=lower, ymax=upper,colour=study)) +
   		#geom_pointrange(position=position_dodge(.1), stat = "identity") + #,fatten=0.8
    	##scale_x_discrete (limits =exp_stages) +
    	##scale_colour_manual(values=c("black","#f8766d","#00ba38","#619cff","#c77cff","")) +
    	##geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
    	#coord_flip() +  # flip coordinates (puts labels on y axis)
    	#xlab(xlab) + ylab(ylab) +
    	##theme(axis.text.x = element_text(colour=a) ) + 
    	#theme_bw((base_size = 11)) + # use a white background
    	#theme(legend.position="bottom")
	#return(fp)
#}

#kk_plot(for_pj) #+ facet_wrap(~fac,nrow=1,scales="free")
#ggsave("./a11_dexa_date_kk_plot/dexa_kk_plot_pj.pdf",width=7,height=5)
