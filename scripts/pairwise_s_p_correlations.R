library(data.table)
library(psych) # for descriptive statistics
library(ppcor) # this pacakge computes partial and semipartial correlations.
library(ggplot2)
library(corrplot)
library(dplyr)

args<-commandArgs(T)

#args<-c("/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/overlap/pairwise/",FALSE)

wd<-args[1]
st<-args[2]

if(!file.exists(wd)){dir.create(wd)}

setwd(wd)
source("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/pipeline_functions.R")
args<-c("dexa/new","horvath_cpgs","lipidomics","metabolon_metabolomics_new_new","pheno/fewer","hannum_cpgs","igg_glycomics","metabolon_complex_lipids_new","nmr","protein_new") #,"combined"

heading("Get Clocks")
if(!file.exists("./a01_gather_data")){dir.create("./a01_gather_data")}
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

write.table(full,"./a01_gather_data/all_clock_types_resid_training_testing_29_06_2020.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

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
colnames(final_data)<-gsub("protein_new","protein",colnames(final_data))
colnames(final_data)<-gsub("pheno/fewer","pheno",colnames(final_data))
colnames(final_data)<-gsub("_fixed_alpha_p03_make_clock","",colnames(final_data))

#need to remove slashes fromthe names
colnames(final_data)<-gsub("\\/","_",colnames(final_data))

write.table(final_data,"./a01_gather_data/all_clocks_resid_training_testing_age_29_06_2020.tsv",col.names=T,row.names=F,quote=F,sep="\t")

#now

#change names
args<-gsub("dexa/new","dexa",args)
args<-gsub("pheno/fewer","pheno",args)
args<-gsub("protein_new","protein",args)

comps<-t(combn(args, 2))
write.table(comps,"./a01_gather_data/pairs.tsv",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")





heading("Partition Variance")
if(!file.exists("./a02_partition_variance")){dir.create("./a02_partition_variance")}

pairs<-fread("./a01_gather_data/pairs.tsv",header=FALSE,data.table=FALSE)

#for each pair want to calculate
#unique_1 = sr^2 1
#unique_2 = sr^2 2
#overlap (area c in text book) = 1-(a+b+d)
#all the vaiance in Y that can be expained by either x1 or x2 =a+b+c
#unexplained = 1-R2

partition_variance<-function(pair){
	clock_1<-pair[1]
	clock_2<-pair[2]
	#get just those columns from final data
	df<-final_data[,c("age_at_vene",clock_1,clock_2)]
	head(df)
	dim(df)
	#fit bivariate model
	model<-lm(age_at_vene~.,data=df)
	bi<-summary(model)$r.squared
	unex<-1-bi
	df<-df[complete.cases(df),]
	sr1<-spcor.test(df[,"age_at_vene"],df[,clock_1],df[,clock_2])$estimate
	sr1_p<-spcor.test(df[,"age_at_vene"],df[,clock_1],df[,clock_2])$p
	u1<-spcor.test(df[,"age_at_vene"],df[,clock_1],df[,clock_2])$estimate^2
	sr2<-spcor.test(df[,"age_at_vene"],df[,clock_2],df[,clock_1])$estimate
	sr2_p<-spcor.test(df[,"age_at_vene"],df[,clock_2],df[,clock_1])$p
	u2<-spcor.test(df[,"age_at_vene"],df[,clock_2],df[,clock_1])$estimate^2
	overlap<-1-(u1+u2+unex)
	total<-u1+u2+unex+overlap
	return(data.frame(
		clock_1=clock_1,
		clock_2=clock_2,
		bi=bi,
		unexplained=unex,
		sr1=sr1,
		sr1_p=sr1_p,
		u1=u1,
		sr2=sr2,
		sr2_p=sr2_p,
		u2=u2,
		overlap=overlap))
}

#res<-partition_variance(pair)


traits<-colnames(final_data)[!colnames(final_data)%in%c("iid")]

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

final_data<-df

res<-apply(pairs,1,partition_variance) #
res<-do.call(rbind,res)
head(res)

writeLines("Writing out partitioned variance results...")
write.table(res,"./a02_partition_variance/partitioned_variance_results.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


heading("Creating Overlap Barchart")
if(!file.exists("./a03_overlap_barchart")){dir.create("./a03_overlap_barchart")}

for_plot<-reshape2::melt(res[,c("clock_1","clock_2","unexplained","overlap","u1","u2")])
for_plot$Comparison<-paste0(for_plot$clock_1," vs ",for_plot$clock_2)

for_plot$Comparison<-gsub("nmr","NMR Metabolomics",for_plot$Comparison)
for_plot$Comparison<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",for_plot$Comparison)
for_plot$Comparison<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",for_plot$Comparison)
for_plot$Comparison<-gsub("igg_glycomics","UPLC IgG Glycomics",for_plot$Comparison)
for_plot$Comparison<-gsub("protein","PEA Proteomics",for_plot$Comparison)
for_plot$Comparison<-gsub("dexa","DEXA",for_plot$Comparison)
for_plot$Comparison<-gsub("pheno","Clinomics",for_plot$Comparison)
for_plot$Comparison<-gsub("lipidomics","MS Fatty Acids Lipidomics",for_plot$Comparison)
for_plot$Comparison<-gsub("horvath_cpgs","DNAme Horvath CpGs",for_plot$Comparison)
for_plot$Comparison<-gsub("hannum_cpgs","DNAme Hannum CpGs",for_plot$Comparison)
for_plot$Comparison<-gsub("combined","Mega Omics",for_plot$Comparison)

for_plot$Aspect<-for_plot$variable
for_plot$Aspect<-as.character(for_plot$Aspect)
for_plot$clock_1<-as.character(for_plot$clock_1)
for_plot$clock_2<-as.character(for_plot$clock_2)

for_plot$Aspect<-ifelse(for_plot$Aspect=="u1",for_plot$clock_1,for_plot$Aspect)
for_plot$Aspect<-ifelse(for_plot$Aspect=="u2",for_plot$clock_2,for_plot$Aspect)

for_plot$Aspect<-gsub("nmr","NMR Metabolomics",for_plot$Aspect)
for_plot$Aspect<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",for_plot$Aspect)
for_plot$Aspect<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",for_plot$Aspect)
for_plot$Aspect<-gsub("igg_glycomics","UPLC IgG Glycomics",for_plot$Aspect)
for_plot$Aspect<-gsub("protein","PEA Proteomics",for_plot$Aspect)
for_plot$Aspect<-gsub("dexa","DEXA",for_plot$Aspect)
for_plot$Aspect<-gsub("pheno","Clinomics",for_plot$Aspect)
for_plot$Aspect<-gsub("lipidomics","MS Fatty Acids Lipidomics",for_plot$Aspect)
for_plot$Aspect<-gsub("horvath_cpgs","DNAme Horvath CpGs",for_plot$Aspect)
for_plot$Aspect<-gsub("hannum_cpgs","DNAme Hannum CpGs",for_plot$Aspect)
for_plot$Aspect<-gsub("combined","Mega Omics",for_plot$Aspect)
for_plot$Aspect<-gsub("overlap","Overlap",for_plot$Aspect)
for_plot$Aspect<-gsub("unexplained","Unexplained",for_plot$Aspect)

head(for_plot)
for_plot$Aspect<-factor(for_plot$Aspect,levels=c("Overlap","NMR Metabolomics","MS Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","PEA Proteomics","MS Fatty Acids Lipidomics","DNAme Hannum CpGs","DNAme Horvath CpGs","DEXA","Clinomics","Unexplained")) #,"Mega Omics"
ggplot(for_plot,aes(x=Comparison,y=value,fill=Aspect)) +
	geom_bar(stat="identity", position = position_stack(reverse = TRUE)) + #, position = position_stack(reverse = TRUE) ,alpha=0.8
	scale_fill_manual(values=c("#8c6bb1","#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4","#bdbdbd")) + #"#ffffcc",
	theme_classic(base_size = 16) +
	theme(legend.position = "right",axis.text.x=element_text(angle = -90, hjust = 0,size=8)) + #
	ylab("Variance Explained in Chronological Age")
ggsave("./a03_overlap_barchart/new_overlap_barchart.pdf",width=15,height=7)

heading("Creating Facet Barchart")
if(!file.exists("./a04_overlap_barchart_facet")){dir.create("./a04_overlap_barchart_facet")}

head(for_plot)
dim(for_plot)

get_facet<-function(omic){
	omic_df<-for_plot %>% filter(grepl(omic,Comparison)) #& Aspect==omic
    omic_df$Comparison<-gsub(" vs ","",gsub(omic,"",omic_df$Comparison))
    omic_df$fac<-omic
    return(omic_df)
}

omics<-c("NMR Metabolomics","MS Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","PEA Proteomics","MS Fatty Acids Lipidomics","DNAme Hannum CpGs","DNAme Horvath CpGs","DEXA","Clinomics")
new_df<-lapply(omics,get_facet)
new_df<-do.call(rbind,new_df)

new_df$Aspect<-factor(new_df$Aspect,levels=c("Overlap","NMR Metabolomics","MS Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","PEA Proteomics","MS Fatty Acids Lipidomics","DNAme Hannum CpGs","DNAme Horvath CpGs","DEXA","Clinomics","Unexplained")) #"Mega Omics",
ggplot(new_df,aes(x=Comparison,y=value,fill=Aspect)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) + #, position = position_stack(reverse = TRUE) ,alpha=0.8
    facet_wrap( ~  fac,ncol=5) +
    scale_fill_manual(values=c("#8c6bb1","#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4","#bdbdbd")) + #"#ffffcc",
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",axis.text.x=element_blank()) + #axis.text.x=element_text(angle = 45, hjust = 1,size=8)
    ylab("Variance Explained in Chronological Age")
ggsave("./a04_overlap_barchart_facet/facet_overlap_barchart.pdf",height=6,width=8)

#need to write out table
out<-new_df[,c("clock_1","clock_2","Aspect","value")]

out$clock_1<-gsub("dexa","DEXA",out$clock_1)
out$clock_1<-gsub("horvath_cpgs","DNAme Horvath CpGs",out$clock_1)
out$clock_1<-gsub("hannum_cpgs","DNAme Hannum CpGs",out$clock_1)
out$clock_1<-gsub("lipidomics","MS Fatty Acids Lipidomics",out$clock_1)
out$clock_1<-gsub("pheno","Clinomics",out$clock_1)
out$clock_1<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",out$clock_1)
out$clock_1<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",out$clock_1)
out$clock_1<-gsub("protein","PEA Proteomics",out$clock_1)
out$clock_1<-gsub("nmr","NMR Metabolomics",out$clock_1)


out$clock_2<-gsub("dexa","DEXA",out$clock_2)
out$clock_2<-gsub("horvath_cpgs","DNAme Horvath CpGs",out$clock_2)
out$clock_2<-gsub("hannum_cpgs","DNAme Hannum CpGs",out$clock_2)
out$clock_2<-gsub("lipidomics","MS Fatty Acids Lipidomics",out$clock_2)
out$clock_2<-gsub("pheno","Clinomics",out$clock_2)
out$clock_2<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",out$clock_2)
out$clock_2<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",out$clock_2)
out$clock_2<-gsub("protein","PEA Proteomics",out$clock_2)
out$clock_2<-gsub("nmr","NMR Metabolomics",out$clock_2)

write.table(out,"./a04_overlap_barchart_facet/facet_overlap_barchart_table.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

################
#New
################
heading("Getting Average Overlap Per Facet")
if(!file.exists("./a06_average_overlap")){dir.create("./a06_average_overlap")}

res<-data.frame(omic=omics,av_overlap=NA)

get_av<-function(omic){
	omic_df<-new_df[new_df$fac==omic & new_df$Aspect=="Overlap",]
	av<-mean(omic_df$value)
	print(av)
	return(av)
}

lapply(omics,get_av)

res$av_overlap<-unlist(lapply(omics,get_av))

write.table(res,"./a06_average_overlap/av_overlap.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#pie chart testing
#x<-reshape2::melt(res[,!grepl("_p|sr|bi",colnames(res))])
#library(RColorBrewer)
#myPalette <- brewer.pal(4, "Dark2") 
#pie(x$value,labels=x$variable,border="white",col=myPalette)

