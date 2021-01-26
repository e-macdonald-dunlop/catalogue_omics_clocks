library(data.table)
library(openxlsx)

setwd("/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/biomarker_dictionaries/")
source("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/pipeline_functions.R")

heading("Proteins")

#need to read in protein dictionary
prot<-read.xlsx("all_protein_dictionary.xlsx")

#read in main protein clock
clock<-fread("../fixed_alpha/protein_new/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/protein_new/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")

sub_1<-fread("../fixed_alpha/CVD_INF/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(sub_1)
dim(sub_1)

sub_1<-sub_1[3:nrow(sub_1),]
names(sub_1)<-c("predictor","effect_size")


sub_2<-fread("../fixed_alpha/EGCUT_panels/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(sub_2)
dim(sub_2)

sub_2<-sub_2[3:nrow(sub_2),]
names(sub_2)<-c("predictor","effect_size")

out<-prot

for(i in 1:nrow(out)){
	print(out[i,"variable_name"])
	out[i,"pass_qc"]<-ifelse(out[i,"variable_name"]%in%clock$predictor,1,0)
	if(out[i,"variable_name"]%in%clock$predictor){
		out[i,"selected"]<-ifelse(clock[clock$predictor==out[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		out[i,"selected"]<-0
	}
	out[i,"available_core"]<-ifelse(out[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(out[i,"variable_name"]%in%clock_core$predictor){
		out[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==out[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		out[i,"selected_core"]<-0
	}
	out[i,"avail_sub_1"]<-ifelse(out[i,"variable_name"]%in%sub_1$predictor,1,0)
	if(out[i,"variable_name"]%in%sub_1$predictor){
		out[i,"selected_sub_1"]<-ifelse(sub_1[sub_1$predictor==out[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		out[i,"selected_sub_1"]<-0
	}
	out[i,"avail_sub_2"]<-ifelse(out[i,"variable_name"]%in%sub_2$predictor,1,0)
	if(out[i,"variable_name"]%in%sub_2$predictor){
		out[i,"selected_sub_2"]<-ifelse(sub_2[sub_2$predictor==out[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		out[i,"selected_sub_2"]<-0
	}
}


wb <- createWorkbook()
sheet_name<-"PEA Proteomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,out)

heading("DNA Methylation")

dname<-read.xlsx("cpg_list.xlsx")

hannum<-fread("../fixed_alpha/hannum_cpgs/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(hannum)
dim(hannum)

hannum<-hannum[3:nrow(hannum),]
names(hannum)<-c("predictor","effect_size")

hannum_core<-fread("../fixed_alpha/hannum_cpgs/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(hannum_core)
dim(hannum_core)

hannum_core<-hannum_core[3:nrow(hannum_core),]
names(hannum_core)<-c("predictor","effect_size")

horvath<-fread("../fixed_alpha/horvath_cpgs/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(horvath)
dim(horvath)

horvath<-horvath[3:nrow(horvath),]
names(horvath)<-c("predictor","effect_size")

horvath_core<-fread("../fixed_alpha/horvath_cpgs/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(horvath_core)
dim(horvath_core)

horvath_core<-horvath_core[3:nrow(horvath_core),]
names(horvath_core)<-c("predictor","effect_size")

for(i in 1:nrow(dname)){
	if(dname[i,"Paper"]=="Hannum"){
		dname[i,"pass_qc"]<-ifelse(dname[i,"Marker"]%in%hannum$predictor,1,0)
		if(dname[i,"Marker"]%in%hannum$predictor){
			dname[i,"selected"]<-ifelse(hannum[hannum$predictor==dname[i,"Marker"],"effect_size"]!=0,1,0)
		}else{
			dname[i,"selected"]<-0
		}
		dname[i,"available_core"]<-ifelse(dname[i,"Marker"]%in%hannum_core$predictor,1,0)
		if(dname[i,"Marker"]%in%hannum_core$predictor){
			dname[i,"selected_core"]<-ifelse(hannum_core[hannum_core$predictor==dname[i,"Marker"],"effect_size"]!=0,1,0)
		}else{
			dname[i,"selected_core"]<-0
		}
	}
	if(dname[i,"Paper"]=="Horvath"){
		dname[i,"pass_qc"]<-ifelse(dname[i,"Marker"]%in%horvath$predictor,1,0)
		if(dname[i,"Marker"]%in%horvath$predictor){
			dname[i,"selected"]<-ifelse(horvath[horvath$predictor==dname[i,"Marker"],"effect_size"]!=0,1,0)
		}else{
			dname[i,"selected"]<-0
		}
		dname[i,"available_core"]<-ifelse(dname[i,"Marker"]%in%horvath_core$predictor,1,0)
		if(dname[i,"Marker"]%in%horvath_core$predictor){
			dname[i,"selected_core"]<-ifelse(horvath_core[horvath_core$predictor==dname[i,"Marker"],"effect_size"]!=0,1,0)
		}else{
			dname[i,"selected_core"]<-0
		}
	}
}

sheet_name<-"DNA Methylation"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,dname)

heading("NMR Metabolomics")

nmr<-read.xlsx("phenotype_dictionary_nmr.xlsx")

clock<-fread("../fixed_alpha/nmr/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/nmr/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")


for(i in 1:nrow(nmr)){
	print(nmr[i,"variable_name"])
	nmr[i,"pass_qc"]<-ifelse(nmr[i,"variable_name"]%in%clock$predictor,1,0)
	if(nmr[i,"variable_name"]%in%clock$predictor){
		nmr[i,"selected"]<-ifelse(clock[clock$predictor==nmr[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		nmr[i,"selected"]<-0
	}
	nmr[i,"available_core"]<-ifelse(nmr[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(nmr[i,"variable_name"]%in%clock_core$predictor){
		nmr[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==nmr[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		nmr[i,"selected_core"]<-0
	}
}

sheet_name<-"NMR Metabolomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,nmr)


heading("Clinomics")

clin<-read.xlsx("clinomics_dictionary.xlsx")

clock<-fread("../fixed_alpha/pheno/fewer/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/pheno/fewer/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")

for(i in 1:nrow(clin)){
	print(clin[i,"variable_name"])
	clin[i,"pass_qc"]<-ifelse(clin[i,"variable_name"]%in%clock$predictor,1,0)
	if(clin[i,"variable_name"]%in%clock$predictor){
		clin[i,"selected"]<-ifelse(clock[clock$predictor==clin[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		clin[i,"selected"]<-0
	}
	clin[i,"available_core"]<-ifelse(clin[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(clin[i,"variable_name"]%in%clock_core$predictor){
		clin[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==clin[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		clin[i,"selected_core"]<-0
	}
}

sheet_name<-"Clinomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,clin)

heading("DEXA")

dexa<-read.xlsx("dexa_dictionary.xlsx")

clock<-fread("../fixed_alpha/dexa/new/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/dexa/new/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")

dexa<-dexa[dexa$variable_name%in%clock$predictor,]

for(i in 1:nrow(dexa)){
	print(dexa[i,"variable_name"])
	dexa[i,"pass_qc"]<-ifelse(dexa[i,"variable_name"]%in%clock$predictor,1,0)
	if(dexa[i,"variable_name"]%in%clock$predictor){
		dexa[i,"selected"]<-ifelse(clock[clock$predictor==dexa[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		dexa[i,"selected"]<-0
	}
	dexa[i,"available_core"]<-ifelse(dexa[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(dexa[i,"variable_name"]%in%clock_core$predictor){
		dexa[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==dexa[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		dexa[i,"selected_core"]<-0
	}
}

sheet_name<-"DEXA"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,dexa)

heading("Metabolon Metabolomics")

met<-read.xlsx("metabolon_metabolomics_dictionary.xlsx")

clock<-fread("../fixed_alpha/metabolon_metabolomics_new_new/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/metabolon_metabolomics_new_new/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")


for(i in 1:nrow(met)){
	print(met[i,"variable_name"])
	met[i,"pass_qc"]<-ifelse(met[i,"variable_name"]%in%clock$predictor,1,0)
	if(met[i,"variable_name"]%in%clock$predictor){
		met[i,"selected"]<-ifelse(clock[clock$predictor==met[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		met[i,"selected"]<-0
	}
	met[i,"available_core"]<-ifelse(met[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(met[i,"variable_name"]%in%clock_core$predictor){
		met[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==met[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		met[i,"selected_core"]<-0
	}
}

sheet_name<-"MS Metabolomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,met)


heading("Metabolon Complex Lipids")

metl<-read.xlsx("metabolon_complex_lipids_dictionary.xlsx")
metl$variable_name<-paste0("X",metl$variable_name)

clock<-fread("../fixed_alpha/metabolon_complex_lipids_new/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/metabolon_complex_lipids_new/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")


for(i in 1:nrow(metl)){
	print(metl[i,"variable_name"])
	metl[i,"pass_qc"]<-ifelse(metl[i,"variable_name"]%in%clock$predictor,1,0)
	if(metl[i,"variable_name"]%in%clock$predictor){
		metl[i,"selected"]<-ifelse(clock[clock$predictor==metl[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		metl[i,"selected"]<-0
	}
	metl[i,"available_core"]<-ifelse(metl[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(metl[i,"variable_name"]%in%clock_core$predictor){
		metl[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==metl[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		metl[i,"selected_core"]<-0
	}
}

sheet_name<-"MS Complex Lipidomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,metl)


heading("MS Fatty Acids Lipidomics")

fa<-read.xlsx("eurospan_fatty_acids_dictionary.xlsx")

clock<-fread("../fixed_alpha/lipidomics/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/lipidomics/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")


for(i in 1:nrow(fa)){
	print(fa[i,"variable_name"])
	fa[i,"pass_qc"]<-ifelse(fa[i,"variable_name"]%in%clock$predictor,1,0)
	if(fa[i,"variable_name"]%in%clock$predictor){
		fa[i,"selected"]<-ifelse(clock[clock$predictor==fa[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		fa[i,"selected"]<-0
	}
	fa[i,"available_core"]<-ifelse(fa[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(fa[i,"variable_name"]%in%clock_core$predictor){
		fa[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==fa[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		fa[i,"selected_core"]<-0
	}
}

sheet_name<-"MS Fatty Acids Lipidomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,fa)


heading("UPLC IgG Glycomics")

igg<-read.xlsx("igg_glycomics_dictionary.xlsx")

clock<-fread("../fixed_alpha/igg_glycomics/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/igg_glycomics/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")


for(i in 1:nrow(igg)){
	print(igg[i,"variable_name"])
	igg[i,"pass_qc"]<-ifelse(igg[i,"variable_name"]%in%clock$predictor,1,0)
	if(igg[i,"variable_name"]%in%clock$predictor){
		igg[i,"selected"]<-ifelse(clock[clock$predictor==igg[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		igg[i,"selected"]<-0
	}
	igg[i,"available_core"]<-ifelse(igg[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(igg[i,"variable_name"]%in%clock_core$predictor){
		igg[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==igg[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		igg[i,"selected_core"]<-0
	}
}

sheet_name<-"UPLC IgG Glycomics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,igg)


heading("Mega Omics")

meg<-read.xlsx("mega_omics_dictionary.xlsx")

clock<-fread("../fixed_alpha/combined_new/p03_make_clock/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock)
dim(clock)

clock<-clock[3:nrow(clock),]
names(clock)<-c("predictor","effect_size")

clock_core<-fread("../fixed_alpha/combined_new/p06_core_model_prediction/st03_model_coefficients_1.tsv",data.table=FALSE)
head(clock_core)
dim(clock_core)

clock_core<-clock_core[3:nrow(clock_core),]
names(clock_core)<-c("predictor","effect_size")


for(i in 1:nrow(meg)){
	print(meg[i,"variable_name"])
	meg[i,"pass_qc"]<-ifelse(meg[i,"variable_name"]%in%clock$predictor,1,0)
	if(meg[i,"variable_name"]%in%clock$predictor){
		meg[i,"selected"]<-ifelse(clock[clock$predictor==meg[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		meg[i,"selected"]<-0
	}
	meg[i,"available_core"]<-ifelse(meg[i,"variable_name"]%in%clock_core$predictor,1,0)
	if(meg[i,"variable_name"]%in%clock_core$predictor){
		meg[i,"selected_core"]<-ifelse(clock_core[clock_core$predictor==meg[i,"variable_name"],"effect_size"]!=0,1,0)
	}else{
		meg[i,"selected_core"]<-0
	}
}


meg$variable_name<-gsub("dexa_","",meg$variable_name)
meg$variable_name<-gsub("igg_","",meg$variable_name)
meg$variable_name<-gsub("lip_","",meg$variable_name)
meg$variable_name<-gsub("mm_","",meg$variable_name)
meg$variable_name<-gsub("mcl_","",meg$variable_name)
meg$variable_name<-gsub("mm_","",meg$variable_name)

sheet_name<-"Mega Omics"
addWorksheet(wb=wb, sheetName=sheet_name)
writeData(wb,sheet=sheet_name,meg)

saveWorkbook(wb, "Biomarkers_All_Assays.xlsx", overwrite = TRUE)




