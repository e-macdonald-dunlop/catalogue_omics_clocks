

#======================
source("../scripts/helper_fns.R")
source(wp("/apps_by_us/gen_useful.R"))

library(RODBC)
library(reshape2)
library(data.table)

options(width=240)

docodes_to_icd10 <- ym_int_to_r(199604)
dolast_smr<-ym_int_to_r(201712)
#imp_disease_chapters <- c("C","E","F","I","J")
imp_disease_chapters <- c("C","E","I","J")





pheno_file <- wp("/data/base_data/orcades/phenotypes/orcades_base_phenotypes.tsv")
icd_10_def_rds_file<-"../p01_prep_icd_defs/st01_01_icd10_def.RDS"

#resids_file<-wp("/projects/prj_086_omics_clocks/final/between_clock_correlations/all_clocks_resid.tsv")

wanted_ph_cols<-c("iid","sex","dob","dovene","height","weight")
poss_dates <- c("vene_date" ,     "meas_date"  ,    "eye_date" ,      "bone_date"  ,    "memory_date"   )


#date convention is that dates are R dates, unless suffixed by _ym_int in which case per above
#cant remmeber how figured this out - maybe length of icd code - it is the year smr moved to icd10


icd10_def <- readRDS(icd_10_def_rds_file)
smr_wide <- read_pj("../p02_read_smr/st02_01_smr01.txt")



smr <- reshape2::melt(smr_wide,id.vars=c("iid","table_id","dt1"),measure.vars=c("c1","c2","c3","c4","c5","c6"),value.name="icd")

smr <- rename_df_cols(smr,c("variable","dt1","icd"),c("cond_num","doa_ym_int","icd_is"))
smr$cond_num <- as.integer(gsub("c","",smr$cond_num))
smr<-smr[smr$icd_is !="_",]

smr<-rename_df_cols(smr,"doa_ym_int","doa")
smr$doa <- ym_int_to_r(smr$doa)




smr$icd_type<-ifelse(smr$doa < docodes_to_icd10 ,9,10)

                                                              

smr$icd<-parse_isd_icd(smr$icd_is)

smr <- smr[smr$icd_type==10,]
smr <-smr[grepl("ORCA",smr$iid),]

smr$icd_short <- gsub("\\..*","",smr$icd)



icd10_code_to_block_map_df <- icd10_def$code_to_block

smr <- merge(smr,icd10_code_to_block_map_df,all.x=T)

#k64.9 - unspecified hemorrhoids somehow not mapping to table in 14 rows?? should check - just delete for now
#there are also some nasty xxx_x conditions
smr <- smr[!is.na(smr$icd_block),]


heading("ICD10 Admissions")
heading("subject counts
	")
length(unique(smr$iid))
tab_smr(smr)

heading("Selected ICD10 Admissions")
smr<- smr[smr$icd_letter %in% imp_disease_chapters,]
tab_smr(smr)




#----
# First admission only
#----

smr_first<-NULL

for (group_level in c("all","chapter","block")){
	 heading(paste("group_level",group_level))
	 smr_t<-smr
	if (group_level=="all") smr_t$group<-"all"
	if (group_level=="chapter") smr_t$group<-smr$icd_letter
	if (group_level=="block") smr_t$group<-smr$icd_block
	smr_t$group_level<-group_level
	smr_first1 <- first_of_joint_v1_v2 (smr_t,"doa","iid","group")
	print(with(smr_first1,proctabulate(iid,group,1,countf)))
	smr_first<-rbind(smr_first,smr_first1)
}
smr_first$group <- make.names(smr_first$group)
cols_wanted<-c( "iid","doa","cond_num","icd","group","group_level", "meaning"   )
smr_first <- smr_first[,cols_wanted]

#next bit is a wee bit clever - because the group values are distinct at diffecernt levels, we should get first of each group_level
smr_first <- first_of_joint_v1_v2 (smr_first,"doa","iid","group")

print(with(smr_first,proctabulate(iid,group,1,countf)))




heading("ICD10 first admissions - ORCADES")
tab_smr(smr_first)



pheno<-read.table(pheno_file,header=T,as.is=T,sep="\t")
pheno<-pheno[!grepl("NIMS",pheno$iid),]
dim(pheno)




pheno1 <- pheno
pheno1 <-pheno1[grepl("ORCA",pheno$iid),]

pheno1<-rename_df_cols(pheno1,c("date_of_birth"),c("dob"))


pheno1$dovene <- as.Date(NA)
#dovene missing too often opverwrite with a plausible non missing note use of data.table::fifelse to preserve date class
for (poss_date in poss_dates){
pheno1$dovene <- fifelse(is.na(pheno1$dovene),as.Date(pheno1[,poss_date],"%d/%m/%Y"),pheno1$dovene)
}
pheno1<-pheno1[!is.na(pheno1$dovene),]

pheno1$dob <- as.Date(pheno1$dob,"%d/%m/%Y")
pheno1[pheno1$vene_date=="",c(poss_dates,"dovene","dob")]

pheno1<-pheno1[,wanted_ph_cols]


heading("here's the sample")

dim(pheno1)
with(pheno1,proctabulate(iid,1,sex,countf))



smr1<-smr_first

smr1<- smr1[,c("iid","group","doa")]
 table(smr1$group)


head(smr1)
#following is like a cross multiply mayebe theres a smarter way, bbut want a iid * icd_code for every iid/code
 framework_table <- merge(pheno1,unique(smr1$group))
framework_table <- rename_df_cols(framework_table,"y","group")

smr2<- merge(framework_table,smr1,all.x=T)

#I variables are T or NA
smr2$prevalent <- ifelse(!is.na(smr2$doa) & (smr2$doa < smr2$dovene),T,F)
smr2$not_prevalent_i <- ifelse(!is.na(smr2$doa) & (smr2$doa < smr2$dovene),NA,T)

smr2$do_end <- with(smr2,fifelse(is.na(doa),dolast_smr, doa))
smr2$incid_event <- with(smr2, ifelse(is.na(doa) & !smr2$prevalent,0,1))
smr2$age_start <- with(smr2,diff_dates_in_yrs(dob,dovene))
smr2$age_end <- with(smr2,diff_dates_in_yrs(dob,do_end))
smr2$t<-smr2$age_end-smr2$age_start

smr2<-smr2[,c(1,2,3,12,15,13,14,4,5,11,10,6,7,8,9)]

#now make a variable which says the analysis group eg gp_all, gp_C, gp_I20-25, .... if we add + gp_I20-25 to our anaysis plan will pick out only these rows :)

for (group in unique(smr2$group)){
	is_from_group_var<-paste("gp",group,"i",sep="_")
smr2[,is_from_group_var]<-ifelse(smr2$group==group,T,NA)
}

names(smr2)<-tolower(names(smr2))



write_pj1(smr2,"st02_02_smr.txt")
write_pj1(smr2[!is.na(smr2$gp_all_i),],"st02_02_smr_all.txt")

incid_vs_prev <- with(smr2,proctabulate(iid,group,paste(incid_event,not_prevalent_i),countf))
colnames(incid_vs_prev)<- c("free","prevalent","incident","total")

incid_vs_prev<-as.data.frame(incid_vs_prev)

incid_vs_prev$icd_block<-row.names(incid_vs_prev)
incid_vs_prev$icd_block<-gsub("\\.","-",incid_vs_prev$icd_block)
incid_vs_prev<-incid_vs_prev[!incid_vs_prev$icd_block=="group",]
#h25-28 has 4 rows H25, H26, ..
block_def<-unique(icd10_code_to_block_map_df[,c("icd_block","meaning")])
incid_vs_prev<-merge(incid_vs_prev,block_def,all.x=T)
print(incid_vs_prev)


analysis_plan<-data.frame(anal_id='standard_covariates',	formula="not_prevalent_i + sex",	analysis_type=as.character(NA),	residual_transform=as.character(NA),fit_mixed_model=as.character(NA),stringsAsFactors =F)
analysis_fla_template<-"Surv(t,incid_event)~ gp_XXXX_i+ age_start + standard_covariates"
for (group in tolower(unique(smr2$group)) ) {
#for (group in "all"){
analysis_formula <- gsub("XXXX",group,analysis_fla_template)
analysis_plan1<-data.frame(anal_id=group,	formula=analysis_formula,	analysis_type="survival",	residual_transform="none",fit_mixed_model=FALSE,stringsAsFactors = F)

analysis_plan<-rbind(analysis_plan,analysis_plan1)
}

write.table(analysis_plan,"analysis_plan.txt",sep="\t",row.names=F)

analysis_plan$formula <-  gsub("\\+ age_start","",analysis_plan$formula)

write.table(analysis_plan,"analysis_plan_wout_age.txt",sep="\t",row.names=F)



q()

#for (group in unique(smr2$group)){


analysis_plan<-data.frame(anal_id='standard_covariates',	formula="not_prevalent_I + sex",	analysis_type=as.character(NA),	residual_transform=as.character(NA),fit_mixed_model=as.character(NA),stringsAsFactors =F)
analysis_fla_template<-"Surv(t,incid_event)~ gp_XXXX_I+ age_start + standard_covariates"
for (group in "all"){
analysis_formula <- gsub("XXXX",group,analysis_fla_template)
analysis_plan1<-data.frame(anal_id=group,	formula=analysis_formula,	analysis_type="survival",	residual_transform="none",fit_mixed_model=FALSE,stringsAsFactors = F)

analysis_plan<-rbind(analysis_plan,analysis_plan1)
}

write.table(analysis_plan,"analysis_plan_all.txt",sep="\t",row.names=F)

date()
heading("Done a03_01_process_smr")









#write_yaml(analysis_plan,"analysis_plan.yml")

mod1<- coxph(Surv(age_start,age_end,incid_event)~gp_C_I+not_prevalent_I+ sex,data=smr3)


ln -s ../../p03_process_smr/st02_02_smr.txt st02_01_qcd_phenotypes.txt 



library(survival)
 
mod1<- coxph(Surv(age_start,age_end,incid_event)~sex,data=smr3)
summary(mod1)

mod2<- coxph(Surv(t,incid_event)~age_start+sex,data=smr3)
summary(mod2)

pdf("st02_02_doa.pdf",paper="a4r")


dev.off()
date()
heading("Done a03_01_process_smr")


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
require('survival')
require('survminer')


library(ggfortify)
library(survival)

mod1<- survfit(Surv(t,incid_event)~sex,data=smr3)
summary(mod1)
autoplot(mod1)


