

# From David B
#Wed 29/01/2020 15:50
# Dear all,

 

# You should now have access to the following SMR tables:

 

# dead

# diabetes

# ecg_clean

# pis

# simd

# smr00

# smr01
# smr01c

# smr02

# smr04

# smr06

# spec

# viking_orcades_recruit_dates

 

# --

 

# In the smr01 table patients can have up to six ICD codes assigned at each visit.  smr01c is a flattened version of the table where each code is a separate record.

 

# ICD10 codes generally start with a letter, ICD9 codes with a number.

 

# Some ICD codes have been used in both ICD9 and ICD10 schemes such as V252.

 

# ICD-10 was used from1995 onward so we can use this to work out the correct code.

 

# I have avoided this problem by filtering out the non-disease specific codes of R,S,T,U,V,W,X,Y,Z

 

#======================
source("/exports/igmm/eddie/wilson-lab/apps_by_us/gen_useful.R")
source("../scripts/helper_fns.R")
library(RODBC)
library(reshape2)
library(data.table)

options(width=240)


odbc_string<-'driver=libmsodbcsql-17.4.so.2.1;server=igmm-store.igmm.ed.ac.uk;database=VikingData;UID=vikingdata;Port=1433;PWD=afr87e21r9!'
table <- 'smr01'
#orca_whole_body_data
#orca_weight_results


pheno_file <- "/exports/igmm/eddie/wilson-lab/data/base_data/orcades/phenotypes/orcades_base_phenotypes.tsv"
icd_10_def_rds_file<-"../p01_prep_icd_defs/st01_01_icd10_def.RDS"
docodes_to_icd10 <- ym_int_to_r(199604)
resids_file<-"/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/between_clock_correlations/all_clocks_resid.tsv"

wanted_ph_cols<-c("iid","sex","dob","dovene","height","weight")

poss_dates <- c("vene_date" ,     "meas_date"  ,    "eye_date" ,      "bone_date"  ,    "memory_date"   )
dolast_smr<-ym_int_to_r(201712)


#date convention is that dates are R dates, unless suffixed by _ym_int in which case per above
#cant remmeber how figured this out - maybe length of icd code - it is the year smr moved to icd10


icd10_def <- readRDS(icd_10_def_rds_file)



smr_wide<-get_table(odbc_string,table)
class(smr_wide$dt1)<-"integer";class(smr_wide$dt2)<-"integer"

smr_wide <- rename_df_cols (smr_wide,"id","iid")


smr_wide <-  get_rid_of_empty_strings_and_spaces(smr_wide)


write_pj1(smr_wide,"st02_01_smr01.txt")


heading ("smr written to st02_01_smr01.txt ")

print(head(smr_wide))

date()
heading("Done")


q()







setwd("../p03_process_smr")

smr_wide <- read_pj("../p02_read_smr/st02_01_smr01.txt")
smr <- melt(smr_wide,id.vars=c("iid","table_id","dt1"),measure.vars=c("c1","c2","c3","c4","c5","c6"),value.name="icd")

smr <- rename_df_cols(smr,c("variable","dt1","icd"),c("cond_num","doa_ym_int","icd_is"))
smr$cond_num <- as.integer(gsub("c","",smr$cond_num))
smr<-smr[smr$icd_is !="",]

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
smr <- smr[!is.na(smr$icd_block),]


heading("ICD10 Admissions")
heading("subject counts
	")
length(unique(smr$iid))
tab_smr(smr)

#----
# First admission only
#----


smr$icd_block <- smr$icd_letter

smr_first <- first_of_joint_v1_v2 (smr,"doa","iid","icd_block")






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
proctabulate(iid,1,sex,countf)



smr1<-smr_first

smr1<- smr1[,c("iid","icd_block","doa")]
 table(smr1$icd_block)


head(smr1)
#following is like a cross multiply mayebe theres a smarter way, bbut want a iid * icd_code for every iid/code
 framework_table <- merge(pheno1,unique(smr1$icd_block))
framework_table <- rename_df_cols(framework_table,"y","icd_block")

smr2<- merge(framework_table,smr1,all.x=T)


smr2$prevalent <- ifelse(!is.na(smr2$doa) & (smr2$doa < smr2$dovene),T,F)


smr2$do_end <- with(smr2,fifelse(is.na(doa),dolast_smr, doa))
smr2$incid_event <- with(smr2, ifelse(is.na(doa) & !smr2$prevalent,0,1))
smr2$age_start <- with(smr2,diff_dates_in_yrs(dob,dovene))
smr2$age_end <- with(smr2,diff_dates_in_yrs(dob,do_end))



smr3<-smr2[smr2$icd_block=="I" & !smr2$prevalent,]

write_pj1(smr3,"st02_02_incid_I20_25.txt")


library(survival)
 
mod1<- coxph(Surv(age_start,age_end,incid_event)~sex+height+weight,data=smr3)

head(smr2)

pdf("st02_02_doa.pdf",paper="a4r")


dev.off()

heading("Done a02_01_make_smr")




