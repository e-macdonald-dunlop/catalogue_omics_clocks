


source("/exports/igmm/eddie/wilson-lab/apps_by_us/gen_useful.R")
source("../scripts/helper_fns.R")
source("../scripts/helper_fns_plot.R")


file_outcome_map<- "../p01_prep_icd_defs/st01_01_icd10_def.RDS"


icd_def<-readRDS(file_outcome_map)
out_def<-icd_def$blocks
out_def<- out_def[,c(1,2)]
out_def<-rename_df_cols(out_def,"icd_block","outcome")


write_pj1(out_def,file="st03a_01_out_defs.tsv")


