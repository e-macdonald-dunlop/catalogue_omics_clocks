


make_block_def_icd10<-function(file_in) { 
  icd_list<-readRDS(file_in)
  icd10_code_list_raw<-icd_list$icd10codes
  
  icd10_block_def<-icd10_code_list_raw[which(grepl("Block",icd10_code_list_raw$code)),]
  exception<-icd10_block_def$meaning=="U50 Special purpose ONS coding"
  icd10_block_def$code<-substr(icd10_block_def$code,7,13)
  icd10_block_def[exception,"code"]<-"U50-U50"
  
  names(icd10_block_def)[1] <- "code_block"
  
  
  assault<-which(icd10_block_def$code_block=="X85-Y09")
  icd10_block_def[assault,"code_block"]<- "X85-X99"
  assault_line1<- icd10_block_def[assault,]
  assault_line1[,"code_block"]<- "Y00-Y09"
  icd10_block_def <- rbind(icd10_block_def,assault_line1)
  
  icd10_block_def$letter<-left(icd10_block_def$code)
  icd10_block_def$start<-as.integer(substr(icd10_block_def$code,2,3))
  icd10_block_def$end<-as.integer(right(icd10_block_def$code,2))
  
  icd10_code_list<-icd10_code_list_raw[!grepl("Block",icd10_code_list_raw$code),]
  icd10_code_list$full_code <- parse_isd_icd(icd10_code_list$code)
  icd10_code_list$icd_short<-strip_post_dot(icd10_code_list$full_code)
  icd10_code_list$letter<-left(icd10_code_list$icd_short)
  icd10_code_list$pre_dec_num<-as.integer(right(icd10_code_list$icd_short,2))
  
  icd10_code_pre_dec_list <- icd10_code_list[!duplicated(icd10_code_list$icd_short),4:6]
  
  
  res<-NULL
  for (row in row.names(icd10_code_pre_dec_list)){
    #  for (row in {16140:1700} ){
    row_to_process <- icd10_code_pre_dec_list[row,]
    block_of_int<-icd10_block_def[icd10_block_def$letter==row_to_process$letter,]
    block_of_int<-block_of_int[block_of_int$start <= row_to_process$pre_dec_num,]
    block_of_int<-block_of_int[block_of_int$end >= row_to_process$pre_dec_num,]
    block_of_int <- block_of_int[,c(1,2)]
    ans1<-cbind(row_to_process,block_of_int)
    res<- rbind(res,ans1)
  }
  
  res <- rename_df_cols(res,c("letter","code_block"),c("icd_letter","icd_block"))  
  icd10_block_def <-  rename_df_cols(icd10_block_def,c("letter","code_block"),c("icd_letter","icd_block") )
  return(list(blocks=icd10_block_def, code_to_block=res))
}



#start======================================

options(width=180)
source("../scripts/helper_fns.R")
source(wp("/apps_by_us/gen_useful.R"))



icd10_def<-make_block_def_icd10("../from_mort_morb/icd_list.RDS")

lapply(icd10_def,head)
saveRDS(icd10_def,file="st01_01_icd10_def.RDS")

date()
heading("Done Make Blockdef")

