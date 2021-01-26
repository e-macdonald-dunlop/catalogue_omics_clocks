
wp<-function(file){

wilson_path <- "/exports/igmm/eddie/wilson-lab"
mf <- Sys.info()["nodename"] == "muckleflugga"
if (mf) wilson_path <- "/opt/working/wilson"
  wp<-paste(wilson_path,file,sep="/")
  return(wp)
}



head1 <-function(df,iids="ORCA1487",n=30){

print(head(df[df$iid %in% iids,],n=n))
}

extract_year <- function(dates){
  year<-as.integer(format(dates,"%Y"))
  return(year)
}



empty_string_to_underscore <- function(chars){
  chars<-ifelse(chars=="","_",chars)
  return(chars)
}

get_rid_of_empty_strings_and_spaces<-function(df){
  for (column_name in names(df)){
    df[, column_name]<- empty_string_to_underscore(df[, column_name])
    df[, column_name]<- gsub(" ","_",(df[, column_name]))

  }
  return(df)
}



diff_dates_in_yrs <- function (date1,date2){
  as.numeric((date2-date1)/365.25)
}


tab_smr <-function(smr){
smr$yoa<-extract_year(smr$doa)
smr$meaning<-left(smr$meaning,50)
smr$decade<-floor(smr$yoa/10)*10
print(with(smr,proctabulate(iid,meaning,decade,countf)))
}





left <- function(string,n=1){
  left_bit <- substring(string,1,n)
  return(left_bit)
}

right <- function(string,n=1){
  right_bit <- substring(string,nchar(string)-n+1,nchar(string))
  return(right_bit)
}



strip_post_dot <- function(x){
  y<-gsub("\\..*","",x)
  return(y)
}





read_smr<-function(file_in){
  
  smr0 <- read_in_and_clean(file_in,n=-1)
  smr0$icd9<-parse_isd_icd9(smr0$main_condition)
  smr0$icd10<-parse_isd_icd10(smr0$main_condition)
  
  return(smr0)
}
  





parse_isd_icd <- function(code){
  code <- ifelse(!is.na(code),paste(left(code,3),right(code,nchar(code)-3),sep="."),NA)
  return(code)
  
}



ym_int_to_years_float <- function(ym){
  year<- floor(ym/100) 
  year<- round(year + (ym-year*100)/12-0.5,2)
  return(year)
}



factors_to_char<- function(df){
factors_df<-which(unlist(lapply(df,class))=="factor")

for (col in factors_df){
	df[,col]<-as.character(df[,col])

}
	return(df)
}


get_table<-function(odbc_string,table){
table_query<-paste("select * from", table,sep=" ")
dbhandle <- odbcDriverConnect(odbc_string)
print(table_query)
 res<- sqlQuery(dbhandle, table_query,as.is=T)
odbcClose(dbhandle)

print(head(res))

return(res)
}


rename_df_cols <- function(df,old,new){

for (count in 1:length(old)){
  old1<-old[count]
  new1<-new[count]
names(df)[names(df)==old1]  <-  new1
}
return(df)
}


first_row_nums<-function(df,pick_var){
  ans<-match(unique(df[,pick_var]),df[,pick_var])

}

pick_first <- function(df,order_var,pick_var,decreasing = FALSE){
#picks first (or last) occurnce of each value of pick_var, where order_var defines first
  df<-df[order(df[,order_var],decreasing=decreasing),]

  df_first <- df[first_row_nums(df,pick_var) ,]

  return(df_first)
}


first_of_joint_v1_v2 <- function(df,order,v1,v2){
df$compound_var<- paste(df[,v1],df[,v2],sep="_")
df_first_only<-pick_first(df,order,"compound_var")
df_first_only<-df_first_only[,names(df) != "compound_var"]
return(df_first_only)
}




ym_int_to_r<-function(dt_n){
dt_n<-as.integer(dt_n)*100+15
dt_n <- as.Date(as.character(dt_n),"%Y%m%d")
return(dt_n)
}




writeLines("helper_fns made")
