
stdise <- function (x){
	z<-(x-mean(x))/sd(x)
	return(z)
}

sum_and_r2<-function(model){

	print(summary(model))
	return(summary(model)$r.squared)
}


analyse_bv <- function(df0,outcome_name,pred1_name,pred2_name){

df<-df0[,c(outcome_name,pred1_name,pred2_name)]
df<-na.omit(df)	
y<-df[,outcome_name]
x1<-df[,pred1_name]
x2<-df[,pred2_name]
x1<-stdise(x1)
x2<-stdise(x2)
y<-stdise(y)


m1<-lm(y~x1)
m2<-lm(y~x2)
m<-lm(y~x1+x2)

m1_r2<-sum_and_r2(m1)
m2_r2<-sum_and_r2(m2)
m_r2<-sum_and_r2(m)


yhat1<-predict(m1)
yhat2<-predict(m2)
yhat<-predict(m)


m12<-lm(y-yhat1~x2)
m21<-lm(y-yhat2~x1)

m21_r2<-sum_and_r2(m21)
m12_r2<-sum_and_r2(m12)
exp_r2<-1-(1-m1_r2)*(1-m2_r2)
min_r2<-max(m1_r2,m2_r2)
max_r2<-(m1_r2+m2_r2)
xs_overlap<-(exp_r2-m_r2)/(exp_r2-min_r2)
res1<-data.frame(pred1_name,pred2_name,m1_r2,m2_r2,m_r2,m12_r2,m21_r2,min_r2,exp_r2,max_r2,xs_overlap)
return(res1)
}



dendro_order<-function(melt_df0,x,y,value.var){

org_names<-c(x,y,value.var)
melt_df<-melt_df0[,org_names]
names(melt_df)<-c("x","y","value")
wide_df<-dcast(melt_df,x~y,value.var="value")
mat<-wide_df
rownames(mat)<-mat[,1]
mat<-mat[-1]


dendro <- as.dendrogram(hclust(d = dist(x = mat)))
dendro.plot <- ggdendrogram(data = dendro, rotate = TRUE)


order_indexes <- order.dendrogram(dendro)
levs<-wide_df$x[order_indexes]
return(list(levs=levs,plot_obj=dendro.plot))
}


hm_n_dendro <- function(df,x,y,value.var){

dendro_res<-dendro_order(df,x,y,value.var)
levs <- dendro_res$levs

df[,x] <- factor(x = df[,x],
                               levels = levs, 
                               ordered = TRUE)

df[,y] <- factor(x = df[,y],
                               levels = levs, 
                               ordered = TRUE)




gg<-ggplot(data = df, aes_string(x=x, y=y))+
   geom_tile(aes_string(fill=value.var))+
 scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
   midpoint = 0,space = "Lab")  +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1,size=12,hjust = 1),axis.text.y=element_text(size=12))+
 coord_fixed() 


return(list(hm=gg,dendro=dendro_res$plot_obj))
}

mapper<-function(old_name){
  df<-fread("../data/clock_names.txt",data.table=FALSE)
  new_name<-df[df$short==old_name,"long"]
  return(new_name)
}

nice_names<-function(df){
  names_want<-unique(df$pred1_name)
  for(clock in names_want){
    df$pred1_name<-gsub(clock,mapper(clock),df$pred1_name)
    df$pred2_name<-gsub(clock,mapper(clock),df$pred2_name)
  }
  return(df)
}


#=============================







library(data.table)
library(ggplot2)
library(reshape2)
library(ggdendro)
library("grid")
#library(patchwork)
library("gridExtra")


#options(width=200)

source("/exports/igmm/eddie/wilson-lab/apps_by_us/gen_useful.R")

 file_in<- "../data/all_clocks_resid_training_testing_age_09_07_2020.tsv"
 file_out<- "st01_01_bv.tsv"

ages<-read_pj1(file_in)

res<-NULL
pred_names<-names(ages[2:12])
for(pred1_name in pred_names){
	for(pred2_name in pred_names){
		res<-rbind(res,analyse_bv(ages,"age_at_vene",pred1_name,pred2_name))
	}
}


print(res)

heading(file_out)
write_pj1(res,file_out)

heading("done")
date()
