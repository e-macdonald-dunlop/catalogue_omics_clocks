

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
  df<-fread("../data/clock_names_2.txt",data.table=FALSE)
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
 
 file_in<- "../p01_readin_analyse/st01_01_bv.tsv"


res<-read_pj(file_in)

#need to get rid of age at vene from pred1_name and pred2_name
res<-res[res$pred1_name!="age_at_vene" & res$pred2_name!="age_at_vene",]
#need to get nice names for plots
res<-nice_names(res)


pdf("st02_02_plots.pdf",paper="a4r")
plots<-hm_n_dendro(res,"pred1_name","pred2_name","xs_overlap")
print(plots$hm)
print(plots$dendro)
grid.newpage()
print(plots$hm + 
  xlab("") +
  ylab("") +
  guides(fill=guide_legend(title="Excess Overlap")) +
  theme(axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "top"),
  vp = viewport(x = 0.4, y = 0.5, width = .8, height = 1))
print(plots$dendro + theme_dendro() ,
  vp = viewport(x = 0.9, y = 0.54, width = .2, height = 0.55))
dev.off()


heading("done")
date()


