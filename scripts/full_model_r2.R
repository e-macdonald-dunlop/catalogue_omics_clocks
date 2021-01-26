library(data.table)
library(psych) # for descriptive statistics
library(ppcor) # this pacakge computes partial and semipartial correlations.
library(ggplot2)
library(corrplot)

args<-commandArgs(T)
#args<-c("/exports/igmm/eddie/wilson-lab/projects/prj_086_omics_clocks/final/overlap/",FALSE)

wd<-args[1]

if(!file.exists(wd)){dir.create(wd)}

setwd(wd)
source("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/pipeline_functions.R")

heading("Read in Data")
if(!file.exists("./a05_figure")){dir.create("./a05_figure")}
#read in full results
data<-fread("./a04_semi_partial_correlations/full_model_sp_correlations_merged.tsv",data.table=FALSE)
head(data)
dim(data)

full_r2<-fread("./a02_full_model/full_model_r2.tsv",data.table=FALSE)
head(full_r2)

unexplained<-1-full_r2$x
u_sum<-sum(data$sr2)
o_lap<-1-unexplained-u_sum

df<-data[,c("predictor","sr2")]
df2<-data.frame(predictor=c("unexplained","overlap"),sr2=c(unexplained,o_lap))

df<-rbind(df,df2)


##################
# Need to read in all resids
##################

resids<-fread("./a01_gather_resid/all_clocks_resid_training_testing_age_09_07_2020.tsv",data.table=FALSE)
clocks<-names(resids)[!names(resids)%in%c("iid","age_at_vene")]

new<-lapply(clocks,function (x){
  model<-lm(as.formula(paste0("age_at_vene ~",x)),data=resids)
  r2<-summary(model)$r.squared
  return(data.frame(clock=x,r2=r2))
  })

new<-do.call(rbind,new)
names(new)<-c("predictor","r2")

new$predictor<-gsub("nmr","NMR Metabolomics",new$predictor)
new$predictor<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",new$predictor)
new$predictor<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",new$predictor)
new$predictor<-gsub("igg_glycomics","UPLC IgG Glycomics",new$predictor)
new$predictor<-gsub("protein_new","PEA Proteomics",new$predictor)
new$predictor<-gsub("dexa","DEXA",new$predictor)
new$predictor<-gsub("pheno","Clinomics",new$predictor)
new$predictor<-gsub("lipidomics","MS Fatty Acids Lipidomics",new$predictor)
new$predictor<-gsub("horvath_cpgs","DNAme Horvath CpGs",new$predictor)
new$predictor<-gsub("hannum_cpgs","DNAme Hannum CpGs",new$predictor)
new$predictor<-gsub("combined","Mega Omics",new$predictor)
new$predictor<-gsub("overlap","Overlap",new$predictor)
new$predictor<-gsub("unexplained","Unexplained",new$predictor)


#doughnut plot testing
df$ymax = cumsum(df$sr2)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2


df$predictor<-gsub("nmr","NMR Metabolomics",df$predictor)
df$predictor<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",df$predictor)
df$predictor<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",df$predictor)
df$predictor<-gsub("igg_glycomics","UPLC IgG Glycomics",df$predictor)
df$predictor<-gsub("protein_new","PEA Proteomics",df$predictor)
df$predictor<-gsub("dexa","DEXA",df$predictor)
df$predictor<-gsub("pheno","Clinomics",df$predictor)
df$predictor<-gsub("lipidomics","MS Fatty Acids Lipidomics",df$predictor)
df$predictor<-gsub("horvath_cpgs","DNAme Horvath CpGs",df$predictor)
df$predictor<-gsub("hannum_cpgs","DNAme Hannum CpGs",df$predictor)
df$predictor<-gsub("combined","Mega Omics",df$predictor)
df$predictor<-gsub("overlap","Overlap",df$predictor)
df$predictor<-gsub("unexplained","Unexplained",df$predictor)

df$predictor<-factor(df$predictor,levels=c("Overlap","NMR Metabolomics","MS Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","PEA Proteomics","MS Fatty Acids Lipidomics","DNAme Hannum CpGs","DNAme Horvath CpGs","DEXA","Clinomics","Unexplained")) #"Mega Omics",

#want to write out to file
out<-merge(df[,c("predictor","sr2")],new,by="predictor",all.x=TRUE)
write.table(out,"./a05_figure/unique_ve_table.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

heading("Creating Doughnut Plot")
ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=predictor)) +
     geom_rect() +
     #geom_text( x=2, aes(y=labelPosition, label=predictor, color=predictor), size=6) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     scale_fill_manual(values=c("#8c6bb1","#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4","#bdbdbd")) + #"#ffffcc","#bcbddc",
     xlim(c(1, 4)) + # Try to remove that to see how to make a pie chart
     theme_void()
ggsave("./a05_figure/doughnut.pdf")

p<-ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=1, fill=predictor)) +
     geom_rect() +
     #geom_text( x=2, aes(y=labelPosition, label=predictor, color=predictor), size=6) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     scale_fill_manual(values=c("#8c6bb1","#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4","#bdbdbd")) + #"#ffffcc","#bcbddc",
     #xlim(c(1, 4)) + # Try to remove that to see how to make a pie chart
     theme_void()
print(p)
ggsave("./a05_figure/pie.pdf")

heading("Creating Barchart")

df<-data[,c("predictor","sr2")]
df2<-data.frame(predictor=c("unexplained","overlap"),sr2=c(unexplained,o_lap))

df<-rbind(df,df2)
df<-df[!grepl("overlap|unexplained",df$predictor),]

df$predictor<-gsub("nmr","NMR Metabolomics",df$predictor)
df$predictor<-gsub("metabolon_metabolomics_new_new","MS Metabolomics",df$predictor)
df$predictor<-gsub("metabolon_complex_lipids_new","MS Complex Lipidomics",df$predictor)
df$predictor<-gsub("igg_glycomics","UPLC IgG Glycomics",df$predictor)
df$predictor<-gsub("protein_new","PEA Proteomics",df$predictor)
df$predictor<-gsub("dexa","DEXA",df$predictor)
df$predictor<-gsub("pheno","Clinomics",df$predictor)
df$predictor<-gsub("lipidomics","MS Fatty Acids Lipidomics",df$predictor)
df$predictor<-gsub("horvath_cpgs","DNAme Horvath CpGs",df$predictor)
df$predictor<-gsub("hannum_cpgs","DNAme Hannum CpGs",df$predictor)
df$predictor<-gsub("combined","Mega Omics",df$predictor)
df$predictor<-gsub("overlap","Overlap",df$predictor)
df$predictor<-gsub("unexplained","Unexplained",df$predictor)

df$predictor<-factor(df$predictor,levels=c("NMR Metabolomics","MS Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","PEA Proteomics","MS Fatty Acids Lipidomics","DNAme Hannum CpGs","DNAme Horvath CpGs","DEXA","Clinomics")) #"Mega Omics",


empty_bar <- 10
to_add <- matrix(NA, empty_bar, ncol(df))
colnames(to_add) <- colnames(df)
df <- rbind(df, to_add)
df$id <- seq(1, nrow(df))

ggplot(df,aes(x=as.factor(predictor),y=sr2,fill=predictor)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4")) + #"#ffffcc","#bcbddc",
    ylim(-0.001,0.01) +
    theme_minimal() +
  	theme(
    	axis.text = element_blank(),
    	axis.title = element_blank(),
    	panel.grid = element_blank(),
    	plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  	) +
  	coord_polar(start = 0)
ggsave("./a05_figure/bar.pdf")


#need to sort the order
#sort_order<-df[order(df$sr2),"predictor"]
#df$predictor<-factor(df$predictor,levels=sort_order)
#c("NMR Metabolomics","Metabolon Metabolomics","Metabolon Complex Lipids","IgG Glycomics","Proteins","Lipidomics","DNA Methylation Hannum CpGs","DNA Methylation Horvath CpGs","DEXA new","DEXA u","Clinomics"))
#c("#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4")
df<-df[!is.na(df$predictor),]

df$predictor<-factor(df$predictor,levels=c("NMR Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","MS Fatty Acids Lipidomics","DEXA","MS Metabolomics","Clinomics","DNAme Horvath CpGs","DNAme Hannum CpGs","PEA Proteomics"))
#"#dd3497","#bfd3e6","#02818a","#ece7f2","#4292c6","#c7e9b4","#66c2a4","#993404","#fe9929","#7a0177"
b<-ggplot(df,aes(x=predictor,y=sr2,fill=predictor)) +
  geom_bar(stat="identity") +
  #scale_fill_manual(values=c("#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4")) + #"#ffffcc","#bcbddc",
  scale_fill_manual(values=c("#02818a","#dd3497","#bfd3e6","#ece7f2","#4292c6","#c7e9b4","#66c2a4","#993404","#fe9929","#7a0177")) + #"#ffffcc","#bcbddc",
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_blank(),legend.position = "none")
print(b)
ggsave("./a05_figure/bar_straight.pdf",width=3,height=2)

#try to plot them together
#library(patchwork)

#p + theme(legend.position = "bottom") + b + theme(legend.position = "bottom")

library(ggpubr)

ggarrange(p + theme(legend.text = element_text(size=8),legend.title = element_blank()),b, labels=c("A","B"),nrow = 1,common.legend = TRUE, legend="bottom") #align = "h", + theme(legend.text = element_text(size=8),legend.title = element_blank())
ggsave("./a05_figure/combined.pdf",width=6,height=4)
ggsave("./a05_figure/combined.png",width=6,height=4)

library(patchwork)

combined <- b + p & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("./a05_figure/patchwork.pdf",width=6,height=4)

library(cowplot)
prow <- plot_grid( b + theme(legend.position="none"),
           p + theme(legend.position="none"),
           align = 'vh',
           labels = c("A", "B"),
           hjust = -1,
           nrow = 1
           )

legend_b <- get_legend(p + theme(legend.position="bottom"))
x <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
x
ggsave("./a05_figure/cowplot.pdf",width=6,height=4)


#library(gridExtra)
#grid.arrange(p, b,ncol = 2, nrow = 1)

#library(RColorBrewer)


#myPalette <- brewer.pal(12, "Paired") 
#pie(df$sr2,labels=df$predictor,border="white",col=myPalette)


#ggplot(df, aes(x="", y=sr2, fill=predictor)) +
  #geom_bar(stat="identity", width=1, color="white") +
  #coord_polar("y", start=0) +
  #theme_void() 

#ggplot(df, aes(x=predictor, y=sr2, fill=predictor)) +
  #geom_bar(stat="identity", width=1, color="white")
  

new_data<-merge(df,new,by="predictor",all.x=TRUE)
long_data<-reshape2::melt(new_data)

long_data$variable<-factor(long_data$variable,levels=c("r2","sr2"))

new_names <- c(`sr2` = "Unique Variance Explained",`r2` = "R Squared")


new_data$r2_new<-new_data$r2/100
new_data$predictor<-factor(new_data$predictor,levels=c("NMR Metabolomics","MS Complex Lipidomics","UPLC IgG Glycomics","MS Fatty Acids Lipidomics","DEXA","MS Metabolomics","Clinomics","DNAme Horvath CpGs","DNAme Hannum CpGs","PEA Proteomics"))


b<-ggplot(new_data,aes(x=predictor,fill=predictor)) +
  geom_bar(aes(x=predictor,y=sr2),position="dodge",stat="identity") +
  #facet_wrap(~variable,scales="free",ncol=1, labeller = as_labeller(new_names)) +
  geom_point(aes(x=predictor,y=r2_new)) +
  #scale_fill_manual(values=c("#02818a","#c7e9b4","#dd3497","#bfd3e6","#7a0177","#ece7f2","#fe9929","#993404","#4292c6","#66c2a4")) + #"#ffffcc","#bcbddc",
  #scale_x_discrete()+
  #scale_fill_manual(values=c("#bfd3e6","#02818a","#dd3497","#ece7f2","#993404","#4292c6","#66c2a4","#7a0177","#c7e9b4","#fe9929")) + #"#ffffcc","#bcbddc",
  scale_fill_manual(values=c("#02818a","#dd3497","#bfd3e6","#ece7f2","#4292c6","#c7e9b4","#66c2a4","#993404","#fe9929","#7a0177")) + #"#ffffcc","#bcbddc",
  scale_y_continuous(name = "sr2" , limits=c(0,0.01),sec.axis = sec_axis(trans=~.*100,name="r2")) + #
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_blank(),legend.position = "none") #,strip.background = element_blank(),strip.text.x = element_blank()
print(b)


ggarrange(p + theme(legend.text = element_text(size=8),legend.title = element_blank()),b, labels=c("A","B"),nrow = 1,common.legend = TRUE, legend="bottom") #align = "h", + theme(legend.text = element_text(size=8),legend.title = element_blank())
ggsave("./a05_figure/combined.pdf",width=6,height=4)
ggsave("./a05_figure/combined.png",width=6,height=4)

