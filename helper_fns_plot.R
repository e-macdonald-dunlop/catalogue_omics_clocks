

remove_df_cols<-function(df,cols_to_remove){


df <- df[,!colnames(df) %in% cols_to_remove]
return(df)
}

passed_fdr <- function(df,thresh){
  
  df$q <- p.adjust(df$p)

  df1<-df[df$q<thresh,]
  return(df1)
}


block_n_norm<-function(n,prop_block,mu1,mu2,sd2){

heading("Gnerating a block+normal distn")
writeLines(paste("prop block",prop_block,"at",mu1))
writeLines(paste("prop normal",1-prop_block,"at",mu2,"sd",sd2) )  
d1<-rep(mu1,n*prop_block)
d2<-rnorm(n*(1-prop_block),mu2,sd2)

d3<-c(d1,d2)

return(d3)
}





#https://stats.stackexchange.com/questions/13306/how-can-i-compute-a-posterior-density-estimate-from-a-prior-and-likelihood

bayes_est_beta_1<-function(prior,obs1_beta,obs1_se,plots=TRUE){

#prior is for beta, obs1erved beta and se


likfun <- function(beta,obs1_beta,obs1_se) {
 sapply( beta, function(t) prod( dnorm(obs1_beta, t, obs1_se) ) )
 
} 



#xmin<-min(prior)
#xmax<-max(prior)
#range<-xmax-xmin
#x<-seq(xmin,xmax,range/100)
#lik <- likfun(x,0,1)
#plot(x,lik)

tmp<-likfun(prior,obs1_beta,obs1_se)
post <- sample( prior, 100000, replace=TRUE, prob=tmp )

if(plots){
xmin<-min(post,obs1_beta-2*obs1_se)
xmax<-max(post,obs1_beta+2*obs1_se)

hist(post,xlim=c(xmin,xmax))
abline(v=obs1_beta,col="green")
abline(v=obs1_beta+2*obs1_se,col="green",lty="dashed")
abline(v=obs1_beta-2*obs1_se,col="green",lty="dashed")

abline(v=mean(post),col="blue")
abline(v=mean(post)+2*sd(post),col="blue",lty="dashed")
abline(v=mean(post)-2*sd(post),col="blue",lty="dashed")
}

return(post)

}


bayes_est_beta<-function(prior,obs_beta,obs_se,summary=TRUE,plots=TRUE){

post<-NULL
for (count in {1:length(obs_beta)}){


post1<-bayes_est_beta_1(prior,obs_beta[count],obs_se[count],plots=plots)

if(summary) {beta=mean(post1)
    se=sd(post1)
  post1<-data.frame(obs_beta=obs_beta[count],obs_se=obs_se[count],beta,se)  
  } else {
  post1<-data.frame(obs_beta=obs_beta[count],obs_se=obs_se[count],t(post1))
  } 

post<-rbind(post,post1)
}

return(post)
}



bayes_est_beta_simple_prior<-function(prior_beta,prior_se,   obs_beta,obs_se){




prior_iv <- 1/prior_se^2
obs_iv <- 1/obs_se^2

iv<-prior_iv+obs_iv

beta<- (prior_iv * prior_beta + obs_iv * obs_beta)/iv
se <- 1/sqrt(iv)

post<- data.frame(obs_beta,obs_se,beta,se)


return(post)
}







dev.on<- function(file_out,format,height=8.27,width=11.69){



  if(format=="pdf") pdf(paste(file_out,format,sep="."),height=height,width=width)

  if(format=="jpg") jpeg(paste(file_out,format,sep="."),height=800,width=1200)
return()
}

kk_plot<-function(df,x="x",attribute="attribute",xlab="SNP",ylab="Beta (95% CI)",vline=0){
 

df$x<-df[,x]
df$attribute<-df[,attribute]

df$lower<-df$beta-2*df$se
df$upper<-df$beta+2*df$se

fp <- ggplot(data=df, aes(x=x, y=beta, ymin=lower, ymax=upper,colour=attribute)) +
        geom_pointrange(position=position_dodge(.05), stat = "identity") + 
        geom_hline(yintercept=vline, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab(xlab) + ylab(ylab) +
        theme_bw()  # use a white background


return(fp)



}



add_rownames_as_first_col<-function(df){
  df<-cbind(row.names(df),df)
  return(df)
}






library(tidyverse)
 library(stringr)


split_into_multiple <- function(df,col_name, pattern = ", ", new_col_names){
  cols <- str_split_fixed(df[,col_name], pattern, n = Inf)
  # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
  cols[which(cols == "")] <- NA
  cols <- as.data.frame(cols, stringsAsFactors=F)
  # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m' 
  # where m = # columns of 'cols'
  

  names(cols) <- new_col_names
  old_col_num<-which(names(df)==col_name)
  df<-cbind(df[,1:(old_col_num-1),drop=F],cols,df[,(old_col_num+1):ncol(df),drop=F] )

  return(df)

}


read_in_assoc<-function(file,one_sided=F,quant_trait=F,clock_names,traits_to_reverse,traits_to_exc,age_only=F) {

if(is.na(file)) return(NULL)

a<-read.table(file,header=T,as.is=T)
heading(paste("Reading in",file))
print(head(a))


#emds messy names
for (row_num in 1:nrow(clock_names)){
  one_line<- clock_names[row_num,]
  a$score <- gsub(one_line$syst_name,one_line$clock,a$score)
}

#NB sub i.e. first _ only
a$score<-sub("_","-",a$score)
b<-split_into_multiple(a,"score","-",c("clock","type"))


if (any(grepl("_p0",b$type))) {
  b<-(split_into_multiple(b,"type","_p0",c("method","features")))
  b$features <- ifelse(b$features=="3_make_clock","full","core")
   } else {
  b<-cbind(b[,1:3],"full",b[4:ncol(b)])
  names(b)[3:4] <- c("method","features")
}

b$method<-sub("fixed_alpha_","fixed_alpha-",b$method)
b$method<-sub("at_vene","at_vene-",b$method)
b<-split_into_multiple(b,"method","-",c("method","submethod"))
b$submethod <-ifelse(b$submethod==b$method,"all",b$submethod)
b$submethod <-ifelse(is.na(b$submethod),"all",b$submethod)


res<-b[!is.na(b$p)  ,]


# lose the clock regression type - we're just using one


res<-rename_df_cols(res,c("anal_id","score"),c("outcome","clock"))








#res$q <- p.adjust(res$p)
#better to calculate q on the fly - is affected by subsetting
res$outcome<-gsub("\\.","-",res$outcome)
res$outcome<-toupper(res$outcome)
res$outcome<-left(res$outcome,7)

if(quant_trait) res <- standardise_assoc(res,vars="trait")


res<-res[!(res$trait %in% traits_to_exc),]


if (age_only) res<-res[res$clock =="age",] else res<-res[res$clock !="age",]

reversible_cols<-c("beta","z","trait_mean")
if (!is.null(traits_to_reverse)) res[res$trait %in% traits_to_reverse,reversible_cols] <- -res[res$trait %in% traits_to_reverse,reversible_cols]


if (one_sided) res$p<-pnorm(-res$z) else res$p <- 2* pnorm(-abs(res$z))

return(res)
}



pj_one_sided_sign_test_not_vec<- function(x,y){
p<-	binom.test(x, y,alternative="greater")$p.value
return(p)
}

pj_one_sided_sign_test_vec <- function (xs,ys){
#xs are fails
ps<-mapply(pj_one_sided_sign_test_not_vec,xs,ys)

return(ps)
}


sign_and_show <- function (df,split){


 
 sign_test_0 <-as.data.frame(table(df[,split],df$beta>0))
 sign_test <-dcast(sign_test_0,Var1~Var2)
 names(sign_test) <- c(split,"fail","success")
 sign_test$n <-  sign_test$success + sign_test$fail 

  sign_test$p  <- pj_one_sided_sign_test_vec(sign_test$success,sign_test$n)
return(sign_test)
}


pick_some_methods <- function(df,list_of_methods){
  var_ord <- names(df)
  df1<-merge(df,list_of_methods)
  df1<-df1[,var_ord]
  return(df1)
}



display_matching_cond<-function(df,cond,heading_str){
heading(heading_str)

print(proctabulate(df$outcome,cond,1,countf))

df<-df[cond,]
print(df)


meaning_outs<-paste(df$meaning_out,collapse=",")
return(meaning_outs)
}


ivm_mu <- function(beta,se){
  library(rmeta)
  mu <-  meta.summaries(beta,se)$summary
  return(mu)
}

ivm_se <- function(beta,se){
  library(rmeta)
  mu <-  meta.summaries(beta,se)$se.summary
  return(mu)
}



ivm_t_test <- function(df,gp){

mus <- ddply(df, gp, summarise, beta=ivm_mu(beta,se))

ses <-ddply(df, gp, summarise, se=ivm_se(beta,se))

res<-merge(mus,ses)

all_row <- data.frame(a=gp,beta=with(res,ivm_mu(beta,se)),se=with(res,ivm_se(beta,se)))
names(all_row)[1]<- gp

res <- rbind(res,all_row)
res$t<-res$beta/res$se
res$p<- pnorm(-abs(res$t)) *2
return(res)
}

write_pj1<-
function (df,file_out,sep="\t"){

df <- format(df,digits=3)
  write.table(df,file_out,col.names=T,row.names=F,sep=sep)

return()
}


show_beta_n_se_tables <- function(df){



clock_res_ave_outcome<- ivm_t_test(df,"clock")
out_res_ave_clock<- ivm_t_test(df,"outcome")

betas<-with(df,proctabulate(beta,outcome,clock,meanf))
#each cell is only one entry so mean is that single value, except the
#"all" col/row - which is arithmetic mean - want IVW
betas[nrow(betas),] <- t(clock_res_ave_outcome$beta)
betas[,ncol(betas)] <- out_res_ave_clock$beta

colnames(betas)[ncol(betas)] <- "ALL"
rownames(betas)[nrow(betas)] <- "ALL"

ses<-with(df,proctabulate(se,outcome,clock,meanf))
ses[nrow(ses),] <- t(clock_res_ave_outcome$se)
ses[,ncol(ses)] <- out_res_ave_clock$se

colnames(ses)[ncol(ses)] <- "ALL"
rownames(ses)[nrow(ses)] <- "ALL"

ts<-betas/ses
ps<-pnorm(-abs(ts))*2


return( list(betas=betas,ses=ses,ps=ps))
}



annot_sig_n_hm <- function(df,fdr_thresh=0.1,nom_sig=.05,facets=NULL,root="stx",height=8.27,width=11.69){
  df<-annotate_sig(df,fdr_thresh,nom_sig)
  heatmap(df,root=root,facetw=facets,height=height,width=width)




}


annotate_sig <- function(df,fdr_thresh,nom_sig){

df$q<-p.adjust(df$p)
df$sig_level <- ifelse(df$q< fdr_thresh ,"*",ifelse(df$p< nom_sig , "+",""))

heading("tests performed")
print(with(df,proctabulate(n,meaning_out_s,clock,countf)))
heading("tests *")
print(with(df[df$sig_level=="*",],proctabulate(n,meaning_out_s,clock,countf)))
heading("tests +")
print(with(df[df$sig_level=="+",],proctabulate(n,meaning_out_s,clock,countf)))

heading("tests +|*")
print(with(df[df$sig_level %in% c("+","*"),],proctabulate(n,meaning_out_s,clock,countf)))


return(df)
}


heatmap<-function(df,root="st05_04_X",facetw=NULL,faceth=NULL,height=8.27,width=11.69,gt_size_sig=NULL){



dev.on(paste0(root,"heatmap"),graph_format,height=height,width=width)

text_file<-paste0(root,".txt")
write_pj1(df,text_file)

gg<-ggplot(data = df, aes(meaning_clock, meaning_out_s))+
   geom_tile(aes(fill=beta))+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0,space = "Lab", 
     name="beta")
 #  name=eval(parse(text="beta[age_acc]~acc.")) ) 

gg <- gg    +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() 


gg <- gg+ labs(x="Clock",y="Outcome")


if(!is.null(facetw)) gg<-gg+facet_wrap(as.formula(paste("~", facetw)))

if(!is.null(faceth)) gg<-gg+facet_grid(as.formula(paste(".~", faceth)))


if(!is.null(df$sig_level)) {

  if(is.null(gt_size_sig)) {
    gg<- gg + geom_text(data=df,aes(meaning_clock, meaning_out_s, label = sig_level))
    } else {
    gg<- gg + geom_text(data=df,aes(meaning_clock, meaning_out_s, label = sig_level),size=gt_size_sig)
    }
  }




print(gg)



heading(paste(root,"beta near one:-"))
print(df[order(abs(df$beta-1)),][1,c("clock","outcome","beta")])



dev.off()


}

heatmap_old<-function(df,root="st05_04_X",se_type_thresh=99){

#this crazy code puts two hms on top of each other


heading("Testing whether p is uniformly distributed")
print(ks.test(df$p,qunif))
dev.on(paste0(root,"qq"),graph_format)
qq(df$p,main="Association of age accelaration and outcome is enriched",sub="Distribution of p values for association of age acc. with outcome")
dev.off()

df$q<-p.adjust(df$p)
df$sig_level <- ifelse(df$q< fdr_thresh_heatmap ,"*",ifelse(df$p< nom_sig , "+",""))

#note just counting lines next not adding up n
heading("tests performed")
print(with(df,proctabulate(n,meaning_out_s,clock,countf)))
heading("tests *")
print(with(df[df$sig_level=="*",],proctabulate(n,meaning_out_s,clock,countf)))
heading("tests +")
print(with(df[df$sig_level=="+",],proctabulate(n,meaning_out_s,clock,countf)))

dev.on(paste0(root,"heatmap"),graph_format)


gg<-ggplot(data = df[df$se <= se_type_thresh,], aes(clock, meaning_out_s))+
	 geom_tile(aes(fill=beta))+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0,space = "Lab", 
   name=eval(parse(text="beta[age_acc]~acc.")) ) 


 if (se_type_thresh<99) {gg<-gg+
  new_scale("fill") +
  	 geom_tile(data=df[df$se > se_type_thresh,],aes(fill=beta))+
 scale_fill_gradient2(low = "lightsteelblue", high = "rosybrown4", mid = "white", 
   midpoint = 0 ,space = "Lab", 
   name=eval(parse(text="beta[age_acc]~acc."))  )}

gg <- gg    +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() 



gg<- gg + geom_text(data=df,aes(clock, meaning_out_s, label = sig_level)) +
 labs(x="Clock",y="outcome Group")



print(gg)

dev.off()


}



rescale_clock<-function(res,standardise_scores,one_age_eff,res_age=NULL,beta_age_all=NULL){

#here we are going to rescale the 
if(standardise_scores) {
  factor_mult_beta<-res$score_sd
  res<- rescale_ass(res,factor_mult_beta, changing="score")
} else {

  if(!one_age_eff) {
    key_vars <- c("outcome" , "beta")
    res_age <- res_age[,key_vars]
    res_age$factor_mult_beta<-1/res_age$beta
    res_age <- remove_df_cols(res_age,"beta")
    res<-merge(res,res_age)
#res now has a column factor_mult_beta= 1/beta$age
  }
  if(one_age_eff) res$factor_mult_beta<-1/beta_age_all
  

  res<- rescale_ass(res,res$factor_mult_beta, changing="trait")

#res <- remove_df_cols(res,"factor_mult_beta")
}


print(head(res))
return(res)
}


rescale_ass<- function(assoc_df,factor_mult_beta, changing){
#think the changing had sd 3 and we want sd 1

  if (changing=="score"){
    assoc_df$score_mean <- assoc_df$score_mean / factor_mult_beta
    assoc_df$score_sd <- assoc_df$score_sd / factor_mult_beta

  }


  if (changing=="trait"){
    assoc_df$trait_mean <- assoc_df$trait_mean * factor_mult_beta
    assoc_df$trait_sd <- assoc_df$trait_sd * factor_mult_beta

  }


    assoc_df <- ratio_beta_n_se(assoc_df,1/factor_mult_beta)



  return(assoc_df)


}





ratio_beta_n_se<-function(df,value){
  df$beta<-df$beta/value
  df$se<-df$se/abs(value)
  return(df)
}



ratio_n_se <- function(df){
beta<-df$beta.x/df$beta.y
se<-sqrt( (df$se.x/df$beta.y)^2 + (df$se.y/df$beta.x)^2)
return(data.frame(beta=beta,se=se))
}



standardise_assoc <- function(assoc_df,vars=c("trait","score")){

  if ("trait" %in% vars){
    assoc_df <- ratio_beta_n_se(assoc_df,assoc_df$trait_sd)
    
    assoc_df$trait_mean <- assoc_df$trait_mean / assoc_df$trait_sd
    assoc_df$trait_sd <- 1

  }

  if ("score" %in% vars){
     assoc_df <- ratio_beta_n_se(assoc_df, 1 / assoc_df$score_sd)
    assoc_df$score_mean <- assoc_df$score_mean / assoc_df$score_sd
    assoc_df$score_sd <- 1

  }

  return(assoc_df)


}


check_confounder <- function(res,res_sm,outcomes,plot_file){

 res_bth<-merge(res,res_sm,by=c("clock","outcome"))


res_bth <- res_bth[res_bth$outcome %in% outcomes,]

res_bth <- cbind(res_bth,ratio_n_se(res_bth))
#beta is now the ratio

heading("ratio of effect to effect when confounder fitted i.e. beta/beta_with_sm \n H0 : beta =1 ")
ratios<-show_beta_n_se_tables(res_bth)

ratios_minus1<-ratios;ratios_minus1$betas<-ratios_minus1$betas-1
ts <- ratios_minus1$betas/ratios_minus1$ses
ratios_minus1$ps<-pnorm(-abs(ts))*2
rm(ts)
ratios$ps <- ratios_minus1$ps

print(ratios)


heading("distribution of p values for beta beta/beta_with_confounder \n H0 : beta =1 p rounded down")
print(table(floor(ratios$p/.05)*.05))

dev.on(plot_file,graph_format)

ggp<-ggplot(data = res_bth[res_bth$se.x<.5& res_bth$se.y<.5,],aes(x = beta.x,y = beta.y,colour=clock)) + 
    geom_point(alpha=.2) + 
    geom_errorbar(aes(ymin = beta.y- 2*se.y,ymax = beta.y+2*se.y),alpha=0.05) + 
    geom_errorbar(aes(xmin = beta.x- 2*se.x,xmax = beta.x+2*se.x),alpha=0.05) +
    geom_abline(slope=1,intercept=0,alpha=0.2)+
    coord_fixed(ratio=1)+
    labs(
  title = "Age accelaration",
  subtitle = "not mediated through smoker status")+
xlab("Beta_age_acc")+
ylab("Beta_age_acc with smoker covariate")
print(ggp)
dev.off()
}





show_age_effects <- function(res_age,file_out){

detected_disease <- any(grepl("-",res_age$outcome))

if(detected_disease){
res_age$category<-left(res_age$meaning_out,1)
all_eff<-res_age$beta[res_age$outcome=="ALL"]
} else {
  res_age$category<-"Q"
all_eff<- ivm_t_test(res_age,"category")$beta[1]
}

dev.on(file_out,graph_format)
gg_obj <- ggplot(res_age, aes(x=meaning_out_s,y=beta,fill=category)) +
geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=beta-2*se, ymax=beta+2*se), width=.2,
                 position=position_dodge(.9)) +
   geom_hline(aes(yintercept=all_eff),linetype="dashed" ) +
  coord_flip() + xlab("Outcome") + ylab("Effect of chronologocal age (years) on outcome")
  print(gg_obj)
  dev.off()

res_age <- res_age[order(res_age$beta),]
heading("Least Association with age")
print(head(res_age,n=5))
heading("Most Association with age")
print(tail(res_age,n=5))

display_matching_cond(res_age,res_age$beta<0,"outcomes where beta_chron_age<0")

p_sugg<-.05
z_th<-abs(qnorm(p_sugg/2))
title_disp<-  paste0("outcomes where beta_chron_age\n suggestively(p<",p_sugg,") differs from beta_chron_age_all=",all_eff)
display_matching_cond(res_age,abs((res_age$beta-all_eff)/res_age$se) > z_th,title_disp)


return(all_eff)

}



move_col_x_to_after_y <- function(df,x,y){

n_to_move <- which(names(df)==x)


df0<-df[,-n_to_move]
new_loc<-which(names(df)==y)+1

df2<-df[,n_to_move,drop=F]
df1<-df0[,1:new_loc-1]
df3<-df0[,new_loc:ncol(df0)]

df_ans <- cbind(df1,df2,df3)

return(df_ans)

}

melt_mat <- function(mat,names_vars){
df<-as.data.frame(mat)
df$id<-row.names(df)
df1<-melt(df,id_vars="id")
names(df1)<-names_vars
return(df1)
}


add_meaning_out <- function(df,definitions,length_short=25){
definitions<-rename_df_cols(definitions,"meaning","meaning_out")
df1<-merge(definitions,df,by="outcome", all.y=T)
df1$meaning_out <- ifelse(is.na(df1$meaning_out),df1$outcome,df1$meaning_out)
df1$meaning_out_s <- gsub("Malignant neoplasm","MN",df1$meaning_out)
df1$meaning_out_s <- left(df1$meaning_out_s,length_short)
return(df1)
}


add_meaning_clock <- function(df,definitions){
definitions<-rename_df_cols(definitions,c("meaning","meaning_s"),c("meaning_clock","meaning_clock_s"))
df1<-merge(definitions,df,by="clock", all.y=T)
return(df1)
}

xform_outcome_map <- function(file_in,res) {

if (grepl("RDS",file_in)){
icd_def<-readRDS(file_in)
out_def<-icd_def$blocks
out_def<-rename_df_cols(out_def,"icd_block","outcome")
out_def<-out_def[out_def$icd_letter %in% unique(res$outcome),]
out_def<- out_def[,c(1,2)]


names(out_def)<- c("Block","Title")
out_def$Title <- substr(out_def$Title,8,nchar(out_def$Title))
write_pj1(out_def,file="st05_00_out_defs.tsv")
names(out_def)<- c("outcome","meaning_out")
out_def$meaning_out <- paste(out_def$outcome,out_def$meaning_out)
print(out_def)

} else {

  out_def <- read_pj(file_in)
}


return(out_def)
}








covers_95_ci_not <-function(df,value){



heading(paste("Estimated effects 95 CI that do not not cover", value))



df$lc<-df$beta-2*df$se
df$uc<-df$beta+2*df$se

out<-df$lc>value | df$uc<value
return(df[out,])
}



qc<-function(df,crit,desc){
  heading(paste("Qcing out",desc))
  print(table(crit))
  print(df[!crit,])
  df<-df[crit,]
  print(dim(df))
  return(df)

}



display_beta<-function(res_all,compact_out=FALSE){

if( ! compact_out){
heading("Count betas observed")
print(with(res_all,proctabulate(beta,outcome,clock,countf)))


heading("Number of times age residual associates postively with outcome")
heading("Count b>0")
print(with(res_all,proctabulate(beta>0,outcome,clock,sumf)))
heading("Mean (beta>0) i.e. %freq with which beta>0")
print(with(res_all,proctabulate(beta>0,outcome,clock,meanf)))
}


heading("beta for ungrouped outcomes")
just_blocks<-res_all[nchar(res_all$outcome)>1 & res_all$outcome != "ALL",]
res_blocks_xtabbed<-show_beta_n_se_tables (just_blocks)
print(res_blocks_xtabbed)
overall_ave_clock_beta <- res_blocks_xtabbed$betas["ALL","ALL"]





bar_plot_pj <- function(df){

p<-ggplot(data=df,aes(x=meaning,y=beta))+geom_bar(stat='identity',fill="steelblue") +
  geom_text(aes(label=formatC(beta, digits = 2, format = "f")), vjust=-0.3, size=3.5)+
    geom_errorbar( aes(x=meaning, ymin=beta-2*se, ymax=beta+2*se),
     width=0.4,colour="grey")    +
  theme(axis.text.x = element_text(angle = 90))+labs(x="outcome",y="beta (~95% CI) ")
print(p)
}

heading("t1")
y<-res_blocks_xtabbed$betas["ALL",]
se<-res_blocks_xtabbed$ses["ALL",]
bp_data<-data.frame(clock=names(y),beta=y,se=se)
bp_data<-merge(bp_data,rbind(clock_map,data.frame(clock="ALL",meaning="ALL")))
bar_plot_pj(bp_data)
heading("t2")
y<-res_blocks_xtabbed$betas[,"ALL"]
se<-res_blocks_xtabbed$ses[,"ALL"]
bp_data<-data.frame(outcome=names(y),beta=y,se=se)
bp_data<-merge(bp_data,rbind(outcome_map,data.frame(outcome="ALL",meaning="ALL")))
bp_data$meaning<-gsub("Malignant neoplasm","MN",bp_data$meaning)
bp_data$meaning<-left(bp_data$meaning,20)
bar_plot_pj(bp_data)
heading("t3")



if( ! compact_out){
write_pj1(add_rownames_as_first_col(res_blocks_xtabbed$betas),"st05_06a_beta_xtab.tsv")
write_pj1(add_rownames_as_first_col(res_blocks_xtabbed$ses),"st05_06a_se_xtab.tsv")

heading("Mean beta for grouped groups")
print(show_beta_n_se_tables (res_all[nchar(res_all$outcome)==1,]))

heading("And ALL")
print(show_beta_n_se_tables (res_all[res_all$outcome=="ALL",]))


just_betas<- melt_mat(res_blocks_xtabbed$betas,names=c("outcome","clock","beta"))


excluding_alls<-just_betas[just_betas$clock != "ALL" & just_betas$outcome != "ALL",]
outcome_ivw_aves<-just_betas[just_betas$clock != "ALL" & just_betas$outcome == "ALL",]

p<-ggplot(excluding_alls, aes(x=beta)) +
  geom_density() + facet_wrap(.~clock,nrow=2) + 
  xlim(c(-0.5,2)) +  
    geom_vline(aes(xintercept=0)) +
        geom_vline(aes(xintercept=1),colour="blue") +

  geom_vline(data=outcome_ivw_aves, aes(xintercept=beta),linetype="dashed") 


  print(p)

}
return(overall_ave_clock_beta)
}


