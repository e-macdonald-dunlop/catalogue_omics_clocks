
source("/exports/igmm/eddie/wilson-lab/apps_by_us/gen_useful.R")
source("../scripts/helper_fns.R")
source("../scripts/helper_fns_plot.R")

#=================================

options(width=400)
options(digits=4)
library(reshape2)
library(ggplot2)

library(qqman)
library(plyr)
library(ggnewscale)
library(stringr)

#https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/

date()
heading("start")

args<-commandArgs(T)
quant_trait <- as.logical(args[2])
standardise_scores<-as.logical(args[3])
standardise_scores<-ifelse(is.na(standardise_scores),FALSE,standardise_scores)


unhelpful_colnames<-c("formula", "factor_mult_beta","meaning_out","method" ,"submethod" ,"features" ,"beta_type")



heading(paste("Analysing assoc QT is ",quant_trait))
heading(paste("Analysing assoc stdise scores is ",standardise_scores))


normal_analysis_type=data.frame(method="fixed_alpha",submethod="all",features="full")


fdr_thresh_heatmap <- 0.1
nom_sig <- 0.05
one_sided <- T
se_max_for_hm <- 0.5
chapters_of_int <- c("C","E","I","J")
meaning_out_s_chars <- 35
graph_format <- "pdf"


std_analysis_type <- data.frame(method="fixed_alpha",submethod="all",features="full",beta_type="full")
std_shrunk_analysis_type <- std_analysis_type
std_shrunk_analysis_type$beta_type <- "shrunken"
all_shrunk_analysis_type <- remove_df_cols(std_shrunk_analysis_type,"submethod")

#bayes prior assumptions for distn of beta
#mixture distn of point mu1 and rnorm mu2,

n_prior<-100000
prop_block<-0.0001
mu1<-0
mu2<-0
sd2<-ifelse(standardise_scores,.25,1)
#if scores are not standardised they are on age scale -> beta ~ 0-1
#if  standardised they are on scale free  -> beta ~ 0-0,05


file_clock_map<-"../data/clock_names.tsv"

if (quant_trait) {

input_res_pref<-"../data/clock_pheno_assoc/final/"
res_values <- c("age","no_age","smoking/age")
file_outcome_map<- "../data/qt_mapper.txt"

traits_to_exc<-c("cigsperday","edu","hdl","age_at_vene","ascvd")
traits_to_exc_hm<-c("totchol")
traits_to_reverse <- c("cortisol_nmol_l","fev1")

fdr_thresh_age <- NA
emp_case_thresh_age <- NA


} else {
input_res_pref<-"../p04_assoc_"
res_values <- c("age","no_age","age_smk")
file_outcome_map<- "..//p03a_process_icd_defs/st03a_01_out_defs.tsv"


traits_to_exc_hm<-NULL
traits_to_exc<-NULL
traits_to_reverse <- NULL

fdr_thresh_age <- 0.1
emp_case_thresh_age <- 5


}


gt_size_sig<-NULL
if (quant_trait) gt_size_sig <-10

clock_names<-   read_pj(file_clock_map)

sub_folder_suff<-"/p82_comb_assoc/st82_01_assoc_res.txt"
res_files<-paste0(input_res_pref,res_values,sub_folder_suff)

assc_res_file <- res_files[1]
assc_age_file <- res_files[2]
assc_smk_file <- res_files[3]

res<-read_in_assoc(assc_res_file,one_sided=T,quant_trait,clock_names,traits_to_reverse,traits_to_exc,age_only=F)
res_sm<-read_in_assoc(assc_smk_file,one_sided=T,quant_trait,clock_names,traits_to_reverse,traits_to_exc,age_only=F)
res_age <-read_in_assoc(assc_age_file,one_sided=T,quant_trait,clock_names,traits_to_reverse,traits_to_exc,age_only=T)




outcome_map<-   read_pj(file_outcome_map)
outcome_map<-outcome_map[outcome_map$outcome %in% unique(res_age$outcome),]


if(!quant_trait) {
out_def<-outcome_map
names(out_def)<- c("Block","Title")
out_def$Title <- substr(out_def$Title,8,nchar(out_def$Title))
write_pj1(out_def,file="st05_00_out_defs.tsv")
}

clock_map <- clock_names[,-1]
clock_map<-clock_map[clock_map$clock %in% unique(res$clock),]
clock_map$meaning<-gsub("_"," ",clock_map$meaning)



heading("Association of Chron Age with outcome")

res_age$q <- p.adjust(res_age$p)
res_age <- add_meaning_out(res_age,outcome_map,length_short=meaning_out_s_chars)

res_age_to_out <- remove_df_cols(res_age,"formula")
write_pj1(res_age_to_out,"st05_01_res_age.tsv")


res_age_qcd <- res_age


if(!is.na(fdr_thresh_age))res_age_qcd<-qc(res_age_qcd,res_age_qcd$q<fdr_thresh_age,paste("res_age$q<",fdr_thresh_age,"\n exclusions"))
if(!is.na(emp_case_thresh_age)) res_age_qcd<-qc(res_age_qcd,res_age_qcd$cases>=emp_case_thresh_age,paste("Subsetting res_age to those n cases >",emp_case_thresh_age,"\n exclusions"))


#we are reversing here nut no need for main analysis as we divide by beta so b_cc/B_age reverses itself 
heading("Unqcd chron age results")
show_age_effects(res_age,"st05_02_age_eff_all")




heading("Qcd chron age results - i.e. only look at traits where age has an effect")
beta_age_all <- show_age_effects(res_age_qcd,"st05_03_age_eff_qcd_only")

heading("Does it look like the effect of age is common to all traits")
covers_95_ci_not(res_age_qcd,beta_age_all)



res_std_only<-pick_some_methods(res,normal_analysis_type)
res_smk_std_only <- pick_some_methods( res_sm,normal_analysis_type)

heading("Are results a smoking effect?")
if(!is.null(res_sm)) {
  check_confounder(res_std_only,res_smk_std_only,res_age_qcd$outcome,"st05_04_betaxtab_sm_nonsm")
} else {
  heading("no smoking here ;)")
}



#heading("Pre qc standard results goodish q value")
#add_meaning_out(passed_fdr(res_std_only,0.15),outcome_map)















#following step would limit qt to those where res_age applies, but not for stdised scores
#so make explicit for those


heading("Now taking Beta_age_acc, for outcomes where beta_age passed qc")
#note not dropping awkward traits from stdise but are otherwise
if (standardise_scores) {
  res_all <- rescale_clock(res,standardise_scores=T)
  res_all <- res_all[res_all$outcome %in% res_age_qcd$outcome,]
  } else  {
  if (quant_trait)  res_all <- rescale_clock(res,standardise_scores=F,one_age_eff=FALSE,res_age=res_age_qcd) else  {
    res_all <- res[res$outcome %in% res_age_qcd$outcome,]
    res_all  <- rescale_clock(res_all,standardise_scores=F,one_age_eff=TRUE, beta_age_all=beta_age_all)
  }
}


#if(quant_trait) one_age <- F else one_age<-T
#res_all<-rescale_clock(res_all,standardise_scores,one_age,res_age_qcd,beta_age_all)



heading("Bayes - prior expected beta = 0 - standard clock method")

#prior<-block_n_norm(n_prior,prop_block,mu1,mu2,sd2)
#res<-bayes_est_beta(prior,res_all$beta,res_all$se,summary=TRUE,plots=F)

res<-bayes_est_beta_simple_prior(mu2,sd2,res_all$beta,res_all$se)


res_all_b_0<-cbind(res_all[,!(names(res_all) %in% c("beta",     "se"))],res)
res_all_b_0<-res_all_b_0[,names(res_all)]
res_all$beta_type <- "full"
res_all_b_0$beta_type <- "shrunken"
#note the p-values here are still frequentist

res_all <- rbind(res_all,res_all_b_0)
saveRDS(res_all,"st05_05_res_all.RDS")











#---------------- restart



res_all<-readRDS("st05_05_res_all.RDS")

 res_all<- add_meaning_out(res_all,outcome_map)
  res_all<- add_meaning_clock(res_all,clock_map)




res_all <- move_col_x_to_after_y(res_all,"beta_type","features")

res_std_only <- pick_some_methods(res_all,std_analysis_type)



heading("tests run")
with(res_std_only,proctabulate(n,clock,outcome,meanf))


heading("cases run")
with(res_std_only,proctabulate(cases,clock,outcome,meanf))



heading("Standarddisation")


heading("trait_sd")
with(res_std_only,proctabulate(trait_sd,clock,outcome,meanf))

heading("score_sd")
with(res_std_only,proctabulate(score_sd,clock,outcome,meanf))



dev.on("st05_06a_betas",graph_format)
beta_all<- display_beta(res_std_only)

heading(paste("overall ave beta",beta_all))
dev.off()

heading("Sign test of beta by outcome")
a<-sign_and_show(res_std_only,"outcome")
a$q<-p.adjust(a$p)
a<-a[order(a$p),]
print(a)



heading("Sign test of beta by clock")
a<-sign_and_show(res_std_only,"clock")
a$q<-p.adjust(a$p)
a<-a[order(a$p),]
print(a)






res_std_only <- annotate_sig(res_std_only,fdr_thresh_heatmap,nom_sig)


res_to_pub <- res_std_only 
res_to_pub <- remove_df_cols (res_to_pub,unhelpful_colnames)
res_to_pub<-res_to_pub[order(res_to_pub$p),]
heading("Base analysis - passed FDR")
print(res_to_pub[res_to_pub$sig_level=="*",])
heading("Saving ase analysis")
write_pj1(res_to_pub,"st05_06_res_std_clock_full_beta.tsv")


#chosen_analysis_type_temp <-std_analysis_type 
#heading(paste("counts of nominal significance of beta - one_sided=",one_sided))
#heading("dimensions inc final row.col tots")
#for(subm in unique(res_all$submethod)){
#	chosen_analysis_type_temp$submethod<-subm
#heading(subm)
#res_to_analyse <- pick_some_methods( res_all,chosen_analysis_type_temp)
#res_to_analyse<-res_to_analyse[res_to_analyse$p<nom_sig,]
#print(dim((with(res_to_analyse, proctabulate(outcome,outcome,clock,countf)))))
#print(with(res_to_analyse, proctabulate(outcome,outcome,clock,countf)))
#}



heading("heatmapping std res - all ses")
heatmap(res_std_only,root="st05_07_01_res_all")
# heatmap(res_std_only[!(res_std_only$trait %in% traits_to_exc_hm),],root="st05_07_01a_res_all_some_excs")





#heading("heatmapping std res - all ses")
#annot_sig_n_hm (res_std_only,fdr_thresh=fdr_thresh_heatmap,nom_sig=nom_sig,root="st05_07_01_res_all")
#annot_sig_n_hm (res_std_only[!(res_std_only$trait %in% traits_to_exc_hm),],fdr_thresh=fdr_thresh_heatmap,nom_sig=nom_sig,root="st05_07_01a_res_all_some_excs")


heading("powered res")

res_powered<- res_std_only[res_std_only$se< se_max_for_hm,]


heatmap(res_powered,root="st05_07_02_res_all_powered",gt_size_sig=gt_size_sig)




#weird <- abs(res_std_only$beta)>1
#print(res_std_only[weird,])
#res_std_only <- res_std_only[!weird,] 

#heatmapping very large SE data was giving nonsense -> impose a moderate prior to shrink noisy ests
#alternative to favour a block at zero not used.



#but noisy data was pulled to zero. Perhaps our best estimate of prior is what we observed from powered?

heading("Bayes - prior expected beta = 0 - standard clock method")
res_std_shrunk <- pick_some_methods(res_all,std_shrunk_analysis_type)
#this sgould be identical to unshrunk analysis - P is unchanged
res_std_shrunk <- annotate_sig(res_std_shrunk ,fdr_thresh_heatmap,nom_sig)


heatmap(res_std_shrunk,root="st05_08_01_bayes_prior_mean_0_all",gt_size_sig=gt_size_sig)

heading("Bayes - prior expected beta = 0 - standard clock method some traits only")
res_bayes_excs <- res_std_shrunk [!(res_std_shrunk$trait %in% traits_to_exc_hm),]
heatmap(res_bayes_excs,root="st05_08_02_bayes_prior_mean_0_for_pub_some_traits_exc",gt_size_sig=gt_size_sig)

#look atr explanatiry power
if(standardise_scores){



heading("Bayes - prior expected beta = 0 - faceted clock method some traits only  BMI further exc")

faceted_data <- pick_some_methods(res_all,all_shrunk_analysis_type)
faceted_data <-faceted_data[!(faceted_data$trait %in% c(traits_to_exc_hm,"bmi")),]


fd_pcs <- faceted_data[grepl("pcs",faceted_data$submethod) | grepl("all",faceted_data$submethod),]

#drop the clonomic clock - its not got lots of pcs and the storng signal distorts the colouring
fd_pcs <- fd_pcs[fd_pcs$clock != "pheno",]

fd_pcs$sm_pc_ord <- ordered(fd_pcs$submethod,levels=c("3_pcs",  "5_pcs",  "10_pcs", "20_pcs",  "all"))

if(quant_trait) {
	heatmap(fd_pcs ,facetw="sm_pc_ord",
root="st05_08_03_bayes_prior_mean_0_for_pub_by_pcmethod_ord",height=10,width=20) } else
	{
	heatmap(fd_pcs ,faceth="sm_pc_ord",
root="st05_08_03_bayes_prior_mean_0_for_pub_by_pcmethod_ord",height=10,width=20)}





all_unshrunk_analysis_type <- remove_df_cols(std_shrunk_analysis_type,"submethod")
all_unshrunk_data <- pick_some_methods(res_all,all_unshrunk_analysis_type)
pc_all_data<-all_unshrunk_data
pc_all_data<-pc_all_data[grepl("pc",pc_all_data$submethod)|(pc_all_data$submethod=="all"),]
pc_all_data <- pc_all_data[pc_all_data$clock != "pheno",]
pc_all_data <-pc_all_data[!(pc_all_data$trait %in% c(traits_to_exc_hm,"bmi")),]


for (pc_type in c("3_pcs","all")){

  heading(pc_type)
  one_pc_type_data<-pc_all_data[pc_all_data$submethod==pc_type,]
  display_beta(one_pc_type_data,compact_out=T)
}


}



#fd_aminb <- faceted_data[!grepl("pcs",faceted_data$submethod),]
#fd_aminb$sm_aminb_ord <- ordered(fd_aminb$submethod,levels=c("a_only",    "b_only" ,   "a_minus_b",    "all"))

   
#heatmap(fd_aminb ,faceth="sm_aminb_ord",root="st05_08_04_bayes_prior_mean_0_for_pub_by_minusmethod",
#  height=10,width=20)




#annot_sig_n_hm (faceted_data,facets="submethod",
#	fdr_thresh=fdr_thresh_heatmap, nom_sig=nom_sig,root="st05_08_03_bayes_prior_mean_0_for_pub_by_method")

  

date()
heading("done")






q()











