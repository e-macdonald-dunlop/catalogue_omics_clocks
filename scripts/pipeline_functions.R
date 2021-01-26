
heading<-function(text){
	writeLines("")
	writeLines("==================================================")
	writeLines(text)
	writeLines("==================================================")
	writeLines("")
}

get_data_from_yml_file<-function(path){
  yml<-yaml.load_file(path)
  yml<-lapply(yml,function(x){ifelse(x=="NA",NA,x)})
  return(yml)
}

get_paramaters_from_config<-function(){
	if(!file.exists("../control_files/config.yml"))stop("Error: No config file found in expected directory, ckeck location!")
	config<-get_data_from_yml_file("../control_files/config.yml")
	if(!is.na(config$covariate_file)){
		covariate_file<-unlist(strsplit(config$covariate_file," "))
	}else{
		covariate_file<-character(0)
	}
	omics_file<-unlist(strsplit(config$omics_file," "))
	if(!is.na(config$validation_sample)){
		validation_covariate_file<-unlist(strsplit(config$validation_covariate_file," "))
	}else{
		validation_covariate_file<-NA
	}
	steps<-config$steps
	steps<-unlist(strsplit(steps," - "))
	analysis_file<-get_data_from_yml_file(config$analysis_file)
	if(is.na(analysis_file$standard_covariates)){standard_covariates<-character(0)}else{standard_covariates<-analysis_file$standard_covariates}
	outcome<-unlist(strsplit(analysis_file$formula," ~ "))[1]
	additional_covariates<-unlist(strsplit(analysis_file$formula," ~ "))[2]
	additional_covariates<-additional_covariates[additional_covariates!="standard_covariates"]
	additional_covariates<-ifelse(is.na(additional_covariates),character(0),additional_covariates)
	if(!is.na(analysis_file$cols_to_drop)){
		to_drop<-unlist(strsplit(analysis_file$cols_to_drop," - "))
	}else{
		to_drop<-character(0)
	}
	if(!is.na(analysis_file$pre_correction_covariates)){
		pre_correction_covariates<-unlist(strsplit(analysis_file$pre_correction_covariates," \\+ "))
	}else{
		pre_correction_covariates<-character(0)
	}
	return(params<-list(omics_file=omics_file,
		cohort=config$cohort,
		omic_assay=config$omic_assay,
		covariate_file=covariate_file,
		steps=steps,
		outcome=outcome,
		standard_covariates=standard_covariates,
		pre_correction_covariates=pre_correction_covariates,
		additional_covariates=additional_covariates,
		to_drop=to_drop,
		covariate_list=c(outcome,standard_covariates,additional_covariates,pre_correction_covariates),
		bespoke_omics_rscript=config$bespoke_omics_rscript,
		bespoke_pheno_rscript=config$bespoke_pheno_rscript,
		bespoke_omics_qc_rscript=config$bespoke_omics_qc_rscript,
		zscore_raw=config$zscore_raw,
		validation_omics_file=ifelse(!is.na(config$validation_omics_file),config$validation_omics_file,NA),
		validation_covariate_file=validation_covariate_file,
		validation_sample=config$validation_sample,
		zscore_resid=config$zscore_resid,
		iteration=config$iteration,
		method=config$method,
		if_elastic_net=config$if_elastic_net))
}

make_names_r_readable<-function(data){
  names(data)<-gsub("/",".",names(data))
  names(data)<-gsub("%",".",names(data))
  names(data)<-gsub("\\(",".",names(data))
  names(data)<-gsub("\\)",".",names(data))
  names(data)<-gsub("\\+",".",names(data))
  names(data)<-gsub(" ","",names(data))
  names(data)<-gsub(":",".",names(data))
  names(data)<-gsub("-",".",names(data))
  names(data)<-gsub("\\[",".",names(data))
  names(data)<-gsub("\\]",".",names(data))
  names(data)<-tolower(names(data))
  return(data)
}


get_omics_data<-function(path,cols_to_drop){
	if(!file.exists(path[1]))stop("Error: No file Omics file found!")
	data<-fread(paste0(path[1]),header=T,drop=cols_to_drop,stringsAsFactors=F,data.table=F)
	data<-make_names_r_readable(data)
	names(data)[1]<-"iid"
	if(length(path)>1){
		for(file in path[-1]){
			if(!file.exists(file))stop("Error: No file Omics file found!")
			df<-fread(paste0(file),header=T,drop=cols_to_drop,stringsAsFactors=F,data.table=F)
			df<-make_names_r_readable(df)
			names(df)[1]<-"iid"
			data<-merge(data,df,by="iid",all=TRUE)
		}
	}
	return(data)
}

get_covariate_data<-function(path){
	if(!file.exists(path[1]))stop("Error: No file covariate file found!")
	data<-fread(paste0(path[1]),data.table=FALSE)
	data<-make_names_r_readable(data)
	if(!("iid"%in%colnames(data))){data$iid<-data$idcode}
	#names(data)[1]<-"iid"
	if(length(path)>1){
		for(file in path[-1]){
			if(!file.exists(file))stop("Error: No file covariate file found!")
			df<-fread(paste0(file),data.table=FALSE)
			df<-make_names_r_readable(df)
			if(!("iid"%in%colnames(df))){df$iid<-df$id}
			if(!("iid"%in%colnames(df))){df$iid<-df$idcode}
			data<-merge(data,df,by="iid",all=TRUE)
		}
	}
	return(data)
}

get_p00_data<-function(){
	data<-fread("./p00_pipeline_phenotypes/st00_all_pipeline_data.tsv",header=T,stringsAsFactors=F,data.table=F)
	return(data)
}

get_distribution_plot<-function(vec,threshold=NULL){
  if(!missing(threshold)){
    vec[,"status"]<-"Raw" #ifelse(abs(scale(vec[,1],scale=TRUE,center=TRUE))>threshold,"outlier","fine")
    filtered<-vec[abs(scale(vec[,1],scale=TRUE,center=TRUE))<threshold,]
    filtered[,"status"]<-"Outliers Removed"
    if(nrow(vec)>nrow(filtered)){
    	vec<-rbind(vec,filtered)
    	#vec[,"status"]<-factor(vec[,"status"],levels=c("Raw","Outliers Removed"))
   		p <- ggplot(vec, aes(vec[,1],fill=status)) + 
      		geom_density(alpha = 0.5,colour="grey") +
	      	facet_grid(. ~ factor(status,levels=c("Raw","Outliers Removed")),scales="free") +
    	  	theme_light(base_size = 14) +
      		xlab("") +
      		theme(legend.position="bottom")
    }else{
    	p <- ggplot(vec, aes(vec[,1])) + 
      		geom_density(alpha = 0.5,fill="#1c9099",colour="grey") +
      		theme_light(base_size = 14) +
      		xlab("")
    }
  }else{
    p <- ggplot(vec, aes(vec[,1])) + 
      geom_density(alpha = 0.5,fill="#1c9099",colour="grey") +
      theme_light(base_size = 14) +
      xlab("")
  }
    return(p)
}

zscore_filter<-function(vec,threshold){
	filtered<-ifelse(abs(scale(vec,scale=TRUE,center=TRUE))>threshold,NA,vec)
	return(filtered)
}

apply_zscore_filter<-function(df,threshold){
	if(!is.na(threshold)){
		filtered<-apply(df,2,zscore_filter, threshold=threshold)
	}else if(is.na(threshold)){
		filtered<-df
	}
	return(filtered)
}

plot_NA_heatmap<-function(data){
	long_data<-reshape2::melt(data)
	long_data$value<-ifelse(!is.na(long_data$value),1,long_data$value)
	for_plot<-long_data[!complete.cases(long_data),]
	x<-ggplot(for_plot, aes(variable, iid)) +
	  geom_tile(aes(fill = value), colour = "white") +
	  scale_fill_continuous(na.value = 'red')
	return(x)
}

do_bespoke_qc<-function(df,bespoke_script){
		if(!is.na(bespoke_script)){
			heading(paste0("Running ",bespoke_script,"..."))
			source(bespoke_script)
			df<-bespoke(df)	
		}	
	return(df)
}


remove_NAs_validation<-function(validation_df,clean_clock_df){
	validation_df<-validation_df[,colnames(clean_clock_df)]
	validation_df[validation_df == 0] <- NA
	validation_df<-trim_rows(validation_df)
	return(validation_df)
}


prop_nas<-function(vec){
  sum(is.na(vec))/length(vec)
}


trim_cols<-function(df,threshold){
  drops<-apply(df,2,prop_nas)>threshold
  print(names(df)[drops])
  df<-df[,!drops]
  return(df)
}

trim_rows<-function(df,threshold=0){
  drops<-apply(df,1,prop_nas)>threshold
  print(df[drops,])
  df<-df[!drops,]
  return(df)
}


remove_missing_nmr<-function(data_orc,col_threshold,row_threshold=0){
	data_trim<-trim_cols(data_orc,col_threshold)
	data_trim_again<-trim_rows(data_trim)
	df<-data_trim_again
	return(df)
}



missingness_by_predictor<-function(df){
  x<-data.frame(missingness=apply(df,2,function(x) sum(is.na(x))))
  ggplot(x,aes(x=missingness)) +
    geom_density(fill="#1c9099",alpha = 0.5,colour="grey")
  ggplot(x, aes(missingness)) +
    geom_histogram()
  ggsave("missingness_by_predictor.png")
}

missingness_by_sample<-function(df){
  x<-data.frame(missingness=apply(df,1,function(x) sum(is.na(x))))
  ggplot(x,aes(x=missingness)) +
    geom_density(fill="#1c9099",alpha = 0.5,colour="grey")
  ggplot(x, aes(missingness)) +
    geom_histogram()
  ggsave("missingness_by_sample.png")
}


do_bespoke<-function(df,bespoke_script,validation=FALSE){
	if(validation==TRUE){
		if(!is.na(bespoke_script)){
			heading(paste0("Running ",bespoke_script,"..."))
			source(bespoke_script)
			df<-bespoke_validation(df)
		}	
	}	
	if(validation==FALSE){
		if(!is.na(bespoke_script)){
			heading(paste0("Running ",bespoke_script,"..."))
			source(bespoke_script)
			df<-bespoke(df)	
		}	
	}
	return(df)
}

get_distribution_stats<-function(df){
	skewness<-apply(df,2,skewness,na.rm=T)
	kurtosis<-apply(df,2,kurtosis,na.rm=T)
	return(list(skewness=skewness,kurtosis=kurtosis))
}

make_distribution_plots<-function(df,threshold,sample=""){
	rownames(df)<-df[,1]
	data<-df[,-1]
	data<-as.matrix(data)
	####
	#want to make a dataframe with the skewness + kutosis stats for each predictor
	dist_stats<-data.frame(matrix(NA,nrow=ncol(data),ncol=4))
	names(dist_stats)<-c("predictor","skewness","kurtosis","binomial")
	dist_stats[,1]<-colnames(data)
	output<-get_distribution_stats(data)
	dist_stats[,2]<-output$skewness
	dist_stats[,3]<-output$kurtosis
	print(dist_stats)
	write.table(dist_stats,paste0("distribution_statistics_",sample,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
	skew<-dist_stats[order(abs(dist_stats[,"skewness"])),]
	best_skew<-skew[1:2,"predictor"]
	worst_skew<-skew[nrow(dist_stats)-1:nrow(dist_stats),"predictor"]
	median_skew<-skew[(nrow(dist_stats)+1)/2,"predictor"]
	kurt<-dist_stats[order(abs(dist_stats[,"kurtosis"])),]
	best_kurt<-kurt[1:2,"predictor"]
	worst_kurt<-kurt[nrow(dist_stats)-1:nrow(dist_stats),"predictor"]
	median_kurt<-kurt[(nrow(dist_stats)+1)/2,"predictor"]
	#best skew
	png(paste0("st01_densisty_plot_best_skewness_",best_skew,"_",sample,".png"))
  	plot_obj<-get_distribution_plot(data.frame(data[,best_skew]),config$zscore_raw)
  	plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",best_skew))
  	print(plot_obj)
  	dev.off()
  	#worst skew
	png(paste0("st01_densisty_plot_worst_skewness_",worst_skew,"_",sample,".png"))
  	plot_obj<-get_distribution_plot(data.frame(data[,worst_skew]),config$zscore_raw)
  	plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",worst_skew))
  	print(plot_obj)
  	dev.off()
  	#median skew
	png(paste0("st01_densisty_plot_median_skewness_",median_skew,"_",sample,".png"))
  	plot_obj<-get_distribution_plot(data.frame(data[,median_skew]),config$zscore_raw)
  	plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",median_skew))
  	print(plot_obj)
  	dev.off()
  	#best kurt
	png(paste0("st01_densisty_plot_best_kurtosis_",best_kurt,"_",sample,".png"))
  	plot_obj<-get_distribution_plot(data.frame(data[,best_kurt]),config$zscore_raw)
  	plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",best_kurt))
  	print(plot_obj)
  	dev.off()
  	#worst kurt
	png(paste0("st01_densisty_plot_worst_kurtosis_",worst_kurt,"_",sample,".png"))
  	plot_obj<-get_distribution_plot(data.frame(data[,worst_kurt]),config$zscore_raw)
  	plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",worst_kurt))
  	print(plot_obj)
  	dev.off()
  	#median kurt
	png(paste0("st01_densisty_plot_median_kurtosis_",median_kurt,"_",sample,".png"))
  	plot_obj<-get_distribution_plot(data.frame(data[,median_kurt]),config$zscore_raw)
  	plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",median_kurt))
  	print(plot_obj)
  	dev.off()
	####
  	pdf(paste0("st01_densisty_plot_",sample,".pdf"))
  	par(mfrow=c(3,4))
	for(i in 1:ncol(data)){
 	 	print(colnames(data)[i])
  		print(class(data[,i]))
  		plot_obj<-get_distribution_plot(data.frame(data[,i]),config$zscore_raw)
  		plot_obj<-plot_obj + ggtitle(paste0("Raw Distribution ",colnames(data)[i]))
  		print(plot_obj)
	}
  	dev.off()
	return(data)
}

remove_single_value_predictors<-function(df){
	new_df<-df[,which(apply(df, 2, function(x) length(unique(na.omit(x))) != 1))]
	new_df<-new_df[,which(apply(new_df, 2, function(x) length(unique(x)) != 1))]
	return(new_df)
}

make_pcs<-function(df){
	pcs<-prcomp(na.omit(df[,which(colnames(df)!="iid")]),center=TRUE,scale=TRUE)
	return(pcs)
}

make_screeplot<-function(df,sample=""){
	pcs<-make_pcs(df) #,na.action=na.omit()
	png(paste0("screeplot_",sample,".png"))
	screeplot(pcs)
	dev.off()
	#want to output the number of PCs that explain 95% of the variance in omics df
	temp<-summary(pcs)$importance
	temp2<-temp[,which(temp[3,]>=0.95),drop=F]
	pc_95<-colnames(temp2)
	write.table(pc_95,paste0("pc_95_variance_",sample,".txt"),col.names=T,row.names=T,quote=F,sep="\t")
}


make_correlation_heatmap<-function(df,sample=""){
	cor_matrix<-rcorr(as.matrix(df[-1])) #,use="pairwise.complete.obs"
	png(paste0("correlation_heatmap_",sample,".png"),width=800,height=800)
	my_palette <- colorRampPalette(c("red", "white", "blue"))(n=199)
	heatmap.2(cor_matrix$r,Rowv=TRUE,Colv=TRUE,scale ="none",trace="none",symkey=FALSE,col=my_palette)
	dev.off()
}

make_outcome_correlation_heatmap<-function(omics_df,cov_df,outcome,sample=""){
	outcome_df<-cov_df[,c("iid",outcome)]
	big_df<-merge(outcome_df,omics_df,by="iid",all=T)
	r<-apply(big_df[,2:ncol(big_df)],2,function(x) cor(as.numeric(big_df[,outcome]),as.numeric(x),use="pairwise.complete.obs"))
	r<-data.frame(r)
	df<- data.frame(predictors = row.names(r), r)
	df$colour<-ifelse(df$r>0,"Positive Correlation","Negative Correlation")
	df<-df[order(df$r),]
	df<-df[which(df[,"predictors"]!=outcome),]
	df$predictors<-factor(df$predictors,levels=df$predictors)
	df$colour<-factor(df$colour,levels=c("Positive Correlation","Negative Correlation"))
	ggplot(df,aes(x=predictors,y=r,fill=colour)) +
  		geom_bar(stat="identity") +
  		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  		scale_fill_manual(values=c("blue3","#cb181d"),labels = c("Positive Correlation","Negative Correlation"))
  	ggsave(paste0("outcome_correlation_plot_",sample,".png"))
}

make_outcome_predictor_plots<-function(omics_df,cov_df,outcome,sample=""){
	outcome_df<-cov_df[,c("iid",outcome)]
	big_df<-merge(outcome_df,omics_df,by="iid",all=T)
	predictors<-colnames(big_df)[!(colnames(big_df)%in%c("iid",outcome))]
	pdf(paste0("st01_",outcome,"_predictor_plots.pdf"))
	par(mfrow=c(3,4))
	for(predictor in predictors){
		writeLines(predictor)
		model<-lm(big_df[,outcome]~big_df[,predictor],na.action=na.omit)
		r2<-summary(model)$r.squared
		plot(big_df[,predictor],big_df[,outcome],
			main=paste0("r2 = ",round(r2,2)," ",predictor),
			ylab=outcome,
			xlab=predictor)
		abline(model)
	}
	dev.off()
}

make_batch_plot<-function(df){
	#read in vene date
	base<-fread("/exports/igmm/eddie/wilson-lab/data/base_data/orcades/phenotypes/orcades_base_phenotypes.tsv",header=T,stringsAsFactors=F,data.table=F)
	vene_date<-base[,c("iid","vene_date")]
	vene_date<-vene_date[order(as.Date(vene_date$vene_date, format="%d/%m/%Y")),]
	vene_date$join_order<-1:nrow(vene_date)
	big_df<-merge(vene_date[,c("iid","join_order")],df,by="iid",all=T)
	pdf("batch_plots.pdf")
	par(mfrow=c(3,4))
	for(i in 3:ncol(big_df)){
		p<-ggplot(big_df, aes(x=join_order, y=big_df[,i])) + 
			geom_point(alpha = 0.5) + 
      		geom_smooth() +
      		ggtitle(paste0(names(big_df[i])))
      	print(p)
	}
	dev.off()
}

make_profile_plots<-function(df,sample=""){
	long_df<-reshape2::melt(df,id.vars="iid")
	ggplot(long_df,aes(x=variable,y=value,colour=iid)) +
		geom_point(show.legend = FALSE) +
		xlab("Predictors") +
		ylab("Outcome") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		theme_classic()
	ggsave(paste0("st01_profile_plot_",sample,".png"),width=15,height=7)
}



scale_centre_vec<-function(vec){
	vec<-scale(vec,center=T)
	return(vec)
}

scale_centre<-function(df){
	df<-apply(df,2,scale_centre_vec)
	return(df)
}



correct_data<-function(omics_df,cov_df,covariate_list,sample="discovery"){
	df<-merge(omics_df,cov_df,by="iid")
	corr_data<-data.frame(matrix(NA,ncol=ncol(omics_df),nrow=nrow(df)))
	corr_data[,1]<-df[,"iid"]
	names(corr_data)[1]<-"iid"
	for(i in 2:ncol(omics_df)){
		print(i)
  		possibleError <- tryCatch(
      		model<-lm(df[,i] ~ .,data=subset(df,select=covariate_list),na.action=na.omit),
      		error=function(e) e
  		)
		if(!inherits(possibleError, "error")){
    		results<-summary(model)
    		print(summary(model))
			temp<-data.frame(results$residuals)
			corr_data[,i]<-resid(lm(df[,i] ~ .,data=subset(df,select=covariate_list),na.action=na.exclude))
			colnames(corr_data)[i]<-colnames(omics_df)[i]	
  		}else{
  			corr_data[,i]<-NA
			colnames(corr_data)[i]<-colnames(omics_df)[i]
  		}
	}
	return(list(corr_data=corr_data,covariates_chosen=covariate_list)) #covariates_chosen=new_cov_list
}

extract_pcs<-function(omics_df,n,sample=""){
	mat<-transpose_data(na.omit(omics_df))
	pcs<-prcomp(mat,center=TRUE,scale=TRUE)
	data<-data.frame(pcs$rotation)
	data[,"iid"]<-rownames(data)
	df<-data[,c("iid",paste0("PC",seq(1:n)))]
	write.table(df,paste0("omics_data_",n,"_PCs_",sample,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
}

get_testing_from_overlap<-function(df){
	overlap_iids<-read.table("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/overlap_iids.tsv",header=T,stringsAsFactors=F)
	
}

get_testing_from_overlap<-function(){}

split_training_testing<-function(df,iteration){	
	total<-nrow(df)
	no_testing<-round(0.25*total,digits=0)
	#need to read in the overlap iids
	overlap_iids<-read.table("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/overlap_iids.tsv",header=T,stringsAsFactors=F)
	if(iteration=="core"){
		x<-read.table("../p03_make_clock/st03_testing_data_1.tsv",header=T,stringsAsFactors=F)
		overlap_testing_iids<-x[,"iid",drop=F]
		names(overlap_testing_iids)<-"x"
	}else if(iteration=="pcs"){
		x<-read.table("../../p03_make_clock/st03_testing_data_1.tsv",header=T,stringsAsFactors=F)
		overlap_testing_iids<-x[,"iid",drop=F]
		names(overlap_testing_iids)<-"x"
	}else if(iteration=="b_only"){
		x<-read.table("../../p03_make_clock/st03_testing_data_1.tsv",header=T,stringsAsFactors=F)
		overlap_testing_iids<-x[,"iid",drop=F]
		names(overlap_testing_iids)<-"x"
	}else{
		overlap_testing_iids<-read.table("/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/overlap_testing_iids.tsv",header=T,stringsAsFactors=F)
	}
	#need to find out how many extra from the overlap needed
	data_testing<-subset(df,(df$iid %in% overlap_testing_iids$x))
	extra<-no_testing-nrow(data_testing)
	if(extra>0){
		# get a random 25% of the overlap as testing data
		overlap_sample_not_yet_in_testing<-subset(df,((df$iid %in% overlap_iids$x) & !(df$iid%in%overlap_testing_iids$x)))
		extra_testing<-overlap_sample_not_yet_in_testing[sample.int(nrow(overlap_sample_not_yet_in_testing),extra),]
		data_testing<-rbind(data_testing,extra_testing)
		dim(data_testing)
		write.table(data_testing,paste0("st03_testing_data_",iteration,".tsv"),col.names = T,row.names = F,quote = F,sep="\t")
	}else{
		data_testing<-subset(df,(df$iid %in% overlap_testing_iids$x))
		data_testing<-data_testing[sample.int(nrow(data_testing),no_testing),]
		dim(data_testing)
		write.table(data_testing,paste0("st03_testing_data_",iteration,".tsv"),col.names = T,row.names = F,quote = F,sep="\t")
	}
	#need to get the oter 75% as the training data
	data_training<-subset(df,!(df$iid %in% data_testing$iid))
	dim(data_training)
	write.table(data_training,paste0("st03_training_data_",iteration,".tsv"),col.names = T,row.names = F,quote = F,sep="\t")
	return(list(training_data=data_training,testing_data=data_testing))
}

outcome_matrix<-function(df,outcome){
	x<-as.matrix(df[,outcome])
	return(x)
}

predictor_matrix<-function(df,predictors){
	x<-as.matrix(df[,predictors])
	return(x)
}

split_predictor_outcome_matrices<-function(df,outcome,predictors){
	outcome_matrix<-outcome_matrix(df,outcome)
	predictor_matrix<-predictor_matrix(df,predictors)
	return(list(outcome_matrix=outcome_matrix,predictor_matrix=predictor_matrix))
}

cross_validation<-function(omics_data,outcome_data){
	cv<-cv.glmnet(omics_data,outcome_data,nfolds = 10,family="gaussian")
	return(cv)
}

fit_elastic_net<-function(omics_data,outcome_data){
	model<-glmnet(omics_data,outcome_data,family = "gaussian",alpha=0.5)
	return(model)
}

predict_outcome<-function(model,prediction_sample_predictors,lambda){
	predicted_outcome<-predict(model,prediction_sample_predictors,type="response",s=lambda)
	return(predicted_outcome)
}

get_model_coefficients<-function(model,lambda){
	coeffs<-coef(model, s =lambda)
	coeffs<-as.matrix(coeffs)
	coeffs<-cbind(rownames(coeffs), coeffs)
	return(coeffs)
}

write_out_model<-function(model,lambda,iteration){
	coeffs<-get_model_coefficients(model,lambda)
	write.table(coeffs,paste0("st03_model_coefficients_",iteration,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
}


plot_pred_obs<-function(pred_outcome,observed_outcome,sample){
	if(sample=="testing"){
		colour<-"blue"
	}else if(sample=="training"){
		colour<-"red"
	}else{
		colour<-"black"
	}
	x<-cor.test(pred_outcome,observed_outcome)
	plot_obj<-plot(pred_outcome,observed_outcome,
	     col=colour,
	     main=paste("r =",round(x$estimate, 4),", p<",signif(x$p.value, 3)),
	     xlab = "Predicted Outcome",
	     ylab = "Observed Outcome")
	return(plot_obj)
}

create_pred_obs_plots<-function(pred_outcome,observed_outcome,sample,iteration){
	png(paste0("st03_predicted_vs_observed_outcome_",sample,"_",iteration,".png"))
	plot_obj<-plot_pred_obs(pred_outcome,observed_outcome,sample)
	print(plot_obj)
	dev.off()
}

make_pred_obs_resids<-function(observed_outcome,pred_outcome){
	error_data<-cbind(observed_outcome,pred_outcome)
	error_data<-data.frame(error_data)
	names(error_data)<-c("observed_outcome","pred_outcome")
	error_data['resid']<-error_data$pred_outcome-error_data$observed_outcome
	return(error_data)
}

write_out_model_errors<-function(observed_outcome,pred_outcome,sample,iteration){
	data<-make_pred_obs_resids(observed_outcome,pred_outcome)
	write.table(data,paste0("st03_pred_obs_resid_",sample,"_",iteration,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
	#want to make scatter plot
	ggplot(data,aes(y=observed_outcome,x=resid)) +
		geom_point()
	ggsave(paste0("st03_scatter_obs_resid_",sample,"_",iteration,".png"))
}


write_out_model_errors_eddie<-function(observed_outcome,pred_outcome,sample,iteration){
	data<-make_pred_obs_resids(observed_outcome,pred_outcome)
	write.table(data,paste0("st03_pred_obs_resid_",sample,"_",iteration,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
	#want to make scatter plot
	ggplot(data,aes(y=observed_outcome,x=resid)) +
		geom_point()
	ggsave(paste0("st03_scatter_obs_resid_",sample,"_",iteration,".png"))
}

make_res_plots<-function(sample="orcades"){
	df<-data.frame(status=character(),resid=numeric())
	data<-fread(paste0("st03_pred_obs_resid_",sample,"_curr_smokers_1.tsv"),data.table=F)
	df_c<-data.frame(status="Current Smoker",resid=data$resid)
	df<-rbind(df,df_c)
	data<-fread(paste0("st03_pred_obs_resid_",sample,"_ex_smokers_1.tsv"),data.table=F)
	df_e<-data.frame(status="Ex Smoker",resid=data$resid)
	df<-rbind(df,df_e)
	data<-fread(paste0("st03_pred_obs_resid_",sample,"_non_smokers_1.tsv"),data.table=F)
	df_n<-data.frame(status="Non Smoker",resid=data$resid)
	df<-rbind(df,df_n)
	head(df)
	dim(df)
	#need to merge the res tables
	ggplot(df,aes(x=resid,fill=status)) +
		geom_density(alpha=0.6) +
		xlim(range(df$resid))
	ggsave(paste0("resid_by_status",sample,".png"))
}

make_scatter_plots<-function(sample="orcades"){
	df<-data.frame(status=character(),resid=numeric(),observed_outcome=numeric())
	data<-fread(paste0("st03_pred_obs_resid_",sample,"_curr_smokers_1.tsv"),data.table=F)
	df_c<-data.frame(status="Current Smoker",resid=data$resid,observed_outcome=data$observed_outcome)
	df<-rbind(df,df_c)
	data<-fread(paste0("st03_pred_obs_resid_",sample,"_ex_smokers_1.tsv"),data.table=F)
	df_e<-data.frame(status="Ex Smoker",resid=data$resid,observed_outcome=data$observed_outcome)
	df<-rbind(df,df_e)
	data<-fread(paste0("st03_pred_obs_resid_",sample,"_non_smokers_1.tsv"),data.table=F)
	df_n<-data.frame(status="Non Smoker",resid=data$resid,observed_outcome=data$observed_outcome)
	df<-rbind(df,df_n)
	head(df)
	dim(df)
	#need to merge the res tables
	ggplot(df,aes(y=observed_outcome,x=resid,colour=status)) +
		geom_point()
	ggsave(paste0("scatter_resid_by_status",sample,".png"))
}


plot_cv_curve<-function(cv,iteration){
	png(paste0("st03_cv_curve_",iteration,".png"),width=700,height=400)
	plot(cv)
	dev.off()
}

plot_effect_distribution<-function(iteration){
	coeffs<-fread(paste0("st03_model_coefficients_",iteration,".tsv"),data.table=F,header=T)
	coeffs<-coeffs[-1,]
	ggplot(coeffs, aes(coeffs[,2])) + 
      geom_density(alpha = 0.5,fill="#1c9099",colour="grey") +
      theme_light(base_size = 14) +
      xlab("")
    ggsave(paste0("st03_effect_size_distribution_plot",iteration,".png"))
}


penalised_regression<-function(training_predictors,training_outcome,testing_predictors,testing_outcome,outcome,method,if_elastic_net="cv",iteration=1){
	#sort which method we are doing
	if((method%in%c("lasso","ridge")) | (if_elastic_net!="cv")){
		writeLines("Performing cross validation in the training sample...")
		training_data_cv<-cross_validation(training_predictors,training_outcome)
		writeLines("Plotting cross validation curve...")
		plot_cv_curve(training_data_cv,iteration)
		lambda_for_model<-training_data_cv$lambda.min
		writeLines(paste0("Lambda selected for penalised regression = ",lambda_for_model))
		if(method=="lasso"){
			writeLines("Fitting LASSO Regression...")
			model<-glmnet(training_predictors,training_outcome,family = "gaussian",alpha=1)
		}else if(method=="ridge"){
			writeLines("Fitting Ridge Regression...")
			model<-glmnet(training_predictors,training_outcome,family = "gaussian",alpha=0)
		}else if((method=="elastic_net") & (if_elastic_net!="cv")){
			writeLines("Fitting Elastic Net Regression...")
			model<-glmnet(training_predictors,training_outcome,family = "gaussian",alpha=if_elastic_net)
		}
		writeLines("Writing out model coefficients...")
		write_out_model(model,lambda_for_model,iteration)
		writeLines("Plotting effect size distribution...")
		plot_effect_distribution(iteration)
		writeLines("Predicting outcome phenotype in training set...")
		predicted_outcome_training<-predict_outcome(model,training_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the training set...")
		create_pred_obs_plots(predicted_outcome_training,training_outcome,"training",iteration)
		writeLines("Writing out model residuals for the training set...")
		write_out_model_errors(training_outcome,predicted_outcome_training,"training",iteration)
		writeLines("Predicting outcome phenotype in testing set...")
		predicted_outcome_testing<-predict_outcome(model,testing_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the testing set...")
		create_pred_obs_plots(predicted_outcome_testing,testing_outcome,"testing",iteration)
		writeLines("Writing out model residuals for the testing set...")
		write_out_model_errors(testing_outcome,predicted_outcome_testing,"testing",iteration)
	}else if((method=="elastic_net") & (if_elastic_net=="cv")){
		writeLines("Fitting Elastic Net Regression...")
		x<-cbind(data.frame(training_outcome),data.frame(training_predictors))
		names(x)[1]<-outcome
		set.seed(123)
		model <- caret::train(age_at_vene ~ .,data=x, method = "glmnet",trControl = trainControl("cv", number=10),tuneLength=100)
		lambda_for_model<-model$bestTune$lambda
		writeLines("Writing out model coefficients...")
		write_out_model(model$finalModel,model$bestTune$lambda,iteration)
		writeLines("Predicting outcome phenotype in training set...")
		writeLines("Plotting effect size distribution...")
		plot_effect_distribution(iteration)
		writeLines("Predicting outcome phenotype in training set...")
		predicted_outcome_training<-predict_outcome(model$finalModel,training_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the training set...")
		create_pred_obs_plots(predicted_outcome_training,training_outcome,"training",iteration)
		writeLines("Writing out model residuals for the training set...")
		write_out_model_errors(training_outcome,predicted_outcome_training,"training",iteration)
		writeLines("Predicting outcome phenotype in testing set...")
		predicted_outcome_testing<-predict_outcome(model$finalModel,testing_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the testing set...")
		create_pred_obs_plots(predicted_outcome_testing,testing_outcome,"testing",iteration)
		writeLines("Writing out model residuals for the testing set...")
		write_out_model_errors(testing_outcome,predicted_outcome_testing,"testing",iteration)
	}
}


penalised_regression_eddie<-function(training_predictors,training_outcome,testing_predictors,testing_outcome,outcome,method,if_elastic_net="cv",iteration=1){
	#sort which method we are doing
	if((method%in%c("lasso","ridge")) | (if_elastic_net!="cv")){
		writeLines("Performing cross validation in the training sample...")
		training_data_cv<-cross_validation(training_predictors,training_outcome)
		writeLines("Plotting cross validation curve...")
		plot_cv_curve(training_data_cv,iteration)
		lambda_for_model<-training_data_cv$lambda.min
		writeLines(paste0("Lambda selected for penalised regression = ",lambda_for_model))
		if(method=="lasso"){
			writeLines("Fitting LASSO Regression...")
			model<-glmnet(training_predictors,training_outcome,family = "gaussian",alpha=1)
		}else if(method=="ridge"){
			writeLines("Fitting Ridge Regression...")
			model<-glmnet(training_predictors,training_outcome,family = "gaussian",alpha=0)
		}else if((method=="elastic_net") & (if_elastic_net!="cv")){
			writeLines("Fitting Elastic Net Regression...")
			model<-glmnet(training_predictors,training_outcome,family = "gaussian",alpha=if_elastic_net)
		}
		writeLines("Writing out model coefficients...")
		write_out_model(model,lambda_for_model,iteration)
		writeLines("Plotting effect size distribution...")
		plot_effect_distribution(iteration)
		writeLines("Predicting outcome phenotype in training set...")
		predicted_outcome_training<-predict_outcome(model,training_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the training set...")
		create_pred_obs_plots(predicted_outcome_training,training_outcome,"training",iteration)
		writeLines("Writing out model residuals for the training set...")
		write_out_model_errors_eddie(training_outcome,predicted_outcome_training,"training",iteration)
		writeLines("Predicting outcome phenotype in testing set...")
		predicted_outcome_testing<-predict_outcome(model,testing_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the testing set...")
		create_pred_obs_plots(predicted_outcome_testing,testing_outcome,"testing",iteration)
		writeLines("Writing out model residuals for the testing set...")
		write_out_model_errors_eddie(testing_outcome,predicted_outcome_testing,"testing",iteration)
	}else if((method=="elastic_net") & (if_elastic_net=="cv")){
		writeLines("Fitting Elastic Net Regression...")
		x<-cbind(data.frame(training_outcome),data.frame(training_predictors))
		names(x)[1]<-outcome
		set.seed(123)
		model <- caret::train(age_at_vene ~ .,data=x, method = "glmnet",trControl = trainControl("cv", number=10),tuneLength=100)
		lambda_for_model<-model$bestTune$lambda
		writeLines("Writing out model coefficients...")
		write_out_model(model$finalModel,model$bestTune$lambda,iteration)
		writeLines("Predicting outcome phenotype in training set...")
		writeLines("Plotting effect size distribution...")
		plot_effect_distribution(iteration)
		writeLines("Predicting outcome phenotype in training set...")
		predicted_outcome_training<-predict_outcome(model$finalModel,training_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the training set...")
		create_pred_obs_plots(predicted_outcome_training,training_outcome,"training",iteration)
		writeLines("Writing out model residuals for the training set...")
		write_out_model_errors_eddie(training_outcome,predicted_outcome_training,"training",iteration)
		writeLines("Predicting outcome phenotype in testing set...")
		predicted_outcome_testing<-predict_outcome(model$finalModel,testing_predictors,lambda_for_model)
		writeLines("Plotting observed versus predicted outcome in the testing set...")
		create_pred_obs_plots(predicted_outcome_testing,testing_outcome,"testing",iteration)
		writeLines("Writing out model residuals for the testing set...")
		write_out_model_errors_eddie(testing_outcome,predicted_outcome_testing,"testing",iteration)
	}
}

n_variables_selected<-function(iteration){
	data<-read.table(paste0("st03_model_coefficients_",iteration,".tsv"),header=T,stringsAsFactors = F)
	return(paste0(sum(data[,1]!=0)," out of ",nrow(data)-1, " were selection for model inclusion\n"))
}


transpose_data<-function(df){
	rownames(df)<-df[,"iid"]
	tdf<-t(df[-1])
	return(tdf)
}

do_pca<-function(df,outcome_df,sample=""){
	mat<-transpose_data(na.omit(df))
	pcs<-prcomp(mat,center=TRUE,scale=TRUE)
	png(paste0("pcs_variance_",sample,".png"))
	plot(pcs,type="l")
	dev.off()
	head(summary(pcs))
	biplot(pcs)
	#
	data<-data.frame(pcs$rotation)
	data[,"iid"]<-rownames(data)
	outcome_df<-outcome_df[,c("iid","resid")]
	for_plot<-merge(data,outcome_df,by="iid")
	plot(for_plot$PC1,for_plot$PC2)
	ggplot(for_plot,aes(x=PC1,y=PC2,colour=resid)) +
		geom_point() +
		theme_bw()
	ggsave(paste0("pca_plot_",sample,".png"))
	#hierarcical clustering
	cor_matrix<-rcorr(as.matrix(mat)) #,use="pairwise.complete.obs"
	png(paste0("hclust_heatmap_",sample,".png"),width=800,height=800)
	my_palette <- colorRampPalette(c("red", "white", "blue"))(n=199)
	heatmap.2(cor_matrix$r,Rowv=TRUE,Colv=TRUE,scale ="none",trace="none",symkey=FALSE,col=my_palette)
	dev.off()
}

get_model_errors<-function(sample){
	errors<-read.table(paste0("../p04_clock_500_iterations/st03_pred_obs_resid_",sample,"_1.tsv"),header=T,stringsAsFactors = F)
	errors<-data.frame(model_1=errors$resid)
	model_2<-read.table(paste0("../p04_clock_500_iterations/st03_pred_obs_resid_",sample,"_2.tsv"),header=T,stringsAsFactors = F)
	model_2<-data.frame(model_2=model_2$resid)
	errors<-cbind(errors,model_2)
	#now read in data from models 3:500
	for(i in 3:500){
  		model_i<-read.table(paste0("../p04_clock_500_iterations/st03_pred_obs_resid_",sample,"_1.tsv"),header=T,stringsAsFactors = F)
  		cat(paste0("Successfuly read in errors from model ",i,"\n"))
  		model_i<-data.frame(temp=model_i$resid)
  		names(model_i)<-paste0("model_",i)
  		errors<-cbind(errors,model_i)
	}
	return(errors)
}

get_model_errors_single_variable<-function(sample){
	errors<-read.table(paste0("../p04_clock_500_iterations/st03_pred_obs_resid_",sample,"_1.tsv"),header=T,stringsAsFactors = F)
	model_2<-read.table(paste0("../p04_clock_500_iterations/st03_pred_obs_resid_",sample,"_2.tsv"),header=T,stringsAsFactors = F)
	errors<-rbind(errors,model_2)
	#now read in data from models 3:500
	for(i in 3:500){
  		model_i<-read.table(paste0("../p04_clock_500_iterations/st03_pred_obs_resid_",sample,"_1.tsv"),header=T,stringsAsFactors = F)
  		cat(paste0("Successfuly read in errors from model ",i,"\n"))
  		errors<-rbind(errors,model_i)
	}
	return(errors)
}


get_stats<-function(resid_file,study,omic){
	df<-fread(resid_file,data.table=F)
	x<-cor.test(as.numeric(df$observed_outcome),as.numeric(df$pred_outcome))
	r<-x$estimate
	lower<-x$conf.int[1]
	upper<-x$conf.int[2]
	p<-x$p.value
	return(data.frame(omic=omic,study=study,r=r,lower=lower,upper=upper,p=p))
}


slope_extract<-function(resid_file,study,omic){
	data<-fread(resid_file,data.table=F)
	model<-lm(as.numeric(data$observed_outcome)~as.numeric(data$pred_outcome))
	x<-summary(model)
	beta<-x$coefficients[2,1]
	se<-x$coefficients[2,2]
	p<-x$coefficients[2,4]
	r2<-x$r.squared
	r2_adj<-x$adj.r.squared
	return(data.frame(omic=omic,study=study,beta=beta,se=se,p=p,r2=r2,r2_adj=r2_adj))
}

