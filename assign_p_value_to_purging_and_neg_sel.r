

data<-read.table("/proj/snic2020-6-127/private/Analyses/STATS/ALL_AUTOSOMES/Purging_and_negative_selection_all_breeds_outgroup_andean_jackknifewindows_50_sites_SNPs_and_INDELs_pp0_pp7_autosames_210112.txt", header=TRUE, sep="\t")

z<-abs(1-data$weighted.bias.corrected.jackknife.estimator.of.R)/data$weighted.SD

pvalue2sided=2*(1-pnorm(abs(z)))

bonferoni_threshold=0.05/((nrow(data)-1)/2)

threshold<-c()
for (tal in (1:nrow(data))) {
threshold[tal]<-c(bonferoni_threshold)
}

complete_data<-c()
complete_data<-cbind(data, pvalue2sided, threshold)
cdv1<-sapply(complete_data[,1], as.character)
cdv2<-sapply(complete_data[,2], as.character)
cdv3<-sapply(complete_data[,3], as.character)
cdv.tot<-cbind(cdv1,cdv2,cdv3,complete_data[,4:21])

write.table(cdv.tot,"Purging_and_negative_selection_all_breeds_outgroup_andean_jackknifewindows_50_sites_SNPs_and_INDELs_pp0_pp7_autosames_210112_with_pvalues.txt", sep = "\t")


#cd<-as.data.frame(do.call(cbind, as.character(complete_data)))

#s<-c()
#if (pvalue2sided[tal]<bonferoni_threshold)
#s<-c(paste(data$breed_1[tal], data$breed_2[tal], data$Mutation[tal]))
#print (s)
#}
