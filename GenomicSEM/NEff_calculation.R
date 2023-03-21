# to calculate effective sample size for binary meta sum stats

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS")


# here for ANX, ASD,
##read in information about sample size per cohort
ASD<-read.csv("PTSDcohorts.csv",header=TRUE)
#ANX<-read.csv("TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz", header=TRUE) # not cc
#NEU<-read.csv("luciano_neur_ss.txt", header=TRUE) # not cc

#calculate sample prevalence for each cohort
ASD$v<-ASD$cases/(ASD$cases+ASD$controls)

#calculate cohort specific effective sample size
ASD$EffN<-4*ASD$v*(1-ASD$v)*(ASD$cases+ASD$controls)

#calculate sum of effective sample size: 5831.346
sum(PTSD$EffN)