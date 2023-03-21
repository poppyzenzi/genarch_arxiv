#load in required packages
require(data.table)
require(dplyr)

# set wd
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS")

#read in Anxiety summary stats from Purves et al. 2019
ANX<-gzfile("TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz",data.table=FALSE)
# or try
ANX<-read.table("TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz", header = TRUE)


#read in 1000 genomes reference file used to get approximation of SNP MAF
#as MAF not present in the anxiety summary statistics file
ref<-fread("reference.1000G.maf.0.005.txt",data.table=FALSE)

#subset reference file to just SNP and MAF
attach(ref)
ref<-data.frame(SNP,MAF)

#merge Anxiety and reference file
ANX<-inner_join(ANX,ref,by="SNP",all=F)

#calculate effective sample size implied by GWAS summary statistics
ANX$Neff<-4/((2*ANX$MAF*(1-ANX$MAF))*ANX$StdErr^2)

#calculate total effective N to cap backed out Neff
Ncases<-31977
Ncontrols<-82114
v<-Ncases/(Ncases+Ncontrols)
TotalNeff<-4*v*(1-v)*(Ncases+Ncontrols)

#cap at 1.1 of total effective N
ANX$Neff<-ifelse(ANX$Neff > 1.1*TotalNeff, 1.1*TotalNeff, ANX$Neff)

#lower limit of 0.5 of total effective N
ANX$Neff<-ifelse(ANX$Neff < 0.5*TotalNeff, 0.5*TotalNeff, ANX$Neff)

#remove reference panel MAF from file
ANX$MAF<-NULL

#remove sample size column so munge knows to use Neff column
ANX$TotalSampleSize<-NULL

#output the summary stats with the SNP-specific effective sample size column
write.table(ANX, file = "ANX_withNeff.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)