# R script to calculate Neffective size in Purves 2019 anxiety sumstats (maybe also neu)
# see here: https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics

#load in required packages
require(data.table)
require(dplyr)

# set wd
setwd("/exports/eddie/scratch/s2421111/gwas")

# read in Anxiety summary stats from Purves et al. 2019
ANX <- fread("META_UKBB_iPSYCH_anx.sumstats", data.table = FALSE)

# read in 1000 genomes reference file used to get approximation of SNP MAF as MAF not present in the anxiety sumstats
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging")
ref <- fread("reference.1000G.maf.0.005.txt", data.table = FALSE)

# subset reference file to just SNP and MAF
ref <- ref[, c("SNP", "MAF")]

# convert SNP column in ANX to a factor variable with the same levels as ref$SNP
ANX$SNP <- factor(ANX$SNP, levels = ref$SNP)

# merge ANX and reference file
ANX <- inner_join(ANX, ref, by = "SNP", all = FALSE)

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
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS")
write.table(ANX, file = "ANX_withNeff.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)