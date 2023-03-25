
#this is an R script (on EDDIE) to munge files for step of GSEM

library(GenomicSEM)

#set working directory
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")
#setwd("/exports/eddie/scratch/s2421111/gwas")

#create vector of the summary statistics files
# should be 7 base GWAS to munge
files<-c("daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz","pgc_bip_qcd.txt","pgc_scz_qcd.txt",
           "iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz","mdd3_ss_MAF.txt", "luc_neur_ss.txt",
           "ANX_withNeff.txt")

# "TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz", becomes ANXwithNeff from R script calculation

#define the reference file being used to align alleles across summary stats
#here we are using hapmap3
hm3<-"/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/eur_w_ld_chr/w_hm3.snplist"

#name the traits
trait.names<-c("ADHD","BIP","SCZ", "ASD", "MDD", "NEU", "ANX")

#list the sample sizes. for binary traits from meta-analysed multiple cohorts, [see wiki 2.1]
# sum of effective sample sizes across contributing cohorts
# or the sample size col = SNP specific, sum of effective sample sizes
# where sumstats have SNP-specific cols then put NA and it will use this
# When inputting the sum of effective sample sizes, the sample prevalence should then be entered as 0.5
# when running ldsc to reflect the fact that effective sample size already corrects for sample ascertainment.

# ASD is calculated separately, NEU is continuous trait, ANX col was calculated as per section 2.1
N=c(NA,NA,NA,43777.81914,NA,329821,NA)

#define imputation and MAF filters
info.filter=0.8
maf.filter=0.01

#munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)