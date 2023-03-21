#this is an R script (on EDDIE) to munge files for step of GSEM

#set working directory
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS")
#setwd("/exports/eddie/scratch/s2421111/gwas")

#create vector of the summary statistics files
# should be 7 base GWAS to munge
files<-c("daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz","pgc_bip_qcd.txt","pgc_scz_qcd.txt",
           "iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz","mdd3_ss_MAF.txt",
           "TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz", "luciano_neur_ss.txt",

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging")

#define the reference file being used to align alleles across summary stats
#here we are using hapmap3
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging")
hm3<-"eur_w_ld_chr/w_hm3.snplist"

#name the traits
trait.names<-c("ADHD","BIP","SCZ", "ASD", "MDD", "ANX", "NEU")

#list the sample sizes. for binary traits from meta-analysed multiple cohorts,
# this should reflect the sum of effective sample sizes across contributing cohorts or the sample size col
# should reflect the SNP specific, sum of effective sample sizes
# where sumstats have SNP-specific cols then put NA and it will use this
# When inputting the sum of effective sample sizes, the sample prevalence should then be entered as 0.5
# when running ldsc to reflect the fact that effective sample size already corrects for sample ascertainment.

# sample size calculation needed for ASD, ANX and NEU
N=c(NA,NA,NA,ASD ,NA, ANX , NEU)

#define imputation and MAF filters
info.filter=0.8
maf.filter=0.01

#munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)