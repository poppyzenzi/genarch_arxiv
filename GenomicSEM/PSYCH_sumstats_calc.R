# step 3: prep sum stats for GWAS
# model with SNP effects = multivariate GWAS sumstats

require(GenomicSEM)
library(GenomicSEM)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")


files<-c("daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz","pgc_bip_qcd.txt","pgc_scz_qcd.txt",
           "iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz","mdd3_ss_MAF.txt",
           "luc_neur_ss.txt", "ANX_withNeff.txt")

ref= "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/reference.1000G.maf.0.005.txt"

trait.names<-c('ADHD', 'BIP', 'SCZ', 'ASD','MDD', 'NEU', 'ANX')

# whether SE is logistic - check README file
se.logit=c(T,T,T,T,T,T,T)

# Whether the phenotype was a dichotomous outcome for which there are *only Z-statistics* in the
# summary statistics file -or- it was a dichotomous outcome analyzed using an OLS estimator as in UKB for some using hail software
linprob=NULL
# Whether the phenotype was a continuous outcome analyzed using an observed least square (OLS; i.e., linear)
OLS=NULL

info.filter=0.8
maf.filter=0.01

# if ss file has sample size col then can put NA, make sure NEff for dichiotomous traits

N=c(NA,NA,NA,43777.81914,NA,329821,NA)

# cont trait, if MAF col not available
betas=NULL

# sum stats -> use to generate p-factor PRS

PSYCH_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,
                            se.logit=se.logit,OLS=OLS,linprob=linprob,N=N,
                            betas=betas,info.filter=info.filter,maf.filter=maf.filter,
                            keep.indel=FALSE,parallel=FALSE,cores=NULL)


setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/")

save(PSYCH_sumstats, file="PSYCH_sumstats.txt")

# =====================================================================================================================

# this gets the common factor loading of each SNP


LDSCoutput = "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/LDSCoutput.RData"

# run the multivariate GWAS using parallel processing
PSYCH_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats, estimation = "DWLS",
                                    cores = NULL, toler = FALSE, SNPSE = FALSE, parallel = TRUE,
                                    Output = NULL, GC="standard", MPI=FALSE)

save(PSYCH_sumstats, file="PSYCH_sumstats.txt")

# note that the code written above specifies all arguments for completeness, but as many of these arguments
# are set to the package default it could also be written as below and produce the same result:
# PSYCH_factor  <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats)


# calculating sample size for factors (for PRS software)

#restrict to MAF of 40% and 10%
#PSYCH_factor2<-subset(PSYCH_factor, PSYCH_factor$MAF <= .4 & PSYCH_factor$MAF >= .1)
#calculate expected sample size (N_hat)
#N_hat<-mean(1/((2*PSYCH_factor2$MAF*(1-PSYCH_factor2$MAF))*PSYCH_factor2$SE^2))