
# R on eddie for prepping GWAS for common factor GWAS GSEM analysis

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS")

files<-c("daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz","pgc_bip_qcd.txt","pgc_scz_qcd.txt",
           "iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz","mdd3_ss_MAF.txt", "luc_neur_ss.txt",
           "ANX_withNeff.txt")


ref= "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/reference.1000G.maf.0.005.txt"


trait.names<-c("ADHD","BIP","SCZ", "ASD", "MDD", "NEU", "ANX")

# whether SE is logsitic - check README file
se.logit=c(T,T,T,F)

# Whether the phenotype was a dichotomous outcome for which there are only Z-statistics in the
# summary statistics file -or- it was a dichotomous outcome analyzed using an OLS estimator
linprob=c(F,F,F,T)


info.filter=.8
maf.filter=0.01


#Whether the phenotype was a
# continuous outcome analyzed using an observed least square (OLS; i.e., linear)
OLS=NULL

N=c(NA,NA,NA,43777.81914,NA,329821,NA)

betas=NULL

PSYCH_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,
                            se.logit=se.logit,OLS=NULL,linprob=linprob,N=NULL,
                            betas=NULL,info.filter=info.filter,maf.filter=maf.filter,
                            keep.indel=FALSE,parallel=FALSE,cores=NULL)


###### common factor GWAS (combine sumstats and LDSC) #######

#run the multivariate GWAS using parallel processing
LDSCoutput = "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/LDSCoutput"

PSYCH_factor <- commonfactorGWAS(covstruc = LDSCoutput,
                                    SNPs = PSYCH_sumstats, estimation = "DWLS", cores = NULL,
                                    toler = FALSE, SNPSE = FALSE, parallel = TRUE, Output = NULL,
                                    GC="standard",MPI=FALSE)


# calculating sample size for factors

#restrict to MAF of 40% and 10%
#PSYCH_factor2<-subset(PSYCH_factor, PSYCH_factor$MAF <= .4 & PSYCH_factor$MAF >= .1)
#calculate expected sample size (N_hat)
#N_hat<-mean(1/((2*PSYCH_factor2$MAF*(1-PSYCH_factor2$MAF))*PSYCH_factor2$SE^2))