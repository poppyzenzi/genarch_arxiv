# step 4: prep sum stats for GWAS for 3 parcels
require(GenomicSEM)
library(GenomicSEM)

# ======================= mood ============================
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")
files<-c("mdd3_ss_MAF.txt", "luc_neur_ss.txt", "ANX_withNeff.txt")
ref= "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/reference.1000G.maf.0.005.txt"
trait.names<-c('MDD', 'NEU', 'ANX')
se.logit=c(T,T,T)
N=c(NA,329821,NA)

MOOD_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,
                            se.logit=se.logit,OLS=NULL,linprob=NULL,N=N,
                            betas=NULL,info.filter=0.8,maf.filter=maf.0.01,
                            keep.indel=FALSE,parallel=TRUE,cores=NULL)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/parcel_factors/sumstats/")
saveRDS(MOOD_sumstats, file="MOOD_sumstats.rds")

# ======================= neurodev ============================

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")
files<-c("daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz", "iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz")
ref= "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/reference.1000G.maf.0.005.txt"
trait.names<-c('ADHD', 'ASD')
se.logit=c(T,T)
N=c(NA,43777.81914)

NEURODEV_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,
                            se.logit=se.logit,OLS=NULL,linprob=NULL,N=N,
                            betas=NULL,info.filter=0.8,maf.filter=maf.0.01,
                            keep.indel=FALSE,parallel=TRUE,cores=NULL)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/parcel_factors/sumstats/")
saveRDS(NEURODEV_sumstats, file="NEURODEV_sumstats.rds")

# ======================= psychotic ============================

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")
files<-c("pgc_bip_qcd.txt","pgc_scz_qcd.txt")
ref= "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/reference.1000G.maf.0.005.txt"
trait.names<-c('BIP', 'SCZ')
se.logit=c(T,T)
N=c(NA,NA)

PSYCHOTIC_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,
                            se.logit=se.logit,OLS=NULL,linprob=NULL,N=N,
                            betas=NULL,info.filter=0.8,maf.filter=maf.0.01,
                            keep.indel=FALSE,parallel=TRUE,cores=NULL)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/parcel_factors/sumstats/")
saveRDS(PSYCHOTIC_sumstats, file="PSYCHOTIC_sumstats.rds")




