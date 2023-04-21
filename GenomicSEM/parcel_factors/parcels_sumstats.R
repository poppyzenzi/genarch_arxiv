# step 4: prep sum stats for GWAS
require(GenomicSEM)
library(GenomicSEM)

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

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/")

# =====================================================================================================================

LDSCoutput = "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/parcel_factors/LDSCoutput_mood.RData"

#run the multivariate GWAS using parallel processing
MOOD_parcel <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = MOOD_sumstats, estimation = "DWLS",
                                    cores = NULL, toler = FALSE, SNPSE = FALSE, parallel = TRUE,
                                    Output = NULL,GC="standard",MPI=FALSE)


# note that the code written above specifies all arguments for completeness, but as many of these arguments
# are set to the package default it could also be written as below and produce the same result:
MOOD_parcel  <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = MOOD_sumstats)

save(MOOD_parcel, file="MOOD_GWAS.txt")




