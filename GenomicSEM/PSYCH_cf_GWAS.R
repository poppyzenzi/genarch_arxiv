# step 4 : combine sumstats and ldsc output and run common factor GWAS
# this gets the common factor loading of each SNP

require(GenomicSEM)
library(GenomicSEM)
library(readr)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/")

# run the multivariate GWAS using parallel processing
LDSCoutput <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/LDSCoutput.rds")
PSYCH_sumstats <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/PSYCH_sumstats.rds")

PSYCH_factor  <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats)

save(PSYCH_factor, file="PSYCH_factor.txt")
saveRDS(PSYCH_factor, file="PSYCH_factor.rds")




# ==================================================================================================
# with all args
#PSYCH_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats, estimation = "DWLS",
#                                    cores = NULL, toler = FALSE, SNPSE = FALSE, parallel = TRUE,
#                                   Output = NULL, GC="standard", MPI=FALSE)


# calculating sample size for factors (for PRS software)

#restrict to MAF of 40% and 10%
#PSYCH_factor2<-subset(PSYCH_factor, PSYCH_factor$MAF <= .4 & PSYCH_factor$MAF >= .1)
#calculate expected sample size (N_hat)
#N_hat<-mean(1/((2*PSYCH_factor2$MAF*(1-PSYCH_factor2$MAF))*PSYCH_factor2$SE^2))