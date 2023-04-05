# step 4 : combine sumstats and ldsc output and run common factor GWAS
# this gets the common factor loading of each SNP

require(GenomicSEM)
library(GenomicSEM)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/")

# read in sumstats from step 3 and ldsc from step 2
PSYCH_sumstats <- "PSYCH_sumstats"
LDSCoutput <- "/ldsc/LDSCoutput.RData"


# run the multivariate GWAS using parallel processing
PSYCH_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats, estimation = "DWLS",
                                    cores = NULL, toler = FALSE, SNPSE = FALSE, parallel = TRUE,
                                    Output = NULL, GC="standard", MPI=FALSE)

# note that the code written above specifies all arguments for completeness, but as many of these arguments
# are set to the package default it could also be written as below and produce the same result:
# PSYCH_factor  <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats)


# calculating sample size for factors (for PRS software)

#restrict to MAF of 40% and 10%
#PSYCH_factor2<-subset(PSYCH_factor, PSYCH_factor$MAF <= .4 & PSYCH_factor$MAF >= .1)
#calculate expected sample size (N_hat)
#N_hat<-mean(1/((2*PSYCH_factor2$MAF*(1-PSYCH_factor2$MAF))*PSYCH_factor2$SE^2))