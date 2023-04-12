require(GenomicSEM)
library(GenomicSEM)
library(readr)
# munged sum stats in "/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/"


#run multivariable LDSC to create the S and V matrices
psycho <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/LDSCoutput.rds")

#smooth the S matrix for EFA using the nearPD function in the Matrix package.
require(Matrix)
Ssmooth<-as.matrix((nearPD(psycho$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")

#print the loadings
EFA$loadings
