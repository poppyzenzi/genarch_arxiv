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

# confirmatory factor analysis

#Specify the Genomic confirmatory factor model
# MDD suffers Heywoods case - need to prevent negative residuals

CFAofEFA <- 'F1 =~ ANX + NEU + MDD
             F2 =~ a*BIP + a*SCZ
             F3 =~ b*ADHD + b*ASD

MDD~~a*MDD
a > .001

F1~~F2
F2~~F3
F1~~F3'

#run the model
# std.lv removes need for NA* and 1 constraints
Psycho <- usermodel(psycho, estimation = "DWLS", model = CFAofEFA,
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Psycho$modelfit
Psycho$results