require(GenomicSEM)
library(GenomicSEM)
library(readr)

# script to compare multivariate model fits
# first perform EFA to find parcel structure (hypothesise 3 parcels)
# compare fit between common, correlated and higher-order factor


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

model1 <- 'F1 =~ ANX + NEU + MDD
             F2 =~ a*BIP + a*SCZ
             F3 =~ b*ADHD + b*ASD

MDD~~a*MDD
a > .001

F1~~F2
F2~~F3
F1~~F3'

#run the model
# std.lv removes need for NA* and 1 constraints
correlatedfit <- usermodel(psycho, estimation = "DWLS", model = model1,
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

correlatedfit$modelfit
correlatedfit$results

########### compare with common factor model ###########

model2 <- 'F1 =~ ANX + NEU + MDD + BIP + SCZ + ADHD + ASD
MDD ~~ c*MDD
c > 0.001'

common_fit <- usermodel(psycho, estimation = "DWLS", model = model2,
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

common_fit$modelfit
common_fit$results


########### compare with hierarchical model ###########

model3 <- "F1 =~ ANX + NEU + MDD
          F2 =~ BIP + SCZ
          F3 =~ ADHD + ASD

P =~ F1 + F2 + F3

MDD ~~ c*MDD
c > 0.001"

higherfit <- usermodel(psycho, estimation = "DWLS", model = model3,
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

higherfit$modelfit
higherfit$results
