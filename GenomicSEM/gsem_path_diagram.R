# path diagram

# "F1 =~ ?" lists the indicator loadings on the common factors. "V ~~ V" lists the residual variances of the indicators
# after removing variance explained by the common factor. In the standardized case of a common factor model
# ("V ~~ V" + "F1 =~ V"^2) will sum to 1.

library(lavaan)
library(lavaanPlot)

setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/GenomicSEM/")

data <- ("common_factor_result.txt")

HS.model <- 'F  =~  ADHD + BIP + SCZ + ASD + MDD + NEU + ANX'

fit <- cfa(HS.model, data=data)

lavaanPlot(model = fit, edge_options = list(color = "grey"))

