# R script to run multivariable LDSC on eddie for GSEM step 2eddie

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")

#vector of munged summary statistics
traits<-c("ADHD.sumstats.gz","BIP.sumstats.gz","SCZ.sumstats.gz", "ASD.sumstats.gz", "MDD.sumstats.gz",
            "NEU.sumstats.gz", "ANX.sumstats.gz")

# enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
# if Neff available: 0.5, if continuous: NA
# only calculate for asd (case/case+control)
# neuroticism (continuous) = NA
sample.prev<-c(0.5, 0.5, 0.5, 0.396569579, 0.5, NA, 0.5)

# vector of population prevalences [from large scale epi studies, first look in ref GWAS]
# NA if continuous
population.prev<-c(0.05,0.02,0.01,0.012,0.15,NA,0.16)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/")
# LD scores must match ancestry of sumstats
#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("ADHD","BIP","SCZ", "ASD", "MDD", "NEU", "ANX")

# run LDSC
# stand = optional; output correlation and covariance matrices
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,
                    ld=ld,wld=wld,trait.names=trait.names,stand=TRUE)

# optional to save standard errors
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))

#optional command to save the output as a .RData file for later use
# make this folder on eddie
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/")
save(LDSCoutput,file="LDSCoutput.RData")
