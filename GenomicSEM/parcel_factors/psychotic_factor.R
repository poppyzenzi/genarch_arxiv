# =========================== step 1 munge ===============================
require(GenomicSEM)
library(GenomicSEM)
library(MASS) # for matrix

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")

files<-c("pgc_bip_qcd.txt","pgc_scz_qcd.txt")

hm3<-"/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/eur_w_ld_chr/w_hm3.snplist"
trait.names<-c("BIP", "SCZ")
N=c(NA,NA)
info.filter=0.8
maf.filter=0.01
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

# =========================== step 2 LDSC ===============================
traits<-c("ADHD.sumstats.gz","ASD.sumstats.gz")

sample.prev<-c(0.5, 0.5)
population.prev<-c(0.02,0.01)
ld<-"eur_w_ld_chr/"
wld<-"eur_w_ld_chr/"
trait.names<-c("BIP", "SCZ")
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev, ld=ld,wld=wld,trait.names=trait.names,stand=TRUE)

# optional to save standard errors
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/parcel_factors")
save(LDSCoutput,file="LDSCoutput_psychotic.RData")

# ====================== step 3 common factor model =====================
# optional step - not really necessary as done in whole factor model
#To run using DWLS estimation#
CommonFactor_DWLS <- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS

LDSCoutput$modelfit
LDSCoutput$results

# use these output vals to put into path diagram [save]
result <- CommonFactor_DWLS$results
write.table(result,'common_factor_result_psychotic.txt',sep = "\t")

# saving matrix for plotting genetic correlation heatmap
x <- LDSCoutput$S_Stand
write.matrix(x,'gen_cor_matrix_psychotic.txt',sep = "\t")