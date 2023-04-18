# script to run GWAS of user-specified model from EFAtoCFA
require(GenomicSEM)
library(GenomicSEM)
library(readr)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/")

LDSCoutput <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/LDSCoutput.rds")
PSYCH_sumstats <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/PSYCH_sumstats.rds")

covstruct <- LDSCoutput
sumstats <- PSYCH_sumstats

################## parallelising ########################

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
cat(paste("Running task", k, " of ", max_tasks, "\n"))
out_tsv <- paste0('gwas/', k, '.tsv')
dir.create(dirname(out_tsv), showWarnings=FALSE)

# get subset of SNPs to analyze in this task
nsnps <- nrow(sumstats)
snps_per_task <- ceiling(nsnps / max_tasks)

# index of task IDs repeated for number of SNPs per ask, clip the last task# to the actual number of snps
task_allocations <- rep(1:max_tasks, each=snps_per_task)[seq.int(nsnps)]
SNPs <- sumstats[which(task_allocations == k),]

rm(sumstats)
######################################################

#specify the model
# not sure if need to specify MDD negative residual variance here as in CFA model?
model<-"F1 =~ ANX + NEU + MDD
        F2 =~ BIP + SCZ
        F3 =~ ADHD + ASD
F1~~F2
F2~~F3
F1~~F3
F1 ~ SNP
F2 ~ SNP
F3 ~ SNP"

#run the multivariate GWAS using parallel processing
CorrelatedFactors<-userGWAS(covstruc = covstruct, SNPs = SNPs,
                            estimation = "DWLS", model = model, printwarn = TRUE, sub=c("F1~SNP", "F2~SNP, F3~SNP"),
                            cores = 8, toler = FALSE, SNPSE = FALSE, parallel = TRUE, Output = NULL,
                            GC="standard", MPI=FALSE, smooth_check=TRUE)

write_tsv(gwas, out_tsv)



##Calculate Effective Sample Size for next steps
#restrict to MAF of 40% and 10%
CorrelatedFactors1<-subset(CorrelatedFactors, CorrelatedFactors$MAF <= .4 & CorrelatedFactors$MAF >= .1)
# for Factor 1
N_hat_F1<-mean(1/((2*CorrelatedFactors1[[1]]$MAF*(1-CorrelatedFactors1[[1]]$MAF))*CorrelatedFactors1[[1]]$SE^2))
# for Factor 2
N_hat_F1<-mean(1/((2*CorrelatedFactors1[[2]]$MAF*(1-CorrelatedFactors1[[2]]$MAF))*CorrelatedFactors1[[2]]$SE^2))
# for Factor 3
N_hat_F1<-mean(1/((2*CorrelatedFactors1[[3]]$MAF*(1-CorrelatedFactors1[[3]]$MAF))*CorrelatedFactors1[[3]]$SE^2))