# step 4 - final step, run the GWAS
# mood, psych, neurodev
library(GenomicSEM)
library(readr)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/parcel_factors/gwas")

# run the multivariate GWAS using parallel processing
# cov struct [LDSC] and sum stats
# change for each parcel or loop over in SHELL
covstruct <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/parcel_factors/LDSCoutput_mood.rds")
sumstats <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/parcel_factors/MOOD_sumstats.rds")

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
cat(paste("Running task", k, " of ", max_tasks, "\n"))
out_tsv <- paste0('gwas_output/', k, '.tsv')   # make direc in gsem/parcel_factors/gwas_output
dir.create(dirname(out_tsv), showWarnings=FALSE)

# get subset of SNPs to analyze in this
tasknsnps <- nrow(sumstats)
snps_per_task <- ceiling(nsnps / max_tasks)

# index of task IDs repeated for number of SNPs per ask, clip the last task# to the actual number of snps
task_allocations <- rep(1:max_tasks, each=snps_per_task)[seq.int(nsnps)]
SNPs <- sumstats[which(task_allocations == k),] rm(sumstats)
gwas <- commonfactorGWAS(covstruc=covstruct, SNPs=SNPs, smooth_check=TRUE, parallel=TRUE, cores=8, toler=1e-50)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/gwas/output")

save(gwas, file="MOOD_GWAS.txt", out_tsv)

