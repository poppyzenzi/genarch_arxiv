library(magrittr)
library(car)
library(dplyr)

# each pheno column generates a separate PRS output file in PRSice-2
# pheno columns containing only 0s or only 1s are rejected by PRSice-2 (dummy)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs")

# read in .fam file of QC'd target data to get list of IDs
fam <- read.table("data_QC.fam")

# make FID, IID and dummy cols.  
pheno <- data.frame("FID" = fam[,1], "IID" = fam[,1], "no_pheno" = 0)
# Add random 0 and 1s to no_pheno dummy col
pheno$no_pheno <- sample(c(0,1), replace = T, size = nrow(pheno))

###############################################

# smfq score pheno col binary=F

# read in long dat
df <- readRDS('/exports/eddie/scratch/s2421111/alspac/alspac_smfq_long.rds')
df <- df %>% subset(select = c('IID', 'id', 'time', 'age', 'sex', 'dep'))

# here will change depending on which smfq time point you want to model the phenotype against 
df <- subset(df, time == 1)
pheno <- merge(pheno, df, by = c("IID"))

#############################################

# save table in your scratch dir
write.table(pheno, "pheno_file_alspac.txt", row.names = F, quote = F)
