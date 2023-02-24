# Create binary pheno file where 1 indicates precence of disorder and 0 none
# First two columns should be FID (family id) and IID (sampel id)
# each pheno column generates a separate PRS output file in PRSice-2
# third col a column of 0s and 1s as a no_pheno column.
# pheno columns containing only 0s or only 1s are rejected by PRSice-2 making this a dummy pheno file
# column names do not matter

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy")

# read in .fam file of QC'd target data to get list of IDs
fam <- read.table("ABCD3.0_imputed_whiteonly_MAF0.01_unrelated.fam")

#extracts FID and IID as the NDAR IDs from the unrelated sample
# FID contains no multiple values for an unrelated sample. So if unrelated, FID and IID will be different
pheno <- data.frame("FID" = fam[,1], "IID" = fam[,1], "no_pheno" = 0)

# Add random 0 and 1s to no_pheno
pheno$no_pheno <- sample(c(0,1), replace = T, size = nrow(pheno))

# save table in your scratch dir - use this with --pheno-file command in PRSice 
write.table(pheno, "no_pheno.txt", row.names = F, quote = F)