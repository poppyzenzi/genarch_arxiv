library(magrittr)
library(car)
library(dplyr)

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
pheno <- data.frame("IID" = fam[,2], "no_pheno" = 0)

# Add random 0 and 1s to no_pheno
pheno$no_pheno <- sample(c(0,1), replace = T, size = nrow(pheno))


df <- readRDS('/exports/eddie/scratch/s2421111/abcd/abcd_cbcls01.rds')

df <- df %>% subset(select = c('src_subject_id', 'eventname',
                               'cbcl_scr_dsm5_depress_r', 'interview_age', 'sex')) %>%
  set_colnames(c("IID","time","dep", "age", "sex"))

time_recode <- dplyr::recode(df$time,
                             `baseline_year_1_arm_1`="0",
                             `1_year_follow_up_y_arm_1`="1",
                             `2_year_follow_up_y_arm_1`="2",
                             `3_year_follow_up_y_arm_1`="3")


df$time <- as.numeric(time_recode)
df$age <- as.numeric(df$age)/12 #convert age to years
df$dep <- as.numeric(df$dep)

sex_recode <- dplyr::recode(df$sex, #recode sex to binary
                            `M`="0",
                            `F`="1"
)

df$sex <- as.numeric(sex_recode)
df <- subset(df, time == 1)


pheno <- merge(pheno, df, by = c("IID"))

pheno
# save table in your scratch dir - use this with --pheno-file command in PRSice
write.table(pheno, "pheno_file.txt", row.names = F, quote = F)