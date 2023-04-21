```````````````````````````````````````library(magrittr)
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

# Add random 0 and 1s to no_pheno for dummy
pheno$no_pheno <- sample(c(0,1), replace = T, size = nrow(pheno))

# read in pheno smfq data from scratch space 9data copied from datastore)
df <- readRDS('/exports/eddie/scratch/s2421111/abcd/abcd_yssbpm01.rds')

df <- df %>% subset(select = c('src_subject_id', 'eventname','bpm_y_scr_internal_r', 'interview_age', 'sex')) %>%
  set_colnames(c("IID","time","int", "age", "sex"))

time_recode <- dplyr::recode(df$time,
                             `baseline_year_1_arm_1`="0",
                             `6_month_follow_up_arm_1` = "0.5",
                             `1_year_follow_up_y_arm_1`="1",
                             `18_month_follow_up_arm_1` = "1.5",
                             `2_year_follow_up_y_arm_1`= "2",
                             `30_month_follow_up_arm_1` = "2.5",
                             `3_year_follow_up_y_arm_1`= "3")


df$time <- as.numeric(time_recode)
df$age <- as.numeric(df$age)/12 #convert age to years
df$int <- as.numeric(df$int)

sex_recode <- dplyr::recode(df$sex, #recode sex to binary
                            `M`="0",
                            `F`="1")

df$sex <- as.numeric(sex_recode)
df <- subset(df, time == 0.5)


pheno <- merge(pheno, df, by = c("IID"))

pheno
# save table in your scratch dir - use this with --pheno-file command in PRSice
write.table(pheno, "pheno_file_bpm.txt", row.names = F, quote = F)```````````````````````````````````````