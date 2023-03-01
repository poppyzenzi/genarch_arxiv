import os
import pandas as pd
import numpy as np
import pyreadr  # may be issues with this on eddie > pip install pyreadr from command line

# Set working directory
os.chdir("/exports/igmm/eddie/GenScotDepression/users/poppy")

# Read in .fam file of QC'd target data to get list of IDs
fam = pyreadr.readr("ABCD3.0_imputed_whiteonly_MAF0.01_unrelated.fam", header=None, sep=" ")

# Extract FID and IID as the NDAR IDs from the unrelated sample
# FID contains no multiple values for an unrelated sample. So if unrelated, FID and IID will be different
pheno = pd.DataFrame({"IID": fam.iloc[:, 1], "no_pheno": 0})

# Add random 0 and 1s to no_pheno
pheno["no_pheno"] = np.random.choice([0, 1], replace=True, size=len(pheno))

# Read in mental health data
df = pd.pyre("/Volumes/igmm/GenScotDepression/data/abcd/release4.0/iii.data/Mental_Health/abcd_cbcls01.pkl")
df = df[["src_subject_id", "eventname", "cbcl_scr_dsm5_depress_r", "interview_age", "sex"]]
df.columns = ["IID", "time", "dep", "age", "sex"]

df["time"] = df["time"].replace({
    "baseline_year_1_arm_1": 0,
    "1_year_follow_up_y_arm_1": 1,
    "2_year_follow_up_y_arm_1": 2,
    "3_year_follow_up_y_arm_1": 3
})
df["age"] = pd.to_numeric(df["age"]) / 12
df["dep"] = pd.to_numeric(df["dep"])
df["sex"] = df["sex"].replace({"M": 0, "F": 1})

# choose cbcl time point for phenotype
df = df[df["time"] == 3]

# Merge pheno and mental health data on IID
pheno = pd.merge(pheno, df, on="IID")

# Save table in your scratch dir - use this with --pheno-file and --pheno-col command in PRSice
pheno.to_csv("pheno_file.txt", sep="\t", index=False)
