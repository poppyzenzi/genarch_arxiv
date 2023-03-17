#!/bin/bash
#$ -N prsice_abcd_mdd3_unrelated
#$ -l h_rt=1:00:00
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# this is a shell (bash) script to generate PRS for MDD in ABCD using MJA's summary scores on datastore


# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0

# base and target dat are in HOME path - set

HOME=/exports/igmm/eddie/GenScotDepression/users/poppy
SCRATCH=/exports/eddie/scratch/s2421111/gwas

# Create dummy pheno file
# This has just IID which should match target dat and
# Rscript dummy_pheno_test.R


# anxiety SNP pos chr A1 A2 af info BETA SE P
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $SCRATCH/TotAnx_effect_sumstats.gz \
        --target $HOME/abcd/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --binary-target F \
        --pheno-file $HOME/abcd/pheno_file_bpm.txt \
        --snp SNP --chr chr --bp pos --A1 A1 --A2 A2 --stat BETA --pvalue P \
        --beta \
        --pheno-col int \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/abcd/output/abcd_bpm_anx_prs_$(date +%d%m)

# adhd CHR SNP BP A1 A2 FRQ_A_19099 FRQ_U_34194 INFO OR SE P
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $SCRATCH/daner_adhd_meta_filtered.meta.gz \
        --target $HOME/abcd/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --binary-target F \
        --pheno-file $HOME/abcd/pheno_file_bpm.txt \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --pvalue P \
        --or \
        --pheno-col int \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/abcd/output/abcd_bpm_adhd_prs_$(date +%d%m)

# autism CHR SNP BP A1 A2 INFO OR SE P
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $SCRATCH/iPSYCH-PGC_ASD_Nov2017.gz \
        --target $HOME/abcd/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --binary-target F \
        --pheno-file $HOME/abcd/pheno_file_bpm.txt \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --pvalue P \
        --or \
        --pheno-col int \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/abcd/output/abcd_bpm_asd_prs_$(date +%d%m)

# schizophrenia CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE PVAL NCAS NCON NEFF
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $SCRATCH/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz \
        --target $HOME/abcd/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --binary-target F \
        --pheno-file $HOME/abcd/pheno_file_bpm.txt \
        --snp ID --chr CHROM --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --beta \
        --pheno-col int \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/abcd/output/abcd_bpm_scz_prs_$(date +%d%m)


# bipolar

# neuroticism


# --all-score meaning you generate polygenic risk scores for all thresholds
# --print-snp prints a list of snps which were used in the PRS to a .snp file
# --keep this is a list of IID of all the unrelated individuals in the GenScot
# --snp, car, bp, A1, A2, stat, pvalue: columns in the base file
# --beta/OR shows whether it's beta or OR
# --maf filters SNPs based on MAF - check base data QC levels
# --bar-levels PRSice generated bar chart levels
# --binary-target T or F
# --pheno-file give dummy pheno generated in R script
# --out output files with date appended