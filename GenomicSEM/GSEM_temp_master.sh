#!/bin/bash
#$ -N gsem_master_7ss
#$ -l h_rt=24:00:00
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes
#$ -cwd

# load modules
. /etc/profile.d/modules.sh
module unload igmm/apps/R/3.5.1
module load roslin/R/4.1.0

# GenomicSEM without individual SNP effects
# Rscript munging.R munging is done
Rscript multi_LDSC_CFM.R

# common factor GWAS [outputs sum stats for gsem PRS below]
Rscript PSYCH_sumstats_calc.R
Rscript PSYCH_factor_calc.R

############## GSEM P-FACTOR PGS GENERATION #################

# make PRS with GSEM indivudal SNP effects as the GWAS
module load igmm/apps/PRSice/2.1.11

# base data and saving output
HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/
BASE=/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/PSYCH_sumstats
# target data
ALSPAC=/exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs
ABCD=/exports/igmm/eddie/GenScotDepression/users/poppy/abcd

# alspac, check time point
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/ \
        --target $ALSPAC/data_QC \
        --beta \
        --snp  --chr  --bp  --A1  --A2  --stat  --pvalue  \
        --pheno-file $ALSPAC/pheno_file_alspac.txt \
        --pheno-col dep \
        --binary-target F \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/alspac_gsem_prs$(date +%m%d)


# abcd , check time point
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/ \
        --target $ABCD/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --beta \
        --snp  --chr  --bp  --A1  --A2  --stat  --pvalue  \
        --pheno-file $HOME/pheno_file_bpm.txt.txt \
        --pheno-col int \
        --binary-target F \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/abcd_gsem_prs$(date +%m%d)