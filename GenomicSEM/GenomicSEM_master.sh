#!/bin/bash
#$ -N gsem_master_7ss
#$ -l h_rt=24:00:00
#$ -l h_vmem=64G
#$ -pe sharedmem 2
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# load modules
. /etc/profile.d/modules.sh
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0


# GenomicSEM without individual SNP effects
Rscript munging.R
Rscript multi_LDSC.R
Rscript commonfactor_pp.R

# with SNP effects
Rscript EFA_model.R
Rscript CFA_model.R

# common factor GWAS
Rscript PSYCH_ss_prep.R
Rscript commonfactor_gwas.R


############## GSEM P-FACTOR PGS GENERATION #################

# make PRS with GSEM indivudal SNP effects as the GWAS
module load igmm/apps/PRSice/2.1.11

HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/
BASE=/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/
ALSPAC=/exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs
ABCD=/exports/igmm/eddie/GenScotDepression/users/poppy/abcd

# alspac
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/LDSCoutput \
        --target $ALSPAC/data_QC \
        --beta \
        --snp ID --chr CHR --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --pheno-file $ALSPAC/pheno_file_alspac.txt \
        --pheno-col dep \
        --binary-target F \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/alspac_gsem_prs$(date +%m%d)


# abcd
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/LDSCoutput\
        --target $ABCD/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --beta \
        --snp ID --chr CHR --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --pheno-file $HOME/pheno_file_bpm.txt.txt \
        --pheno-col int \
        --binary-target F \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/abcd_gsem_prs$(date +%m%d)