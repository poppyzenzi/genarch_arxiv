#!/bin/bash
#$ -N prsice_alspac_mdd_anx_neu
#$ -l h_rt=2:00:00
#$ -l h_vmem=64G
#$ -pe sharedmem 2
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# this is a shell (bash) script to generate PRS for MDD, ANX and NEU in ALSPAC

# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module unload igmm/apps/R/3.3.2
module load igmm/apps/R/4.1.0


# base and target dat are in HOME path - set
# bed bim fam in here copied from AES datastore

HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs
BASE=/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS
SCRATCH=/exports/eddie/scratch/s2421111/gwas/


######## make dummy pheno file ########

# this is in a separate R script
# for now run before from command line as problems with R packages
# Rscript /alspac/generate_pheno_alspac.R

######################### C+T calculation ########################


# -------------------------------------------------
# Create PRSs
# -------------------------------------------------

# mdd3

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/mdd3_ss_MAF.txt \
        --extract /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/output/abcd_mdd_prs_230302.valid \
        --target $HOME/data_QC \
        --beta \
        --snp ID --chr CHR --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --pheno-col dep \
        --no-default \
        --binary-target F \
        --pheno-col dep \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --clump-r2 0.25 \
        --clump-kb 500 \
        --thread 8 \
        --print-snp \
        --out $HOME/output/alspac_mdd_prs_$(date +%m%d)


# anxiety a2 = effect

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz \
        --target $HOME/data_QC \
        --binary-target F \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --snp SNP --chr chr --bp pos --A1 A2 --A2 A1 --stat BETA --pvalue P \
        --beta \
        --pheno-col dep \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/output/alspac_anx_prs_$(date +%m%d)



# neuroticism

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $SCRATCH/Luciano_2017/SummaryStats.txt \
        --target $HOME/data_QC \
        --binary-target F \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --snp rsid --chr CHR --bp BP --A1 a_1 --A2 a_0 --stat N_res_beta --pvalue p_value \
        --beta \
        --pheno-col dep \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/output/alspac_neu_prs_$(date +%m%d)





