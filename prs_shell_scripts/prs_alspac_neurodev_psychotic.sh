#!/bin/bash
#$ -N prsice_alspac_scz_bip_adhd_asd
#$ -l h_rt=1:00:00
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# this is a shell (bash) script to generate PRS for SCZ, BIP, ASD and ADHD in ALSPAC

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

# scz

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/pgc_scz_qcd.txt \
        --extract /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/output/abcd_mdd_prs_230302.valid \
        --target $HOME/data_QC \
        --beta \
        --snp ID --chr CHROM --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --binary-target F \
        --pheno-col dep \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/alspac_scz_prs_$(date +%m%d)


# bipolar (all)

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/pgc_bip_qcd.txt \
        --extract /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/output/abcd_mdd_prs_230302.valid \
        --target $HOME/data_QC \
        --beta \
        --snp ID --chr CHROM --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --pheno-col dep \
        --binary-target F \
        --pheno-col dep \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/alspac_bip_prs_$(date +%m%d)


# adhd

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz \
        --extract /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/output/abcd_mdd_prs_230302.valid \
        --target $HOME/data_QC \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --pvalue P \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --pheno-col dep \
        --binary-target F \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/alspac_adhd_prs_$(date +%m%d)


# asd

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $BASE/iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz \
        --extract /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/output/abcd_mdd_prs_230302.valid \
        --target $HOME/data_QC \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --pvalue P \
        --pheno-file $HOME/pheno_file_alspac.txt \
        --pheno-col dep \
        --binary-target F \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/output/alspac_asd_prs_$(date +%m%d)





