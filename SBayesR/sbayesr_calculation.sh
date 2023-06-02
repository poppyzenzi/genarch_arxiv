#!/bin/bash
#$ -N sbayesr_calc
#$ -l h_rt=24:00:00
#$ -l h_vmem=64GB
#$ -pe sharedmem 4
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# Conduct SBayesR on GWAS sum stats
cd /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/sbayesr

SBR=/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/sbayesr
SCRATCH=/exports/eddie/scratch/s2421111
OUT=/exports/igmm/eddie/GenScotDepression/users/poppy
LD=/exports/igmm/datastore/GenScotDepression/data/resources/SBayesR_matrices/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.new.mldmlist

# Exclude MHC:


# access SBayesR matrices Eddie: 
# /exports/igmm/datastore/GenScotDepression/data/resources/SBayesR_matrices/ukb_50k_bigset_2.8M 
/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/sbayesr/gctb_2.04.3_Linux/gctb --sbayes R \
--mldm $LD \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--ambiguous-snp \
--impute-n \
--gwas-summary $SBR/mdd3_gwas_sbayesr.ma \
--chain-length 10000 \
--exclude-mhc \
--burn-in 2000 \
--out-freq 10 \
--out $SBR/output/mdd3_sbayesr_ss_$(date +%m%d)
