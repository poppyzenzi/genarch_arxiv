#!/bin/bash
#$ -N 3_parcel_factors
#$ -l h_rt=12:00:00
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

# GenomicSEM steps 1-3 to prep munged sum stats and LDSC for parcel GWAS
Rscript mood_factor.R
Rscript neurodev_factor.R
Rscript psychotic_factor.R

# next step to run GWAS for each