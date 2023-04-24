#!/bin/bash
#$ -N parcel_gwas_script
#$ -l h_rt=4:00:00
#$ -t 1-100
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes
#$ -cwd

# load modules
. /etc/profile.d/modules.sh
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0

Rscript parcels_gwas.R