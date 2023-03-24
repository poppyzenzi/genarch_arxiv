#!/bin/bash
#$ -N munging_7_sumstats
#$ -l h_rt=1:00:00
#$ -l h_vmem=64G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# load modules
. /etc/profile.d/modules.sh
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0


Rscript munging.R