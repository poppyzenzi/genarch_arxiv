#!/bin/bash
#$ -N prsice_alspac_mdd3_unrelated_euro
#$ -l h_rt=24:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 8
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

# this is a shell (bash) script to generate PRS for MDD in ALSPAC using MJA's summary scores on datastore 

# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0

# base and target dat are in HOME path - set 
# bed bim fam in here copied from AES datastore 

HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs

######## make pheno file ########

# this is in a separate R script
# for now run before from command line as problems with R packages 
# Rscript /alspac/generate_pheno_alspac.R


######################### C+T calculation ########################

# . /etc/profile.d/modules.sh
# module add igmm/apps/R/4.0.3

# -------------------------------------------------
# Create PRSs
# -------------------------------------------------

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
	--prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
    	--base /exports/igmm/eddie/GenScotDepression/users/poppy/mdd3_ss_MAF.txt \
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
    	--out $HOME/output/abcd_mdd_prs_$(date +%y%m%d)

# --extract valid SNPs as there was one duplicate identified