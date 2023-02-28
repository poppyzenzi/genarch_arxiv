#!/bin/bash
#$ -N prsice_abcd_mdd3_unrelated
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -M s2421111@ed.ac.uk
#$ -m baes

# this is a shell (bash) script to generate PRS for MDD in ABCD using MJA's summary scores on datastore 

# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0
 
# base and target dat are in HOME path - set 
HOME=/exports/igmm/eddie/GenScotDepression/users/poppy


# Create dummy pheno file, output is 'abcd/no_pheno.txt'
# This has just IID which should match target dat and 
Rscript dummy_pheno_abcd.R

# which <file> to find where PRSice.R is located - easier to go direct path

Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
	--prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
	--base $HOME/mdd3_ss_MAF.txt \
	--target $HOME/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
	--binary-target T \
    	--maf-base MAF,0.01 \
    	--info-base IMPINFO,0.8 \
	--phleno-file $HOME/abcd/no_pheno.txt \
	--snp ID --chr CHR --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
    	--beta \
	--pheno-col no_pheno \
	--ignore-fid \
	--print-snp \
	--fastscore \
	--bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
  --out $HOME/test_run/abcd_mdd_prs_test_$(date +%y%m%d)


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

