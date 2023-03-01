igmm:#!/bin/bash
#$ -N prsice_alspac_mdd3_unrelated
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/PRS/job_logs
#$ -M s2421111@ed.ac.uk
#$ -m baes

# this is a shell (bash) script to generate PRS for MDD in ALSPAC using MJA's summary scores on datastore 


# look at AES master file to see how to run queue multiple shell scripts using qsub

# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0
 
# base and target dat are in HOME path - set 
HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/alspac

############ BGEN > PGEN > MERGE > BED/BIM/FAM #############

cd /exports/cmvm/datastore/scs/groups/ALSPAC/data/genomics/B3421/genetics/1000G_2021-01-25/all1/
cp -r data/ /exports/eddie/scratch/s2421111/

cd /exports/eddie/scratch/s2421111/data/

mv data_chr01.bgen data_chr1.bgen
mv data_chr02.bgen data_chr2.bgen
mv data_chr03.bgen data_chr3.bgen
mv data_chr04.bgen data_chr4.bgen
mv data_chr05.bgen data_chr5.bgen
mv data_chr06.bgen data_chr6.bgen
mv data_chr07.bgen data_chr7.bgen
mv data_chr08.bgen data_chr8.bgen
mv data_chr09.bgen data_chr9.bgen
# need to do this?

# Job to convert BGEN files to PGEN.CH
CHR=$SGE_TASK_ID

SCRATCH=/exports/eddie/scratch/s2421111/
cd $SCRATCH

mkdir -p alspac

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--bgen data/data_chr${CHR}.bgen 'ref-first' \
--sample data/data.sample \
--keep data/unrelated_children.txt \
--rm-dup 'exclude-all' \
--make-pgen 'vzs' \
--out ALSPAC/data_chr${CHR}_QC \
--memory 4000 \
--threads 1

###### Merge PGENs into one autosome file ###################
#############################################################

SCRATCH=/exports/eddie/scratch/$USER/
cd $SCRATCH

# Paste prefix of PGEN file names and save as txt file.
echo ${SCRATCH}ALSPAC/data_chr1_QC > $SCRATCH/mergelist.txt
for CHR in {2..22}; do
	echo ${SCRATCH}ALSPAC/data_chr${CHR}_QC >> $SCRATCH/mergelist.txt
done

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
	--pmerge-list $SCRATCH/mergelist.txt 'pfile-vzs' \
	--make-pgen \
	--out $SCRATCH/ALSPAC/data_QC \
	--threads 1


######## make bed bim fam ###########
SCRATCH=/exports/eddie/scratch/$USER

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile $SCRATCH/ALSPAC/data_QC \
--maf 0.01 \
--mac 100 \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out $SCRATCH/ALSPAC/data_QC \
--threads 4




######## make dummy pheno file ########

# this is in a separate R script

Rscript /alspac/dummy_pheno_alspac.R





################# C+T calculation ##############


cd /exports/igmm/eddie/GenScotDepression/amelia/ALSPAC_inflam_episodes/PRS
    
. /etc/profile.d/modules.sh
module add igmm/apps/R/4.0.3

# -------------------------------------------------
# Create PRSs
# -------------------------------------------------
# will need to change some params

Rscript /exports/igmm/eddie/GenScotDepression/amelia/packages/PRSice_v2.3.3/PRSice.R \
    --dir . \
    --prsice /exports/igmm/eddie/GenScotDepression/amelia/packages/PRSice_v2.3.3/PRSice \
    --base $sumstats \
    --target ALSPAC/1000G/data_QC \
    --beta \
    --A1 A1 \
    --A2 A2 \
    --pvalue p \
    --snp SNP \
    --stat b \
    --no-default \
    --binary-target T \
    --pheno-file alspac/no_pheno.txt \
    --pheno-col no_pheno \
    --fastscore \
    --bar-levels 5e-8,5e-7,5e-6,5e-5,5e-4,0.001,0.05,0.1,0.2,0.5,1 \
    --all-score  \
    --clump-r2 0.25 \
    --clump-kb 500 \
    --thread 8 \
    --print-snp \
    --out Output/$output
