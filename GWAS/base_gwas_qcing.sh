#!/bin/bash
#$ -N prsice_gwas_qcing
#$ -l h_rt=24:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 8
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes


# base data are in scratch space - set
# save in eddie space for analysis [if too big, save in scratch but remember this gets cleaned monthly]

HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS
SCRATCH=/exports/eddie/scratch/s2421111/gwas/

######################### QC STEPS########################



########## filter for MAF > 0.01 (cases+controls, or all if score), and INFO > 0.8 ########

# PCG3_SCZ_european (cc)
gunzip -c $SCRATCH/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz |\
awk 'NR==1 || ($6 > 0.01) && ($7 > 0.01) && ($8 > 0.8) {print}' |\.
gzip  > PGC3_SCZ_wave3_euro_MAF0.01_INFO0.8.gz


# PCG_bip2021 (cc)
gunzip -c $SCRATCH/pgc-bip2021-all.vcf.tsv.gz |\
awk 'NR==1 || ($10 > 0.01) && ($11 > 0.01) && ($12 > 0.8) {print}' |\.
gzip  > pgc_bip2021_all_MAF0.01_INFO0.8.gz


# TotAnx_effect (from Alex) (fs)
gunzip -c $SCRATCH/TotAnx_effect_sumstats.gz |\
awk 'NR==1 || ($6 > 0.01) && ($7 > 0.8) {print}' |\.
gzip  > TotAnx_effect_sumstats_MAF0.01_INFO0.8.gz

# daner_adhd_meta (cc)
gunzip -c $SCRATCH/daner_adhd_meta_filtered.meta.gz |\
awk 'NR==1 || ($6 > 0.01) && ($7 > 0.01) && ($8 > 0.8) {print}' |\.
gzip  > daner_adhd_meta_MAF0.01_INFO0.8.gz

# pgc_asd_17. Only OR col - already been QC’d for MAF?
gunzip -c $SCRATCH/iPSYCH-PGC_ASD_Nov2017.gz |\
awk 'NR==1 || ($6 > 0.8) {print}' |\.
gzip  > iPSYCH-PGC_ASD_Nov2017_INFO0.8.gz



###################### remove duplicates (check documentation) ####################

# PCG3_SCZ_european
gunzip -c $SCRATCH/PGC3_SCZ_wave3_euro_MAF0.01_INFO0.8.gz |\
awk '{seen[$2]++; if(seen[$2]==1){ print}}' |\.
gzip > PGC3_SCZ_wave3_euro_MAF0.01_INFO0.8_nodup.gz


# PCG_bip2021
gunzip -c $SCRATCH/pgc_bip2021_all_MAF0.01_INFO0.8.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\.
gzip > pgc_bip2021_all_MAF0.01_INFO0.8_nodup.gz


# TotAnx_effect
gunzip -c $SCRATCH/TotAnx_effect_sumstats_MAF0.01_INFO0.8.gz |\
awk '{seen[$1]++; if(seen[$1]==1){ print}}' |\.
gzip > TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup.gz


# daner_adhd_meta (cc)
gunzip -c $SCRATCH/daner_adhd_meta_MAF0.01_INFO0.8.gz |\
awk '{seen[$2]++; if(seen[$2]==1){ print}}' |\.
gzip > daner_adhd_meta_MAF0.01_INFO0.8_nodup.gz


# pgc_asd_17
gunzip -c $SCRATCH/iPSYCH-PGC_ASD_Nov2017_INFO0.8.gz |\
awk '{seen[$2]++; if(seen[$2]==1){ print}}' |\.
gzip > iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup.gz


############################# remove ambiguous SNPs ###########################

# PCG3_SCZ_european
gunzip -c $SCRATCH/PGC3_SCZ_wave3_euro_MAF0.01_INFO0.8_nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > PGC3_SCZ_wave3_euro_MAF0.01_INFO0.8_nodup_noambig.gz


# PCG_bip2021
gunzip -c $SCRATCH/pgc_bip2021_all_MAF0.01_INFO0.8_nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > pgc_bip2021_all_MAF0.01_INFO0.8_nodup_noambig.gz


# TotAnx_effect
gunzip -c $SCRATCH/TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > TotAnx_effect_sumstats_MAF0.01_INFO0.8_nodup_noambig.gz


# daner_adhd_meta (cc)
gunzip -c $SCRATCH/daner_adhd_meta_MAF0.01_INFO0.8_nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > daner_adhd_meta_MAF0.01_INFO0.8_nodup_noambig.gz


# pgc_asd_17
gunzip -c $SCRATCH/iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz


######################### mismatched SNPs removed by PRSice directly ########################

######################### what about unrelated? ####################################

######## don’t forget to check if A1 and A2 are effect or non-effect alleles before running prsice

