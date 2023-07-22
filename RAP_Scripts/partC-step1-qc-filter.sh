#!/bin/sh
# This script runs the QC process using PLINK on the merged file generated in 
# partB-merge-files.sh as described in the 

directory="" #Output directory, example: MAH/cohorts
phenotype="" #Example: one_dose_cohort

run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged\
 --keep ${phenotype}.phe --autosome\
 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15\
 --mind 0.1 --write-snplist --write-samples\
 --no-id-header --out  snps_qc_pass_one_${phenotype}"

dx run swiss-army-knife -iin="/${directory}/ukb22418_c1_22_v2_merged.bed" \
   -iin="/${directory}/ukb22418_c1_22_v2_merged.bim" \
   -iin="/${directory}/ukb22418_c1_22_v2_merged.fam"\
   -iin="/${directory}/${phenotype}.phe" \
   -icmd="${run_plink_qc}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory}/" --brief --yes
