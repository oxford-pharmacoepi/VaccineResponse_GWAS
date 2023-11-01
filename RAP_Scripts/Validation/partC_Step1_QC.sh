#!/bin/sh

# It takes about 10mins to run

# How to Run:
# ./partC_Step1_QC.sh in the command line

# What this file does:
# Runs the quality control process using PLINK on the merged file generated in partB_Merge.sh

directory="MAH/one_dose_validation" #Output directory, example: MAH/cohorts
phenotype="one_dose_validation" #Example: one_dose_cohort

run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged\
 --keep ${phenotype}.phe\
 --autosome\
 --maf 0.01\
 --mac 20\
 --geno 0.1\
 --hwe 1e-15\
 --mind 0.1\
 --write-snplist\
 --write-samples\
 --no-id-header\
 --out  snps_qc_pass_${phenotype}"

dx run swiss-army-knife -iin="/${directory}/ukb22418_c1_22_v2_merged.bed" \
   -iin="/${directory}/ukb22418_c1_22_v2_merged.bim" \
   -iin="/${directory}/ukb22418_c1_22_v2_merged.fam"\
   -iin="/${directory}/${phenotype}.phe" \
   -icmd="${run_plink_qc}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory}/" --brief --yes --name="StepC_${phenotype}"
