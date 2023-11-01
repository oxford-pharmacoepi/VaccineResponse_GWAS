#!/bin/sh

# The time it takes to run depends on the chromosome you are extracting the snps from

# How to Run:
# ./partD_Step1_Regenie_SelectSNPs.sh in the command line

# What this file does:
# Extracts the SNPs listed in validation.txt from the Imputation files

directory="MAH/one_dose_validation" #Output directory, example: MAH/cohorts
phenotype="one_dose_validation" #Example: one_dose_cohort

data="/Bulk/Imputation/UKB\ imputation\ from\ genotype/"

for chr in 6; do
run_regenie_step1="plink2 --bgen ukb22828_c${chr}_b0_v3.bgen ref-first\
        --sample ukb22828_c${chr}_b0_v3.sample\
	-extract validation_${chr}.txt\
	--make-bed\
	--out selected_snps_${phenotype}_${chr}\
	--no-pheno\
	--keep ${phenotype}.phe"

dx run swiss-army-knife -iin="${data}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${data}/ukb22828_c${chr}_b0_v3.sample"\
   -iin="/${directory}/validation_${chr}.txt"\
   -iin="/${directory}/${phenotype}.phe"\
   -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --name="StepD.1_${phenotype}_chr${chr}"\
   --destination="/${directory}/" --brief --yes
 
   
done
