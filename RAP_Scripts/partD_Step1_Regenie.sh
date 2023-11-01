#!/bin/sh



#output directory - this should also be where the files in 02-step1-qc-filter.sh end up
directory="MAH/one_dose" #Output directory, example: MAH/cohorts
phenotype="one_dose" #Example: one_dose_cohort
outcome="immuneResponse"


run_regenie_step1="regenie --step 1\
 --lowmem --out ${phenotype}_results --bed ukb22418_c1_22_v2_merged\
 --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
 --extract snps_qc_pass_${phenotype}.snplist --phenoCol ${outcome}\
 --covarCol Sex\
 --covarCol Age\
 --covarCol Genetic_batch\
 --covarCol PC1\
 --covarCol PC2\
 --covarCol PC3\
 --covarCol PC4\
 --covarCol PC5\
 --covarCol PC6\
 --covarCol PC7\
 --covarCol PC8\
 --covarCol PC9\
 --covarCol PC10\
 --bsize 1000 --bt --loocv --gz --threads 16"

dx run swiss-army-knife -iin="/${directory}/ukb22418_c1_22_v2_merged.bed" \
   -iin="/${directory}/ukb22418_c1_22_v2_merged.bim" \
   -iin="/${directory}/ukb22418_c1_22_v2_merged.fam"\
   -iin="/${directory}/snps_qc_pass_${phenotype}.snplist"\
   -iin="/${directory}/${phenotype}.phe" \
   --name="StepD_regenie"\
   -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/${directory}/" --brief --yes
