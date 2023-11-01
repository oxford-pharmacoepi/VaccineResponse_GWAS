#!/bin/bash


#change exome_file_dir and data_field for the newest release
imputed_file_dir="/Imputation/"
directory="" #Output directory, example: MAH/cohorts
phenotype="" #Example: one_dose_cohort
outcome="" # Column where the outcome is in the phenotype file

for chr in {1..22}; do
  run_regenie_cmd="regenie --step 2 --bed ukb22828_c${chr}_b0_v3 --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_snps_qc_pass_${phenotype}.snplist\
    --phenoCol ${outcome}\
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
    --pred ${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bed"\
   -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bim"\
   -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.fam"\
   -iin="/${directory}/c${chr}_snps_qc_pass_${phenotype}.snplist"\
   -iin="/${directory}/${phenotype}.phe"\
   -iin="/${directory}/${phenotype}_results_pred.list"\
   -iin="/${directory}/${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_hdd1_v2_x16"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory}/" --brief --yes
done
