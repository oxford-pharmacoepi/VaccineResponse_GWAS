#!/bin/bash

directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Bed_files"
imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"


for i in {6..6}; do
    run_plink_wes="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first\
      --sample ukb22828_c${i}_b0_v3.sample\
      --no-pheno\
      --make-bed --out c${i}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done