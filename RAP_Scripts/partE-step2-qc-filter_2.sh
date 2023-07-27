#!/bin/bash


imputed_file_dir="/Imputation/"
directory="" #Output directory, example: MAH/cohorts
phenotype="" #Example: one_dose_cohort

for i in {3..17}; do
    run_plink_wes="plink2 --bfile ukb22828_c${i}_b0_v3\
      --no-pheno --keep ${phenotype}.phe --autosome\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_snps_qc_pass_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bed" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bim" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.fam"\
     -iin="/${directory}/${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory}/" --brief --yes
done

