#!/bin/bash

imputed_file_dir="/Auxiliary_genetic_dataset/Imputation"
directory="MAH/one_dose_validation" #Output directory, example: MAH/cohorts
phenotype="one_dose_validation" #Example: one_dose_cohort

for chr in 6; do
        run_plink_wes="plink2 --bfile selected_snps_${phenotype}_${chr}\
            --no-pheno\
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
            --rm-dup force-first\
            --out c${chr}_snps_qc_pass_${phenotype}"

        dx run swiss-army-knife -iin="/${directory}/selected_snps_${phenotype}_${chr}.bed" \
            -iin="/${directory}/selected_snps_${phenotype}_${chr}.bim" \
            -iin="/${directory}/selected_snps_${phenotype}_${chr}.fam"\
            -iin="/${directory}/${phenotype}.phe" \
            -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
            --name "StepE_${phenotype}_chr${chr}"\
            --destination="${directory}/" --brief --yes
done




