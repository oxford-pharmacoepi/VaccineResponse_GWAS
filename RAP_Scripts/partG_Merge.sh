#!/bin/sh

merge_cmd='out_file="twoDose_assoc.regenie.merged.txt"

# Use dxFUSE to copy the regenie files into the container storage
cp /mnt/project/Whole_genome/Breakthrough_gwas/Intermediary_files/*immune_response.regenie.gz .
gunzip *immune_response.regenie.gz

# add the header back to the top of the merged file
echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tINFO\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file

files="./*immune_response.regenie"
for f in $files
do
   tail -n+2 $f | tr " " "\t" >> $out_file
done

# remove regenie files
rm *immune_response.regenie'

directory="Whole_genome/Breakthrough_gwas"
phenotype="twoDose"
outcome="immune_response"

dx run swiss-army-knife -iin="/${directory}/Intermediary_files/${phenotype}_assoc.c1_${outcome}.regenie.gz" \
   -icmd="${merge_cmd}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory}/" --brief --yes\
   --name="StepG_${phenotype}"
