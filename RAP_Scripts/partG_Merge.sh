#!/bin/sh

merge_cmd='out_file="oneDose_assoc.regenie.merged.txt"

# Use dxFUSE to copy the regenie files into the container storage
cp /mnt/project/Whole_genome/Breakthrough_gwas/Intermediary_files/*.regenie.gz .
gunzip *.regenie.gz

# add the header back to the top of the merged file
echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tINFO\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file

files="./*.regenie"
for f in $files
do
# for each .regenie file
# remove header with tail
# transform to tab delimited with tr
# save it into $out_file
   tail -n+2 $f | tr " " "\t" >> $out_file
done

# remove regenie files
rm *.regenie'

directory="Whole_genome/Breakthrough_gwas"
phenotype="oneDose"
outcome="immuneResponse"

dx run swiss-army-knife -iin="/${directory}/Intermediary_files/${phenotype}_assoc.c1_${outcome}.regenie.gz" \
   -icmd="${merge_cmd}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory}/" --brief --yes\
   --name="StepG_${phenotype}"
