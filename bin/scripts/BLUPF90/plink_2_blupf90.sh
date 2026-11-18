#!/bin/bash

# All variables inherited from run_blup90_pipe.sh via export
# TRAIT, TRAIT_DIR, plink_prefix
ml plink/1.90b3.39 
plink --cow --bfile ${plink_prefix} \
     -chr 1-26 \
    --keep ${TRAIT_DIR}/keep.txt \
    --maf 0.05 --geno 0.1 --mind 0.1 \
    --make-bed --keep-allele-order \
    --out ${TRAIT_DIR}/data.temp

plink --cow --bfile ${TRAIT_DIR}/data.temp \
    --recode A --make-bed --keep-allele-order \
    --out ${TRAIT_DIR}/${TRAIT}

cut -c1-150 ${TRAIT_DIR}/${TRAIT}.raw | head -11

# Extract IDs (col 2), pad to 6 chars
awk '{print $2}' ${TRAIT_DIR}/${TRAIT}.raw | sed '1d' > ${TRAIT_DIR}/Ids.temp
awk '{ print length }' ${TRAIT_DIR}/Ids.temp | sort | uniq
awk -v m=6 '{printf("%-6s\n", $0)}' ${TRAIT_DIR}/Ids.temp > ${TRAIT_DIR}/Ids2.temp

# Extract genotypes: drop first 6 cols, remove header, replace NA->5, remove spaces
cut -d' ' -f7- ${TRAIT_DIR}/${TRAIT}.raw | sed -e '1d' -e 's/NA/5/g' -e 's/[[:blank:]]//g' > ${TRAIT_DIR}/geno.temp
awk '{ print length }' ${TRAIT_DIR}/geno.temp | uniq

# Paste IDs and genotypes
paste ${TRAIT_DIR}/Ids2.temp ${TRAIT_DIR}/geno.temp -d " " > ${TRAIT_DIR}/${TRAIT}.geno

# SNP map
awk '{print $2,$1,$4}' ${TRAIT_DIR}/${TRAIT}.bim | sed -e '1i\SNPID CHR POS' > ${TRAIT_DIR}/${TRAIT}.map

# PED file
awk '{print $1,"0","0"}' ${TRAIT_DIR}/${TRAIT}.geno > ${TRAIT_DIR}/${TRAIT}.PED

# Combine logs
cat ${TRAIT_DIR}/data.temp.log ${TRAIT_DIR}/${TRAIT}.log > ${TRAIT_DIR}/log.file

# Cleanup
rm -f ${TRAIT_DIR}/*.nosex
# rm -f ${TRAIT_DIR}/*.log
rm -f ${TRAIT_DIR}/*.temp*
rm -f ${TRAIT_DIR}/data.temp.*
rm -f ${TRAIT_DIR}/Header.map
rm -f ${TRAIT_DIR}/*.raw