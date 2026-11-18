#!/bin/bash

# All variables inherited from run_blup90_pipe.sh via export
# TRAIT, TRAIT_DIR, SETUP_DIR, BLUPF90_SCRIPTS, PHENO, model_csv, MODEL_COVAR, plink_prefix

cd ${SETUP_DIR}
mkdir -p ${TRAIT_DIR}

ml R
Rscript "${BLUPF90_SCRIPTS}/setup_pheno_geno.R" \
        ${PHENO} \
        ${model_csv} \
        "${MODEL_COVAR}" \
        ${plink_prefix}.fam \
        $TRAIT \
        ${TRAIT_DIR} 

ml plink
bash "${BLUPF90_SCRIPTS}/plink_2_blupf90.sh"


TRAITS_HEAD=$(echo "#`head -1 $PHENO`")
    sed -e "s/#YOUR_TRAITS/${TRAIT}/g" \
    -e "s/YOUR_PED/${TRAIT}.PED/g" \
    -e "s/YOUR_GENO/${TRAIT}.geno/g" \
    -e "s/YOUR_MAP/${TRAIT}.map/g" \
    $AIREML_param > ${TRAIT_DIR}/AIREML.parameter    
 
TRAITS_HEAD=$(echo "#`head -1 $PHENO`")
    sed -e "s/#YOUR_TRAITS/${TRAIT}/g" \
    -e "s/YOUR_PED/${TRAIT}.PED/g" \
    -e "s/YOUR_GENO/${TRAIT}.geno/g" \
    -e "s/YOUR_MAP/${TRAIT}.map/g" \
   $GBLUP_param > ${TRAIT_DIR}/GBLUP.parameter    