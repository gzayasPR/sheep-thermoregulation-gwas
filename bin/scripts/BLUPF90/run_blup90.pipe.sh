#!/bin/bash

env_file="/blue/frezende/gzayas97/For_UF/Sheep_Project/bin/project_env.sh"
source $env_file

# --- Trait & paths ---
export TRAIT="AVE"
export PHENO="${my_data}/UFL_Sheep_master_10.29.2026_matched_pcs.csv"
export plink_prefix="${my_data}/UFL_Sheep_plink"
export BLUPF90_SCRIPTS="${my_bin}/scripts/BLUPF90"

# --- Parameter files ---
export AIREML_param="${BLUPF90_SCRIPTS}/AIREML.parameter"
export GBLUP_param="${BLUPF90_SCRIPTS}/GBLUP.parameter"
export model_csv="${BLUPF90_SCRIPTS}/sheep_model.csv"

# --- Model ---
export MODEL_COVAR="~ REP + PC1 + PC2"
ARRAY_FILE="/blue/frezende/gzayas97/For_UF/Sheep_Project/data/UFL_Sheep_MapInfo_Oar4.xlsx"
# --- Output dirs ---
export out_dir="${my_results}/BLUPF90"
export SETUP_DIR="${out_dir}/Setup"
export AIREML_DIR="${out_dir}/aireml"
export GBLUP_GWAS_DIR="${out_dir}/GBLUP_GWAS"
export TRAIT_DIR="${SETUP_DIR}/${TRAIT}"

mkdir -p $out_dir $AIREML_DIR $GBLUP_GWAS_DIR $SETUP_DIR $TRAIT_DIR

# bash "${BLUPF90_SCRIPTS}/BLUPF90.setup.geno.sh"

# bash "${BLUPF90_SCRIPTS}/AIREML.SNP.sh"

# bash "${BLUPF90_SCRIPTS}/BLUPF90.GWAS.geno.add.sh"
RUN_ENRICHMENT=TRUE
bash "${BLUPF90_SCRIPTS}/BLUPF90.GWAS.summarize.sh"