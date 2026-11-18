#!/bin/bash
# =============================================================================
# popgen.prep.sh
#
# QC + LD pruning of the GWAS genotypes for ADMIXTURE and PCA.
# All variables inherited from run_popgen.pipe.sh via export:
#   plink_prefix, SPECIES, AUTOSOMES, KEEP_FILE, MAF, GENO, MIND,
#   LD_WIN, LD_STEP, LD_R2, REL_CUTOFF, PREP_DIR, POP_PREFIX
#
# Output: ${PREP_DIR}/${POP_PREFIX}.pruned.{bed,bim,fam}
# =============================================================================

ml plink/1.90b3.39 
cd ${PREP_DIR}

echo "============================================================"
echo " PopGen prep: QC + LD pruning"
echo "============================================================"

KEEP_ARG=""
if [[ -n "${KEEP_FILE}" ]]; then
    if [[ ! -f "${KEEP_FILE}" ]]; then
        echo "ERROR: KEEP_FILE not found: ${KEEP_FILE}"; exit 1
    fi
    KEEP_ARG="--keep ${KEEP_FILE}"
    echo " Restricting to animals in: ${KEEP_FILE}"
fi

# --- 1. Autosomes + QC (unnamed SNPs get CHR:POS IDs so pruning works) ---
plink --${SPECIES} --bfile ${plink_prefix} \
    --chr ${AUTOSOMES} ${KEEP_ARG} \
    --set-missing-var-ids @:# \
    --maf ${MAF} --geno ${GENO} --mind ${MIND} \
    --make-bed --keep-allele-order \
    --out ${POP_PREFIX}.qc || exit 1

# --- 2. LD pruning ---
plink --${SPECIES} --bfile ${POP_PREFIX}.qc \
    --indep-pairwise ${LD_WIN} ${LD_STEP} ${LD_R2} \
    --out ${POP_PREFIX}.prune || exit 1

plink --${SPECIES} --bfile ${POP_PREFIX}.qc \
    --extract ${POP_PREFIX}.prune.prune.in \
    --make-bed --keep-allele-order \
    --out ${POP_PREFIX}.pruned.tmp || exit 1

# --- 3. Optional relatedness filter (GRM on pruned SNPs) ---
if [[ -n "${REL_CUTOFF}" ]]; then
    echo " Removing one of each pair with relatedness > ${REL_CUTOFF}"
    plink --${SPECIES} --bfile ${POP_PREFIX}.pruned.tmp \
        --rel-cutoff ${REL_CUTOFF} \
        --out ${POP_PREFIX}.unrel || exit 1
    plink --${SPECIES} --bfile ${POP_PREFIX}.pruned.tmp \
        --keep ${POP_PREFIX}.unrel.rel.id \
        --make-bed --keep-allele-order \
        --out ${POP_PREFIX}.pruned || exit 1
else
    for ext in bed bim fam; do
        mv ${POP_PREFIX}.pruned.tmp.${ext} ${POP_PREFIX}.pruned.${ext}
    done
fi

rm -f ${POP_PREFIX}.pruned.tmp.* *.nosex

echo "------------------------------------------------------------"
echo " Animals after QC : $(wc -l < ${POP_PREFIX}.qc.fam)"
echo " SNPs after QC    : $(wc -l < ${POP_PREFIX}.qc.bim)"
echo " Animals (final)  : $(wc -l < ${POP_PREFIX}.pruned.fam)"
echo " SNPs (pruned)    : $(wc -l < ${POP_PREFIX}.pruned.bim)"
echo "============================================================"
