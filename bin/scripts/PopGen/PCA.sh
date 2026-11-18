#!/bin/bash
# =============================================================================
# PCA.sh
#
# PCA on the LD-pruned genotypes with PLINK 1.9 (GRM-based).
# All variables inherited from run_popgen.pipe.sh via export:
#   SPECIES, PREP_DIR, PCA_DIR, POP_PREFIX, N_PCS, THREADS
#
# Output (in ${PCA_DIR}):
#   ${POP_PREFIX}.pca.eigenvec   FID IID PC1..PCn (with header)
#   ${POP_PREFIX}.pca.eigenval   eigenvalues
# =============================================================================

ml plink/1.90b3.39
cd ${PCA_DIR}

echo "============================================================"
echo " PCA: top ${N_PCS} PCs"
echo "============================================================"

plink --${SPECIES} --bfile ${PREP_DIR}/${POP_PREFIX}.pruned \
    --pca ${N_PCS} header tabs \
    --threads ${THREADS} \
    --out ${POP_PREFIX}.pca || exit 1

rm -f *.nosex
echo " PCA done: ${PCA_DIR}/${POP_PREFIX}.pca.eigenvec"
echo "============================================================"
