#!/bin/bash
# =============================================================================
# ADMIXTURE.sh
#
# Runs ADMIXTURE for K = K_MIN..K_MAX on the LD-pruned set with
# cross-validation, then collects CV error per K.
# All variables inherited from run_popgen.pipe.sh via export:
#   PREP_DIR, ADMIX_DIR, POP_PREFIX, K_MIN, K_MAX, CV_FOLDS, SEED, THREADS
#
# Output (in ${ADMIX_DIR}):
#   ${POP_PREFIX}.pruned.K.Q / .K.P   ancestry proportions / allele freqs
#   log.K<K>.out                      ADMIXTURE log per K
#   CV_error.txt                      K  CV_error
# =============================================================================

ml admixture
cd ${ADMIX_DIR}

echo "============================================================"
echo " ADMIXTURE: K = ${K_MIN}..${K_MAX}, ${CV_FOLDS}-fold CV, ${THREADS} threads"
echo "============================================================"

# ADMIXTURE writes outputs to the working directory, named after the .bed
for ext in bed bim fam; do
    ln -sf ${PREP_DIR}/${POP_PREFIX}.pruned.${ext} .
done

for K in $(seq ${K_MIN} ${K_MAX}); do
    echo " Starting K=${K}"
    admixture --cv=${CV_FOLDS} -s ${SEED} -j${THREADS} \
        ${POP_PREFIX}.pruned.bed ${K} > log.K${K}.out 2>&1
    if [[ $? -ne 0 ]]; then
        echo " WARNING: ADMIXTURE failed for K=${K}; see log.K${K}.out"
    fi
done

# Collect CV error: lines look like "CV error (K=3): 0.51234"
echo "K CV_error" > CV_error.txt
grep -h "CV error" log.K*.out \
    | sed -E 's/CV error \(K=([0-9]+)\): ([0-9.eE+-]+)/\1 \2/' \
    | sort -n -k1,1 >> CV_error.txt

echo "------------------------------------------------------------"
cat CV_error.txt
BEST_K=$(awk 'NR>1' CV_error.txt | sort -g -k2,2 | head -1 | awk '{print $1}')
echo " Lowest CV error at K=${BEST_K}"
echo "============================================================"
