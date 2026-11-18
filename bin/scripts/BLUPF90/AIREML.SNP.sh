#!/bin/bash
# All variables inherited from run_blup90.pipe.sh via export
# TRAIT, SETUP_DIR, AIREML_DIR
#
# Runs airemlf90 (and remlf90 as a backup), then writes the estimated
# genetic and residual variances into ${SETUP_DIR}/${TRAIT}/GBLUP.parameter
# so the GBLUP/GWAS step (BLUPF90.GWAS.geno.add.sh) uses them.

mkdir -p ${AIREML_DIR}/${TRAIT}
cd ${AIREML_DIR}/${TRAIT}
echo "Starting ${TRAIT}"

cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.geno .
cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.PED .
cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.map .
cp ${SETUP_DIR}/${TRAIT}/AIREML.parameter .
cp ${SETUP_DIR}/${TRAIT}/phenotype.txt pheno.txt

ml blupf90/20190531
echo "Starting AIREML"
renumf90 AIREML.parameter > renum.out.txt
airemlf90 renf90.par > ${TRAIT}.airemlf90.txt
grep "Final Estimates" -A 40 ${TRAIT}.airemlf90.txt
echo "REML"
remlf90 renf90.par > ${TRAIT}.remlf90.txt
grep "Final Estimates" -A 40 ${TRAIT}.remlf90.txt

# ---------------------------------------------------------------------------
# Extract variance components from the LAST "Final Estimates" block.
# Single-trait, one random effect:
#   Genetic variance(s) for effect  N
#     <Vg>
#   Residual variance(s)
#     <Ve>
# ---------------------------------------------------------------------------
get_vc() {   # $1 = blupf90 stdout file ; prints "Vg Ve" or nothing
    awk '
        /Final Estimates/                 { fe = 1; g = ""; r = ""; next }
        fe && /Genetic variance\(s\)/     { getline; g = $1; next }
        fe && /Residual variance\(s\)/    { getline; r = $1; next }
        END { if (g != "" && r != "") print g, r }
    ' "$1"
}

VC_SOURCE="airemlf90"
read VG VE < <(get_vc ${TRAIT}.airemlf90.txt)
if [[ -z "$VG" || -z "$VE" ]]; then
    echo "WARNING: no Final Estimates in airemlf90 output — falling back to remlf90"
    VC_SOURCE="remlf90"
    read VG VE < <(get_vc ${TRAIT}.remlf90.txt)
fi
if [[ -z "$VG" || -z "$VE" ]]; then
    echo "ERROR: could not extract variance components for ${TRAIT}."
    echo "       GBLUP.parameter was NOT updated; the GWAS step will stop."
    exit 1
fi

H2=$(awk -v g="$VG" -v r="$VE" 'BEGIN { printf "%.4f", g / (g + r) }')
echo "------------------------------------------------------------"
echo " ${TRAIT} variance components (${VC_SOURCE})"
echo "   Vg = ${VG}   Ve = ${VE}   h2 = ${H2}"
echo "------------------------------------------------------------"
printf "TRAIT\tSOURCE\tVg\tVe\th2\n%s\t%s\t%s\t%s\t%s\n" \
    "$TRAIT" "$VC_SOURCE" "$VG" "$VE" "$H2" > ${TRAIT}.varcomp.txt

# ---------------------------------------------------------------------------
# Write them into the GBLUP parameter file used by the GWAS step.
# Replaces the line after RESIDUAL_VARIANCE and after (CO)VARIANCES, so it
# works on the YOUR_RES_VAR / YOUR_GEN_VAR placeholders and on reruns.
# ---------------------------------------------------------------------------
GBLUP_TRAIT_PAR="${SETUP_DIR}/${TRAIT}/GBLUP.parameter"
awk -v g="$VG" -v r="$VE" '
    set_r { print r; set_r = 0; next }
    set_g { print g; set_g = 0; next }
    /^RESIDUAL_VARIANCE/ { set_r = 1 }
    /^\(CO\)VARIANCES/   { set_g = 1 }
    { print }
' "$GBLUP_TRAIT_PAR" > "${GBLUP_TRAIT_PAR}.tmp" && mv "${GBLUP_TRAIT_PAR}.tmp" "$GBLUP_TRAIT_PAR"

echo "Updated ${GBLUP_TRAIT_PAR}:"
grep -A1 -E "^RESIDUAL_VARIANCE|^\(CO\)VARIANCES" "$GBLUP_TRAIT_PAR"
