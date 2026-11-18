#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Test.it
#SBATCH --output=Test.it_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=7000mb


cd ${GBLUP_GWAS_DIR}
echo "Starting ${TRAIT}"

# Variance components must come from AIREML.SNP.sh
if grep -qE "YOUR_RES_VAR|YOUR_GEN_VAR" ${SETUP_DIR}/${TRAIT}/GBLUP.parameter; then
    echo "ERROR: GBLUP.parameter still has variance placeholders — run AIREML.SNP.sh first."
    exit 1
fi
echo "Using AIREML variances:"
grep -A1 -E "^RESIDUAL_VARIANCE|^\(CO\)VARIANCES" ${SETUP_DIR}/${TRAIT}/GBLUP.parameter
mkdir -p ${GBLUP_GWAS_DIR}/${TRAIT}
cd ${GBLUP_GWAS_DIR}/${TRAIT}
cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.geno .
cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.PED .
cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.map .
cp ${SETUP_DIR}/${TRAIT}/GBLUP.parameter .
cp ${SETUP_DIR}/${TRAIT}/phenotype.txt pheno.txt


ml blupf90/20190531 
ml R
WEIGHTS="others"
WEIGHTS="NONLINEAR"
weight_default="NONLINEAR"
  if [[ "$WEIGHTS" == "$weight_default" ]]
  then
    round=(1 2 3)
  else
    round=(1 2 3)
    fi
  for i in "${round[@]}"
  do
    cd ${GBLUP_GWAS_DIR}/${TRAIT}/
    mkdir -p "round_"$i                #create one folder per round
    cd "round_"$i
    echo "Start round_"$i
    cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.geno .
    cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.PED .
    cp ${SETUP_DIR}/${TRAIT}/${TRAIT}.map .
    cp ${SETUP_DIR}/${TRAIT}/GBLUP.parameter .
    cp ${SETUP_DIR}/${TRAIT}/phenotype.txt pheno.txt
      if (($i == 1))                                  #saving weights to calculate G matrix by blupf90
      then
        awk 'NR>1 {print "1",$1}' ${TRAIT}.map | awk '{print $1}' > wei.txt    #weights=1 for all SNPs for round_1
      else
        awk 'NR>1 {print $7}' ../round_$((i - 1))/snp_sol > wei.txt           #saving weights from snp_sol of previous round for round_2 and round_3
      fi
    cp ${SETUP_DIR}/${TRAIT}/GBLUP.parameter .
    sed -i GBLUP.parameter -e '$aOPTION weightedG wei.txt'
    printf GBLUP.parameter | renumf90 > renum.${i}.${TRAIT}.txt
    echo "renf90.par" | blupf90 > GBLUP.${i}.${TRAIT}.txt
    
      weight_default="NONLINEAR"
      if [[ "$WEIGHTS" == "$weight_default" ]]                                   #saving weights to calculate G matrix by blupf90
      then
        sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION missing -999/#OPTION missing -999/g" | sed "s/solution mean.*/windows_variance_mbp 1/g" | sed "s/sol se/which_weight nonlinearA 1.25/g" | sed -e '$aOPTION Manhattan_plot_R'  > POSTGS.Parameter        #to create plots add: sed -e '$aOPTION Manhattan_plot_R'    #weights=1 for all SNPs for round_1
      else
        sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION missing -999/#OPTION missing -999/g" | sed "s/solution mean.*/windows_variance_mbp 1/g" | sed "s/OPTION sol se/which_weight 1/g"  | sed -e '$aOPTION Manhattan_plot_R'  > POSTGS.Parameter        #to create plots add: sed -e '$aOPTION Manhattan_plot_R'
      fi
    sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION solution mean/#OPTION solution mean/g" | sed "s/OPTION missing -999/#OPTION missing -999/g"| sed "s/OPTION sol se/#OPTION sol se/g"  | sed -e '$aOPTION saveHinv' > PREGS.Parameter 
    printf PREGS.Parameter | preGSf90 > preGSf90.${i}.${TRAIT}.txt
    printf POSTGS.Parameter | postGSf90 > postGSf90.${i}.${TRAIT}.txt
    rm $TRAIT.geno
    Rscript Sft*.R
    Rscript Vft*.R
    Rscript Pft*.R
    echo "Finished round_"$i   
    cd ${GBLUP_GWAS_DIR}/${TRAIT}/
  done
