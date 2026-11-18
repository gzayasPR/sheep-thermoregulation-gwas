#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=PopGen
#SBATCH --output=PopGen_%j.out
#SBATCH --error=PopGen_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
# =============================================================================
# run_popgen.pipe.sh
#
# Population structure from the same PLINK genotypes used for the GWAS:
#   1. popgen.prep.sh      -> QC + LD pruning (+ optional relatedness filter)
#   2. ADMIXTURE.sh        -> ADMIXTURE for K = K_MIN..K_MAX with CV error
#   3. PCA.sh              -> PLINK PCA on the pruned set
#   4. popgen.plots.R      -> PCA plots, scree, CV-error plot, ADMIXTURE bars
# =============================================================================

env_file="/blue/frezende/gzayas97/For_UF/Sheep_Project/bin/project_env.sh"
source $env_file

# --- Inputs (same genotypes / phenotype file as the GWAS) ---
export PHENO="/blue/frezende/gzayas97/For_UF/Sheep_Project/data/UFL_Sheep_master_10.29.2026_matched.csv"
export plink_prefix="${my_data}/UFL_Sheep_plink"
export POPGEN_SCRIPTS="${my_bin}/scripts/PopGen"
export model_csv="${my_bin}/scripts/BLUPF90/sheep_model.csv"   # used to find the GenoID column

# --- Species / chromosome set ---
export SPECIES="sheep"          # PLINK species flag (--sheep); use "cow" for cattle
export AUTOSOMES="1-26"         # ovine autosomes

# --- Sample subset (optional) ---
# Leave empty to use every genotyped animal.
# To match the GWAS animals exactly, point to the keep file from setup, e.g.
#   export KEEP_FILE="${my_results}/BLUPF90/Setup/MAX/keep.txt"
export KEEP_FILE="/blue/frezende/gzayas97/For_UF/Sheep_Project/data/UFL_Sheep_master_10.29.2026_keep.txt"

# --- QC ---
export MAF=0.05
export GENO=0.1
export MIND=0.1

# --- LD pruning (window kb-free: SNP count, step, r2) ---
export LD_WIN=50
export LD_STEP=5
export LD_R2=0.2

# --- Relatedness filter (optional; e.g. 0.125 removes up to ~1st-cousin pairs) ---
# Leave empty to keep all animals. Close relatives can create spurious
# ADMIXTURE clusters / PCs driven by families rather than ancestry.
export REL_CUTOFF=""

# --- ADMIXTURE ---
export K_MIN=2
export K_MAX=10
export CV_FOLDS=10
export SEED=12345
export THREADS=${SLURM_CPUS_PER_TASK:-4}

# --- PCA ---
export N_PCS=20

# --- Plot grouping (optional) ---
# Name from sheep_model.csv (or a PHENO header) used to color PCA points and
# order ADMIXTURE bars (e.g. "BREED", "FARM"). Leave empty to color by
# the dominant ADMIXTURE cluster at the best K.
export GROUP_COL="BREED_LABEL"
export GROUP_COL2="FARM"   # a Name in model.csv or a header in the pheno file
# --- Output dirs ---
export out_dir="${my_results}/PopGen"
export POP_PREFIX="popgen"
export PREP_DIR="${out_dir}/data"
export ADMIX_DIR="${out_dir}/admixture"
export PCA_DIR="${out_dir}/pca"
export PLOT_DIR="${out_dir}/plots"

mkdir -p $out_dir $PREP_DIR $ADMIX_DIR $PCA_DIR $PLOT_DIR

bash "${POPGEN_SCRIPTS}/popgen.prep.sh"   

bash "${POPGEN_SCRIPTS}/ADMIXTURE.sh"     

bash "${POPGEN_SCRIPTS}/PCA.sh"       

ml R
Rscript "${POPGEN_SCRIPTS}/popgen.plots.R"
