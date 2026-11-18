#!/bin/bash
# =============================================================================
# BLUPF90.GWAS.summarize.sh
#
# Summarizes GWAS results from a multi-round BLUPF90/postGSf90 run into a
# single CSV file with columns:
#   SNPID, CHR, POS, snp_effect, var_Explained, neg_log_pval
#
# Sources:
#   SNPID / CHR / POS  ->  AVE.map           (root dir; has header)
#   snp_effect         ->  round_3/chrsnp    (no header)
#   var_Explained      ->  round_3/chrsnpvar (no header)
#   -log(pval)         ->  round_1/chrsnp_pval (no header)
#
# Merge key: CHR + POS + within-position index (to handle duplicate positions)
#
# Window variance (round_3/windows_variance) is also parsed into:
#   GWAS_windows.<TRAIT>.csv               all non-empty windows
#   GWAS_sig_windows.VE<thr>.<TRAIT>.csv   windows >= VE_THRESHOLD % with
#                                          top SNPs + every SNP's info
# Set VE_THRESHOLD in the environment to change the cutoff (default 1).
#
# Step 3 runs find_QTL_genes_OarV4_fromWindows.R (same directory as this
# script): QTL called from round_3/windows_variance, with their variance in
# rounds 1-2 added to the table, a supplementary Manhattan figure with one
# row per round, and genes + KEGG/GO ORA.
# Optional environment overrides:
#   QTL_SCRIPT, QTL_OUT, ANNOT_DIR (GFF/TxDb cache), ARRAY_FILE,
#   RUN_ENRICHMENT (TRUE/FALSE), CLUSTER_GAP_MB, FLANK_KB,
#   QTL_ROUNDS (default "1,2,3"; QTL are called in the last one)
#
# Usage:
#   Called automatically from run_blup90.pipe.sh, or standalone:
#     bash BLUPF90.GWAS.summarize.sh
#
# Requires: R (with data.table; base merge also works — see R script below)
# =============================================================================

# ---------------------------------------------------------------------------
# Inherit environment from the pipeline if sourced from run_blup90.pipe.sh.
# When run standalone, fall back to sensible defaults.
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Root GWAS directory (where AVE.map lives and round_* subdirs are found)
GWAS_ROOT="${GBLUP_GWAS_DIR}/${TRAIT}"

# Output file
SUMMARY_OUT="${GWAS_ROOT}/GWAS_summary.${TRAIT:-TRAIT}.csv"
PLOT_DIR="${GWAS_ROOT}/plots"
echo "============================================================"
echo " BLUPF90 GWAS Summarizer"
echo "============================================================"
echo " GWAS root : ${GWAS_ROOT}"
echo " Output    : ${SUMMARY_OUT}"
echo "------------------------------------------------------------"
ml R/4.4
# ---------------------------------------------------------------------------
# Verify required input files exist
# ---------------------------------------------------------------------------
MAP_FILE="${GWAS_ROOT}/${TRAIT}.map"
SNP_EFFECT_FILE="${GWAS_ROOT}/round_3/chrsnp"
VAR_EXPL_FILE="${GWAS_ROOT}/round_3/chrsnpvar"
PVAL_FILE="${GWAS_ROOT}/round_1/chrsnp_pval"
WIN_VAR_FILE="${WIN_VAR_FILE:-${GWAS_ROOT}/round_3/windows_variance}"
VE_THRESHOLD="${VE_THRESHOLD:-1}"
WIN_OUT="${GWAS_ROOT}/GWAS_windows.${TRAIT:-TRAIT}.csv"
SIGWIN_OUT="${GWAS_ROOT}/GWAS_sig_windows.VE${VE_THRESHOLD}.${TRAIT:-TRAIT}.csv"

if [[ ! -f "$WIN_VAR_FILE" ]]; then
    echo "WARNING: windows_variance not found ($WIN_VAR_FILE) — window summary skipped."
    WIN_VAR_FILE="NONE"
fi

for f in "$MAP_FILE" "$SNP_EFFECT_FILE" "$VAR_EXPL_FILE" "$PVAL_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

echo " All required input files found."
echo "------------------------------------------------------------"

# ---------------------------------------------------------------------------
# Run R to merge and write the summary CSV
# ---------------------------------------------------------------------------
Rscript - \
    "$MAP_FILE" \
    "$SNP_EFFECT_FILE" \
    "$VAR_EXPL_FILE" \
    "$PVAL_FILE" \
    "$SUMMARY_OUT" \
    "$WIN_VAR_FILE" \
    "$VE_THRESHOLD" \
    "$WIN_OUT" \
    "$SIGWIN_OUT" \
<<'EOF'

args <- commandArgs(trailingOnly = TRUE)
map_file       <- args[1]
snpeff_file    <- args[2]
varexpl_file   <- args[3]
pval_file      <- args[4]
out_file       <- args[5]
win_file       <- args[6]
ve_threshold   <- as.numeric(args[7])
win_out        <- args[8]
sigwin_out     <- args[9]

cat("Reading input files...\n")

# ---------------------------------------------------------------------------
# 1. AVE.map — has a header: SNPID CHR POS
# ---------------------------------------------------------------------------
map <- read.table(map_file, header = TRUE, stringsAsFactors = FALSE)
colnames(map) <- c("SNPID", "CHR", "POS")

# Build a within-position index so duplicate CHR+POS rows can be matched to
# the correct SNP name (e.g. two SNPs at the same physical position).


# ---------------------------------------------------------------------------
# 2. chrsnp (round_3) — no header
#    Columns: TRAIT_POS  EFFECT  snp_effect  index  CHR  POS
# ---------------------------------------------------------------------------
snpeff <- read.table(snpeff_file, header = FALSE, stringsAsFactors = FALSE)
colnames(snpeff) <- c("TRAIT_POS", "EFFECT", "snp_effect", "index", "CHR", "POS")
snpeff <- snpeff[, c("CHR", "POS", "index", "snp_effect")]
snpeff$index <- as.integer(snpeff$index)

# ---------------------------------------------------------------------------
# 3. chrsnpvar (round_3) — no header
#    Columns: TRAIT_POS  EFFECT  var_Explained  index  CHR  POS
# ---------------------------------------------------------------------------
varexpl <- read.table(varexpl_file, header = FALSE, stringsAsFactors = FALSE)
colnames(varexpl) <- c("TRAIT_POS", "EFFECT", "var_Explained", "index", "CHR", "POS")
varexpl <- varexpl[, c("CHR", "POS", "index", "var_Explained")]
varexpl$index <- as.integer(varexpl$index)

# ---------------------------------------------------------------------------
# 4. chrsnp_pval (round_1) — no header
#    Columns: TRAIT_POS  EFFECT  -log(pval)  index  CHR  POS
# ---------------------------------------------------------------------------
pval <- read.table(pval_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pval) <- c("TRAIT_POS", "EFFECT", "neg_log_pval", "index", "CHR", "POS")
pval <- pval[, c("CHR", "POS", "index", "neg_log_pval")]
pval$index <- as.integer(pval$index)

# ---------------------------------------------------------------------------
# 5. Merge — spine is the map; left-join everything else by CHR + POS + index
# ---------------------------------------------------------------------------
cat("Merging tables...\n")

merged <- merge(map, snpeff,  by = c("CHR", "POS"), all.x = TRUE)
print(head(merged))
merged <- merge(merged, varexpl, by = c("CHR", "POS", "index"), all.x = TRUE)
merged <- merge(merged, pval,    by = c("CHR", "POS", "index"), all.x = TRUE)

# ---------------------------------------------------------------------------
# 6. Tidy up and write output
# ---------------------------------------------------------------------------
merged <- merged[order(merged$CHR, merged$POS, merged$index), ]
merged <- merged[, c("SNPID", "CHR", "POS", "snp_effect", "var_Explained", "neg_log_pval")]

# Report any SNPs that had no match in the GWAS result files
n_missing_eff  <- sum(is.na(merged$snp_effect))
n_missing_var  <- sum(is.na(merged$var_Explained))
n_missing_pval <- sum(is.na(merged$neg_log_pval))

if (n_missing_eff  > 0) cat("  WARNING:", n_missing_eff,  "SNPs missing snp_effect\n")
if (n_missing_var  > 0) cat("  WARNING:", n_missing_var,  "SNPs missing var_Explained\n")
if (n_missing_pval > 0) cat("  WARNING:", n_missing_pval, "SNPs missing neg_log_pval\n")

write.csv(merged, file = out_file, row.names = FALSE, quote = FALSE)

cat(sprintf("Done. %d SNPs written to: %s\n", nrow(merged), out_file))

# ---------------------------------------------------------------------------
# 7. windows_variance (round_3) — no header
#    Columns: trait  effect  start_SNP  end_SNP  n_SNPs  CHR_start  CHR_end
#             CHR_windowIndex  %var
#    Empty end-of-chromosome windows (0 SNPs, start > end) are dropped.
# ---------------------------------------------------------------------------
if (win_file != "NONE") {
  cat("Summarizing windows_variance...\n")
  w <- read.table(win_file, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("Trait_No", "Effect", "Start_SNP", "End_SNP", "n_SNPs",
                                "Start_ID", "End_ID", "Window_ID", "VE"))
  w <- w[w$n_SNPs > 0, ]
  w$CHR       <- sub("_.*$", "", w$Start_ID)
  w$Win_Start <- as.numeric(sub("^.*_", "", w$Start_ID))
  w$Win_End   <- as.numeric(sub("^.*_", "", w$End_ID))
  w$Win_Size_kb <- (w$Win_End - w$Win_Start) / 1e3
  w <- w[, c("Window_ID", "CHR", "Win_Start", "Win_End", "Win_Size_kb",
             "n_SNPs", "Start_SNP", "End_SNP", "VE")]
  write.csv(w, win_out, row.names = FALSE, quote = FALSE)
  cat(sprintf("  %d windows (total VE = %.2f%%) -> %s\n", nrow(w), sum(w$VE), win_out))

  sig <- w[w$VE >= ve_threshold, ]
  sig <- sig[order(-sig$VE), ]
  if (nrow(sig) == 0) {
    cat(sprintf("  No windows with VE >= %s%%\n", ve_threshold))
  } else {
    # every SNP inside each significant window, with window info attached
    snps <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) {
      s <- merged[as.character(merged$CHR) == sig$CHR[i] &
                  merged$POS >= sig$Win_Start[i] & merged$POS <= sig$Win_End[i], ]
      if (nrow(s) == 0) return(NULL)
      cbind(sig[rep(i, nrow(s)), ], s[, c("SNPID", "POS", "snp_effect",
                                          "var_Explained", "neg_log_pval")])
    }))
    names(snps)[names(snps) == "POS"] <- "SNP_POS"
    snps$Top_P_in_window <- ave(snps$neg_log_pval, snps$Window_ID,
                                FUN = function(x) x == max(x, na.rm = TRUE))
    write.csv(snps, sigwin_out, row.names = FALSE, quote = FALSE)
    cat(sprintf("  %d windows >= %s%% (%d SNPs) -> %s\n",
                nrow(sig), ve_threshold, nrow(snps), sigwin_out))
    print(sig[, c("Window_ID", "CHR", "Win_Start", "Win_End", "n_SNPs", "VE")],
          row.names = FALSE)
  }
}
EOF

echo "============================================================"
echo " Summary complete: ${SUMMARY_OUT}"
[[ "$WIN_VAR_FILE" != "NONE" ]] && echo " Windows         : ${WIN_OUT}" && echo " Sig windows     : ${SIGWIN_OUT}"
echo "============================================================"


# ---------------------------------------------------------------------------
# Step 2: Generate plots
# ---------------------------------------------------------------------------
echo " [2/3] Generating GWAS plots..."
PLOT_SCRIPT="${SCRIPT_DIR}/BLUPF90.GWAS.plots.R"
if [[ -f "$PLOT_SCRIPT" ]]; then
    Rscript "$PLOT_SCRIPT" "$SUMMARY_OUT" "$PLOT_DIR" "${TRAIT:-TRAIT}"
else
    echo "WARNING: Plotting script not found at: ${PLOT_SCRIPT} — plots skipped."
    echo "         Place BLUPF90.GWAS.plots.R in the same directory as this script."
fi

# ---------------------------------------------------------------------------
# Step 3: Call QTL from round-3 windows_variance (+ genes / ORA)
# ---------------------------------------------------------------------------
echo "------------------------------------------------------------"
echo " [3/3] Calling QTL (round-3 windows >= ${VE_THRESHOLD}% genetic variance)..."
QTL_SCRIPT="${QTL_SCRIPT:-${SCRIPT_DIR}/find_QTL_genes_OarV4_fromWindows.R}"
QTL_OUT="${QTL_OUT:-${GWAS_ROOT}/QTL_results_VE${VE_THRESHOLD}}"
ANNOT_DIR="${ANNOT_DIR:-${GBLUP_GWAS_DIR}/annotation}"      # cached GFF + TxDb (shared by traits)
ARRAY_FILE="/blue/frezende/gzayas97/For_UF/Sheep_Project/data/UFL_Sheep_MapInfo_Oar4.xlsx"
RUN_ENRICHMENT="${RUN_ENRICHMENT:-TRUE}"
CLUSTER_GAP_MB="${CLUSTER_GAP_MB:-1}"
FLANK_KB="${FLANK_KB:-0}"
QTL_ROUNDS="${QTL_ROUNDS:-1,2,3}"   # reads round_<k>/windows_variance; QTL called in the last
for r in ${QTL_ROUNDS//,/ }; do
    [[ -f "${GWAS_ROOT}/round_${r}/windows_variance" ]] || \
        echo "WARNING: ${GWAS_ROOT}/round_${r}/windows_variance missing — round ${r} will be skipped."
done

QTL_STATUS="skipped"
if [[ "$WIN_VAR_FILE" == "NONE" ]]; then
    echo "WARNING: no windows_variance file — QTL calling skipped."
elif [[ ! -f "$QTL_SCRIPT" ]]; then
    echo "WARNING: QTL script not found at: ${QTL_SCRIPT} — QTL calling skipped."
else
    if [[ "$RUN_ENRICHMENT" == "TRUE" && ! -f "$ARRAY_FILE" ]]; then
        echo "WARNING: array map not found (${ARRAY_FILE}) — running without gene/ORA step."
        RUN_ENRICHMENT="FALSE"
    fi
    echo " Script     : ${QTL_SCRIPT}"
    echo " Output     : ${QTL_OUT}"
    echo " Enrichment : ${RUN_ENRICHMENT}"
    if Rscript "$QTL_SCRIPT" \
        --summary="$SUMMARY_OUT" \
        --win_dir="$GWAS_ROOT" \
        --rounds="$QTL_ROUNDS" \
        --trait="${TRAIT:-TRAIT}" \
        --out="$QTL_OUT" \
        --ve="$VE_THRESHOLD" \
        --gap_mb="$CLUSTER_GAP_MB" \
        --flank_kb="$FLANK_KB" \
        --enrich="$RUN_ENRICHMENT" \
        --array="$ARRAY_FILE" \
        --annot_dir="$ANNOT_DIR"; then
        QTL_STATUS="done"
    else
        echo "WARNING: QTL script exited with an error (see messages above)."
        QTL_STATUS="failed"
    fi
fi

echo "============================================================"
echo " All outputs ready."
echo "   Summary : ${SUMMARY_OUT}"
echo "   Windows : ${WIN_OUT}"
echo "   Plots   : ${PLOT_DIR}/"
echo "   QTL     : ${QTL_OUT}/  (${QTL_STATUS})"
echo "============================================================"