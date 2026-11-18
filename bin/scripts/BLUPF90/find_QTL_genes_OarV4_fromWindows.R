##############################################################################
# find_QTL_genes_OarV4_fromWindows.R
# Purpose : Build QTL regions from postGSf90 `windows_variance` output (the
#           real window boundaries/variance from BLUPF90), keep windows that
#           explain >= VE_THRESHOLD % of the genetic variance, attach ALL
#           SNP-level GWAS information (effect, SNP var, -log10 p) inside each
#           window, then run gene retrieval + ORA (KEGG + GO BP/MF/CC).
#           Ovis aries Oar_v4.0 (GCF_000298735.2) — NCBI annotation
# Author  : Gabe Zayas Santiago
# Date    : 2026-11-04
#
# Inputs per trait:
#   - GWAS_summary.<TRAIT>.csv  (from BLUPF90.GWAS.summarize.sh)
#       SNPID, CHR, POS, snp_effect, var_Explained, neg_log_pval
#   - windows_variance          (postGSf90 output, one per round; no header)
#       read from <Win_Dir>/round_<k>/windows_variance for k in ROUNDS.
#       QTL are called in the LAST round; earlier rounds are reported
#       alongside (table columns + one figure row per round).
#       trait  effect  start_SNP  end_SNP  n_SNPs  CHR_startPOS  CHR_endPOS
#       CHR_windowIndex  %var
#
# Usage   : (a) interactively / locally — edit the USER SETTINGS block, or
#               Rscript find_QTL_genes_OarV4_fromWindows.R 0.5   (threshold only)
#           (b) from BLUPF90.GWAS.summarize.sh (one trait per call):
#               Rscript find_QTL_genes_OarV4_fromWindows.R \
#                 --summary=GWAS_summary.AVE.csv --win_dir=<GWAS_ROOT> --rounds=1,2,3 \
#                 --trait=AVE --out=QTL_results_VE1 --ve=1 \
#                 --array=UFL_Sheep_MapInfo_Oar4.xlsx --annot_dir=annotation \
#                 --enrich=TRUE --gap_mb=1 --flank_kb=0
##############################################################################

# ── 1. USER SETTINGS (local/interactive defaults) ────────────────────────────
proj_dir   <- "C:/Users/gabbo/OneDrive - UW-Madison/Sheep/Results_11_4_2026"
ARRAY_FILE <- "C:/Users/gabbo/OneDrive - UW-Madison/Sheep/UFL_Sheep_MapInfo_Oar4.xlsx"
ANNOT_DIR  <- proj_dir  # where Oar_v4.0.gff.gz / Oar_v4.0_txdb.sqlite are cached

# One row per trait: label, GWAS summary CSV, folder holding round_<k>/ subdirs
TRAITS <- data.frame(
  Trait   = c("Max BT",                      "Average BT"),
  Summary = c("Max BT/GWAS_summary.MAX.csv", "Average BT/GWAS_summary.AVE.csv"),
  Win_Dir = c("Max BT",                      "Average BT"),
  stringsAsFactors = FALSE
)
ROUNDS      <- 1:3                          # rounds to read; QTL called in max(ROUNDS)
WIN_PATTERN <- "round_%d/windows_variance"  # path of each round's file inside Win_Dir
SNPVAR_PATTERN <- "round_%d/chrsnpvar"      # per-SNP var explained, used for the figures
FIG_ANNOTATE   <- TRUE   # FALSE = figure rows exactly like manhattan_var_explained

VE_THRESHOLD   <- 1.0   # keep windows explaining >= this % of genetic variance
CLUSTER_GAP_MB <- 1     # sig windows on the same CHR closer than this -> one QTL
FLANK_KB       <- 0     # extend each QTL by this many kb on both sides for gene search
RUN_ENRICHMENT <- TRUE  # FALSE = windows/QTL/SNP tables + plots only (no Bioconductor)
out_dir        <- NULL  # NULL = timestamped folder under proj_dir/enrichment_results

# ── 1a. COMMAND-LINE OVERRIDES ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
kv   <- grep("^--[A-Za-z_]+=", args, value = TRUE)
opt  <- setNames(sub("^--[A-Za-z_]+=", "", kv), sub("^--([A-Za-z_]+)=.*$", "\\1", kv))
if (length(opt)) {
  need <- c("summary", "win_dir", "trait", "out")
  if (!all(need %in% names(opt)))
    stop("Pipeline mode needs --summary, --win_dir, --trait and --out")
  TRAITS  <- data.frame(Trait = opt[["trait"]], Summary = normalizePath(opt[["summary"]]),
                        Win_Dir = normalizePath(opt[["win_dir"]]), stringsAsFactors = FALSE)
  if ("rounds"      %in% names(opt)) ROUNDS      <- as.integer(strsplit(opt[["rounds"]], ",")[[1]])
  if ("win_pattern" %in% names(opt)) WIN_PATTERN <- opt[["win_pattern"]]
  if ("annotate"    %in% names(opt)) FIG_ANNOTATE <- toupper(opt[["annotate"]]) %in% c("TRUE", "T", "1", "YES")
  out_dir <- opt[["out"]]
  if ("ve"        %in% names(opt)) VE_THRESHOLD   <- as.numeric(opt[["ve"]])
  if ("gap_mb"    %in% names(opt)) CLUSTER_GAP_MB <- as.numeric(opt[["gap_mb"]])
  if ("flank_kb"  %in% names(opt)) FLANK_KB       <- as.numeric(opt[["flank_kb"]])
  if ("enrich"    %in% names(opt)) RUN_ENRICHMENT <- toupper(opt[["enrich"]]) %in% c("TRUE", "T", "1", "YES")
  if ("array"     %in% names(opt)) ARRAY_FILE     <- opt[["array"]]
  if ("annot_dir" %in% names(opt)) ANNOT_DIR      <- opt[["annot_dir"]]
  if (RUN_ENRICHMENT && !file.exists(ARRAY_FILE))
    stop("Array map not found: ", ARRAY_FILE, " (pass --array=, or --enrich=FALSE)")
} else {
  setwd(proj_dir)
  if (length(args) >= 1 && !is.na(suppressWarnings(as.numeric(args[1]))))
    VE_THRESHOLD <- as.numeric(args[1])
}

if (is.null(out_dir))
  out_dir <- file.path(proj_dir, "enrichment_results",
                       sprintf("%s_VE%s", format(Sys.time(), "%Y%m%d_%H%M%S"),
                               format(VE_THRESHOLD)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Trait(s)     :", paste(TRAITS$Trait, collapse = ", "), "\n")
CALL_ROUND <- max(ROUNDS)
cat("Rounds       :", paste(ROUNDS, collapse = ", "), "(QTL called in round", CALL_ROUND, ")\n")
cat("VE threshold :", VE_THRESHOLD, "%\n")
cat("Cluster gap  :", CLUSTER_GAP_MB, "Mb | Flank:", FLANK_KB, "kb\n")
cat("Enrichment   :", RUN_ENRICHMENT, "\n")
cat("Output dir   :", out_dir, "\n\n")

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
core_pkgs <- c("dplyr", "ggplot2", "writexl")
bio_pkgs  <- c("txdbmaker", "GenomicFeatures", "GenomicRanges", "GenomeInfoDb",
               "AnnotationDbi", "clusterProfiler", "org.Hs.eg.db", "AnnotationHub",
               "readxl")
pkgs <- if (RUN_ENRICHMENT) c(bio_pkgs, core_pkgs) else core_pkgs
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  if (p %in% core_pkgs) install.packages(p, repos = "https://cloud.r-project.org")
  else {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}
# dplyr loaded last so select()/filter()/rename() aren't masked by AnnotationDbi
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

# Small helpers
pick_max <- function(x, v) { i <- which.max(v); if (length(i)) x[i] else x[NA_integer_] }
safemax <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
chr_num <- function(chr) suppressWarnings(as.numeric(ifelse(chr == "X", 27, chr)))

# ── 1b. READ GWAS SUMMARIES + WINDOWS_VARIANCE ───────────────────────────────
read_gwas <- function(path, trait) {
  read.csv(path, stringsAsFactors = FALSE) %>%
    filter(!is.na(CHR), !is.na(POS)) %>%
    mutate(CHR = as.character(CHR), Trait = trait)
}

read_windows <- function(path, trait, round) {
  w <- read.table(path, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("Trait_No", "Effect", "Start_SNP", "End_SNP",
                                "n_SNPs", "Start_ID", "End_ID", "Window_ID", "VE"))
  w <- w %>%
    mutate(CHR         = sub("_.*$", "", Start_ID),
           End_CHR     = sub("_.*$", "", End_ID),
           Win_Start   = as.numeric(sub("^.*_", "", Start_ID)),
           Win_End     = as.numeric(sub("^.*_", "", End_ID)),
           Win_Size_kb = (Win_End - Win_Start) / 1e3,
           Trait       = trait)
  # postGSf90 pads chromosome ends with empty windows (0 SNPs, start > end,
  # 0% variance, re-used Window_ID labels) — drop them
  n_empty <- sum(w$n_SNPs == 0)
  w <- w %>% filter(n_SNPs > 0)
  if (any(w$CHR != w$End_CHR))
    warning(sprintf("[%s] %d windows span two chromosomes", trait, sum(w$CHR != w$End_CHR)))
  if (anyDuplicated(w$Window_ID))
    warning(sprintf("[%s] duplicated Window_IDs among non-empty windows", trait))
  cat(sprintf("  %s round %d: %d windows (%d empty dropped), total VE = %.2f%%\n",
              trait, round, nrow(w), n_empty, sum(w$VE)))
  w %>% mutate(Round = round) %>% select(Trait, Round, Window_ID, CHR, Win_Start, Win_End, Win_Size_kb,
               n_SNPs, Start_SNP, End_SNP, VE)
}

cat("Reading inputs...\n")
gwas <- setNames(lapply(seq_len(nrow(TRAITS)), function(i)
  read_gwas(TRAITS$Summary[i], TRAITS$Trait[i])), TRAITS$Trait)
windows_all <- bind_rows(lapply(seq_len(nrow(TRAITS)), function(i)
  bind_rows(lapply(ROUNDS, function(r) {
    f <- file.path(TRAITS$Win_Dir[i], sprintf(WIN_PATTERN, r))
    if (!file.exists(f)) {
      if (r == CALL_ROUND) stop("windows_variance for the calling round not found: ", f)
      warning("Missing windows_variance (skipped): ", f); return(NULL)
    }
    read_windows(f, TRAITS$Trait[i], r)
  }))))
# QTL are called from the last round only
windows <- windows_all %>% filter(Round == CALL_ROUND) %>% select(-Round)

sig_windows <- windows %>% filter(VE >= VE_THRESHOLD)
cat("\nWindows passing threshold:\n"); print(table(sig_windows$Trait)); cat("\n")
if (nrow(sig_windows) == 0)
  stop(sprintf("No windows with VE >= %s%%. Lower VE_THRESHOLD.", VE_THRESHOLD))

# ── 1c. ALL SNP-LEVEL INFO INSIDE EACH SIG WINDOW ────────────────────────────
snps_in_windows <- bind_rows(lapply(seq_len(nrow(sig_windows)), function(i) {
  w <- sig_windows[i, ]
  gwas[[w$Trait]] %>%
    filter(CHR == w$CHR, POS >= w$Win_Start, POS <= w$Win_End) %>%
    mutate(Window_ID = w$Window_ID, Window_VE = w$VE) %>%
    select(Trait, Window_ID, Window_VE, SNPID, CHR, POS,
           snp_effect, var_Explained, neg_log_pval)
})) %>% arrange(Trait, chr_num(CHR), POS)

win_stats <- snps_in_windows %>%
  group_by(Trait, Window_ID) %>%
  summarise(n_SNPs_summary  = n(),
            Top_P_SNP       = pick_max(SNPID, neg_log_pval),
            Top_P_POS       = pick_max(POS,   neg_log_pval),
            Max_neg_log10p  = safemax(neg_log_pval),
            Top_VE_SNP      = pick_max(SNPID, var_Explained),
            Max_SNP_VE      = safemax(var_Explained),
            Top_Effect_SNP  = pick_max(SNPID, abs(snp_effect)),
            Top_Effect      = pick_max(snp_effect, abs(snp_effect)),
            SNPs            = paste(SNPID, collapse = "; "),
            .groups = "drop")

sig_windows <- sig_windows %>%
  left_join(win_stats, by = c("Trait", "Window_ID")) %>%
  arrange(Trait, desc(VE))

mm <- which(is.na(sig_windows$n_SNPs_summary) | sig_windows$n_SNPs_summary != sig_windows$n_SNPs)
if (length(mm))
  cat("NOTE: SNP count in windows_variance vs GWAS summary differs for:",
      paste(sig_windows$Trait[mm], sig_windows$Window_ID[mm], collapse = ", "), "\n\n")

cat("Significant windows:\n")
print(sig_windows %>% select(Trait, Window_ID, CHR, Win_Start, Win_End, n_SNPs, VE,
                             Top_P_SNP, Max_neg_log10p) %>% as.data.frame())
cat("\n")

# ── 1d. QTL TABLE: merge sig windows within CLUSTER_GAP_MB ───────────────────
qtl <- sig_windows %>%
  arrange(Trait, chr_num(CHR), Win_Start) %>%
  group_by(Trait, CHR) %>%
  mutate(gap = Win_Start - lag(Win_End),
         cluster = cumsum(is.na(gap) | gap > CLUSTER_GAP_MB * 1e6)) %>%
  group_by(Trait, CHR, cluster) %>%
  summarise(Win_Start      = min(Win_Start),
            Win_End        = max(Win_End),
            n_windows      = n(),
            Windows        = paste(Window_ID, collapse = "; "),
            VE_sum         = sum(VE),
            VE_max         = max(VE),
            Lead_Window    = pick_max(Window_ID, VE),
            n_SNPs         = sum(n_SNPs),
            Top_P_SNP      = pick_max(Top_P_SNP, Max_neg_log10p),
            Top_P_POS      = pick_max(Top_P_POS, Max_neg_log10p),
            Max_neg_log10p = safemax(Max_neg_log10p),
            .groups = "drop") %>%
  mutate(Start   = pmax(1, Win_Start - FLANK_KB * 1e3),
         End     = Win_End + FLANK_KB * 1e3,
         Size_kb = (End - Start) / 1e3,
         QTL_ID  = sprintf("%s_chr%s_%.2fMb", gsub(" ", "", Trait), CHR, Win_Start / 1e6)) %>%
  arrange(Trait, chr_num(CHR), Start) %>%
  select(QTL_ID, Trait, CHR, Start, End, Size_kb, Win_Start, Win_End, n_windows,
         Windows, Lead_Window, VE_max, VE_sum, n_SNPs, Top_P_SNP, Top_P_POS,
         Max_neg_log10p) %>%
  as.data.frame()

snps_in_windows <- snps_in_windows %>%
  left_join(qtl %>% select(QTL_ID, Trait, CHR, Win_Start, Win_End), by = c("Trait", "CHR"),
            relationship = "many-to-many") %>%
  filter(POS >= Win_Start, POS <= Win_End) %>%
  select(QTL_ID, everything(), -Win_Start, -Win_End)

# ── 1d-bis. VALUES OF EACH ROUND-CALLED QTL IN EVERY ROUND ─────────────────
# For each QTL (called in CALL_ROUND), sum / max the VE of that round's windows
# whose midpoint falls inside the QTL interval (robust if boundaries shift).
qtl_rounds <- bind_rows(lapply(seq_len(nrow(qtl)), function(i) {
  q <- qtl[i, ]
  windows_all %>%
    filter(Trait == q$Trait, CHR == q$CHR,
           (Win_Start + Win_End) / 2 >= q$Win_Start, (Win_Start + Win_End) / 2 <= q$Win_End) %>%
    group_by(Round) %>%
    summarise(VE_sum = sum(VE), VE_max = max(VE),
              Above = any(VE >= VE_THRESHOLD), .groups = "drop") %>%
    mutate(QTL_ID = q$QTL_ID)
}))
qtl_rounds_wide <- data.frame(QTL_ID = qtl$QTL_ID, stringsAsFactors = FALSE)
for (r in ROUNDS) {
  x <- qtl_rounds %>% filter(Round == r)
  qtl_rounds_wide[[sprintf("VE_sum_R%d", r)]] <- x$VE_sum[match(qtl$QTL_ID, x$QTL_ID)]
  qtl_rounds_wide[[sprintf("VE_max_R%d", r)]] <- x$VE_max[match(qtl$QTL_ID, x$QTL_ID)]
  qtl_rounds_wide[[sprintf("Above_R%d", r)]]  <- x$Above[match(qtl$QTL_ID, x$QTL_ID)]
}

cat("QTL table:\n"); print(qtl %>% select(QTL_ID, CHR, Start, End, Size_kb, n_windows,
                                          VE_max, VE_sum, Top_P_SNP)); cat("\n")

# ── 1e. WINDOW-VARIANCE PLOT WITH THRESHOLD ──────────────────────────────────
win_plot <- windows %>% mutate(mid = (Win_Start + Win_End) / 2, chr_n = chr_num(CHR))
chr_len <- win_plot %>% group_by(CHR, chr_n) %>% summarise(len = max(Win_End), .groups = "drop") %>%
  arrange(chr_n) %>% mutate(offset = cumsum(as.numeric(lag(len, default = 0))))
win_plot <- win_plot %>% left_join(chr_len, by = c("CHR", "chr_n")) %>%
  mutate(cum_pos = mid + offset, sig = VE >= VE_THRESHOLD)
axis_df <- win_plot %>% group_by(CHR, chr_n) %>%
  summarise(center = mean(range(cum_pos)), .groups = "drop") %>% arrange(chr_n)

p_ve <- ggplot(win_plot, aes(cum_pos, VE)) +
  geom_point(aes(colour = factor(chr_n %% 2)), size = 0.9, show.legend = FALSE) +
  geom_point(data = filter(win_plot, sig), colour = "firebrick", size = 1.8) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = filter(win_plot, sig), aes(label = Window_ID),
                               size = 2.8, colour = "firebrick", min.segment.length = 0)
    else geom_text(data = filter(win_plot, sig), aes(label = Window_ID),
                   vjust = -0.8, size = 2.8, colour = "firebrick") } +
  geom_hline(yintercept = VE_THRESHOLD, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("grey30", "steelblue")) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  facet_wrap(~Trait, ncol = 1) +
  labs(x = "Chromosome", y = "% genetic variance explained (window)",
       title = sprintf("postGSf90 window variance - threshold %s%%", VE_THRESHOLD)) +
  theme_bw(base_size = 11) + theme(panel.grid.major.x = element_blank())
ggsave(file.path(out_dir, "Window_VE_threshold.pdf"), p_ve,
       width = 11, height = 3 + 2.5 * nrow(TRAITS))

settings_tbl <- data.frame(
  Setting = c("VE_THRESHOLD_pct", "CALL_ROUND", "ROUNDS", "CLUSTER_GAP_MB", "FLANK_KB",
              "RUN_ENRICHMENT", paste0(TRAITS$Trait, "_summary"), paste0(TRAITS$Trait, "_win_dir")),
  Value   = c(VE_THRESHOLD, CALL_ROUND, paste(ROUNDS, collapse = ","), CLUSTER_GAP_MB, FLANK_KB,
              RUN_ENRICHMENT, TRAITS$Summary, TRAITS$Win_Dir))

out_sheets <- list(Settings        = settings_tbl,
                   QTL             = qtl,
                   Sig_Windows     = as.data.frame(sig_windows),
                   SNPs_in_QTL     = as.data.frame(snps_in_windows),
                   All_Windows     = windows_all)

if (RUN_ENRICHMENT) {

  # ── 0b. SHEEP OrgDb (AnnotationHub → fallback to org.Hs.eg.db) ─────────────
  sheep_orgdb <- tryCatch({
    ah <- AnnotationHub(); q <- query(ah, c("Ovis aries", "OrgDb"))
    if (length(q) > 0) { cat("Sheep OrgDb:", q$title[1], "\n"); q[[1]] } else org.Hs.eg.db
  }, error = function(e) { cat("AHub unavailable — using org.Hs.eg.db proxy\n"); org.Hs.eg.db })

  # ── 2. TxDb: LOAD CACHED OR BUILD FROM NCBI GFF ────────────────────────────
  GFF_URL   <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/298/735/",
                      "GCF_000298735.2_Oar_v4.0/GCF_000298735.2_Oar_v4.0_genomic.gff.gz")
  GFF_FILE  <- file.path(ANNOT_DIR, "Oar_v4.0.gff.gz")
  TXDB_FILE <- file.path(ANNOT_DIR, "Oar_v4.0_txdb.sqlite")

  if (file.exists(TXDB_FILE)) {
    cat("Loading cached TxDb:", TXDB_FILE, "\n")
    txdb <- loadDb(TXDB_FILE)
  } else {
    if (!file.exists(GFF_FILE)) { cat("Downloading GFF...\n"); download.file(GFF_URL, GFF_FILE, mode = "wb") }
    cat("Building TxDb (~5-10 min)...\n")
    txdb <- txdbmaker::makeTxDbFromGFF(GFF_FILE, format = "gff3")
    saveDb(txdb, TXDB_FILE); cat("TxDb saved.\n")
  }

  # ── 3. GENE GRanges + CHR REMAP ─────────────────────────────────────────────
  chr_map <- c(
    NC_019458.2="1",  NC_019459.2="2",  NC_019460.2="3",  NC_019461.2="4",
    NC_019462.2="5",  NC_019463.2="6",  NC_019464.2="7",  NC_019465.2="8",
    NC_019466.2="9",  NC_019467.2="10", NC_019468.2="11", NC_019469.2="12",
    NC_019470.2="13", NC_019471.2="14", NC_019472.2="15", NC_019473.2="16",
    NC_019474.2="17", NC_019475.2="18", NC_019476.2="19", NC_019477.2="20",
    NC_019478.2="21", NC_019479.2="22", NC_019480.2="23", NC_019481.2="24",
    NC_019482.2="25", NC_019483.2="26", NC_040951.1="X"
  )
  genes_gr <- genes(txdb)
  mapped   <- chr_map[seqlevels(genes_gr)]
  genes_gr <- keepSeqlevels(genes_gr, seqlevels(genes_gr)[!is.na(mapped)], pruning.mode = "coarse")
  seqlevels(genes_gr) <- mapped[!is.na(mapped)]
  genes_gr$gene_id <- sub("^gene-", "", genes_gr$gene_id)   # NCBI GFF IDs can carry "gene-"
  cat(sprintf("Genes on named chromosomes: %d\n\n", length(genes_gr)))

  # ── 4. GENE METADATA FROM GFF ───────────────────────────────────────────────
  cat("Extracting gene metadata from GFF...\n")
  gene_info <- read.table(GFF_FILE, sep = "\t", header = FALSE, quote = "",
                          comment.char = "#", stringsAsFactors = FALSE) %>%
    filter(V3 == "gene") %>%
    transmute(
      Gene_Name   = ifelse(grepl("Name=",        V9), sub(".*Name=([^;]+).*",         "\\1", V9), NA_character_),
      entrez_id   = ifelse(grepl("GeneID:",      V9), sub(".*GeneID:([0-9]+).*",      "\\1", V9), NA_character_),
      Description = ifelse(grepl("description=", V9), sub(".*description=([^;]+).*", "\\1", V9), NA_character_),
      Description = gsub("%2C", ",", gsub("%20", " ", Description))
    ) %>%
    filter(!is.na(Gene_Name), !is.na(entrez_id)) %>%
    distinct(Gene_Name, .keep_all = TRUE)
  cat(sprintf("Gene metadata: %d entries\n\n", nrow(gene_info)))

  # ── 5. QTL GENE QUERY ───────────────────────────────────────────────────────
  cat("Querying genes in QTL intervals...\n\n")
  results_list <- lapply(seq_len(nrow(qtl)), function(i) {
    r    <- qtl[i, ]
    hits <- findOverlaps(GRanges(as.character(r$CHR), IRanges(r$Start, r$End)),
                         genes_gr, ignore.strand = TRUE)
    base <- data.frame(QTL_ID = r$QTL_ID, Trait = r$Trait, CHR = r$CHR,
                       QTL_Start = r$Start, QTL_End = r$End, VE_max = r$VE_max,
                       Lead_Window = r$Lead_Window, Top_P_SNP = r$Top_P_SNP)
    if (length(hits) == 0) return(cbind(base, Note = "No genes found"))
    g <- as.data.frame(genes_gr[subjectHits(hits)]) %>%
      rename(Gene_Symbol = gene_id, Gene_Chr = seqnames,
             Gene_Start = start, Gene_End = end, Strand = strand) %>%
      left_join(gene_info, by = c("Gene_Symbol" = "Gene_Name")) %>%
      mutate(Dist_to_TopP_SNP_kb = ifelse(is.na(r$Top_P_POS), NA_real_,
               pmax(0, Gene_Start - r$Top_P_POS, r$Top_P_POS - Gene_End) / 1e3),
             Note = NA_character_)
    cbind(base[rep(1, nrow(g)), ], g) %>%
      select(QTL_ID, Trait, CHR, QTL_Start, QTL_End, VE_max, Lead_Window, Top_P_SNP,
             entrez_id, Gene_Symbol, Gene_Chr, Gene_Start, Gene_End, Strand,
             Dist_to_TopP_SNP_kb, Description, Note)
  })
  all_results <- bind_rows(results_list)
  rownames(all_results) <- NULL

  # ── 6. BACKGROUND FROM ARRAY SNP MAP ────────────────────────────────────────
  cat("Loading array SNP map:", ARRAY_FILE, "\n")
  snp_map <- read_xlsx(ARRAY_FILE) %>%
    filter(!is.na(Chr), !is.na(MapInfo)) %>%
    mutate(Chr = as.character(Chr), MapInfo = as.integer(MapInfo))
  cat(sprintf("Array SNPs: %d\n", nrow(snp_map)))

  snp_gr        <- GRanges(snp_map$Chr, IRanges(snp_map$MapInfo, snp_map$MapInfo))
  snp_hits      <- findOverlaps(snp_gr, genes_gr, ignore.strand = TRUE)
  bg_symbols    <- unique(genes_gr$gene_id[subjectHits(snp_hits)])
  bg_entrez_ids <- gene_info %>% filter(Gene_Name %in% bg_symbols) %>% pull(entrez_id) %>% unique()
  cat(sprintf("Background genes (SNP-covered, with Entrez ID): %d\n\n", length(bg_entrez_ids)))

  # ── 7. ORA: KEGG + GO (BP / MF / CC) ────────────────────────────────────────
  fg_genes <- all_results %>%
    filter(!is.na(entrez_id)) %>%
    group_by(Trait) %>%
    summarise(gene_ids = list(unique(entrez_id)), .groups = "drop")

  ora_results <- list()
  for (i in seq_len(nrow(fg_genes))) {
    trait <- fg_genes$Trait[i]
    fg    <- fg_genes$gene_ids[[i]]
    cat(sprintf("\n── ORA: %s (%d genes) ──\n", trait, length(fg)))

    kegg_res <- tryCatch(
      enrichKEGG(gene = fg, organism = "oas", universe = bg_entrez_ids,
                 pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
      error = function(e) { cat("  KEGG error:", conditionMessage(e), "\n"); NULL })
    if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
      ora_results[[paste0(trait, "_KEGG")]] <- as.data.frame(kegg_res) %>% mutate(Trait = trait, DB = "KEGG")
      cat(sprintf("  KEGG: %d hits\n", nrow(ora_results[[paste0(trait, "_KEGG")]])))
      ggsave(file.path(out_dir, sprintf("ORA_KEGG_%s.pdf", gsub(" ", "_", trait))),
             dotplot(kegg_res, title = paste(trait, "— KEGG")) + theme_bw(base_size = 11),
             width = 8, height = 6)
    } else cat("  KEGG: no hits\n")

    for (ont in c("BP", "MF", "CC")) {
      go_res <- tryCatch(
        enrichGO(gene = fg, OrgDb = sheep_orgdb, keyType = "ENTREZID", ont = ont,
                 universe = bg_entrez_ids, pAdjustMethod = "BH",
                 pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
        error = function(e) { cat(sprintf("  GO %s error: %s\n", ont, conditionMessage(e))); NULL })
      if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
        go_simp <- tryCatch(simplify(go_res, cutoff = 0.7, by = "p.adjust", select_fun = min),
                            error = function(e) go_res)
        ora_results[[paste0(trait, "_GO_", ont)]] <- as.data.frame(go_simp) %>%
          mutate(Trait = trait, DB = paste0("GO_", ont))
        cat(sprintf("  GO %s: %d hits (post-simplify)\n", ont, nrow(as.data.frame(go_simp))))
        ggsave(file.path(out_dir, sprintf("ORA_GO%s_%s.pdf", ont, gsub(" ", "_", trait))),
               dotplot(go_simp, showCategory = 20,
                       title = paste(trait, "—", ont, "GO")) + theme_bw(base_size = 11),
               width = 8, height = 7)
      } else cat(sprintf("  GO %s: no hits\n", ont))
    }
  }

  # ── 8. SUMMARY TABLES ───────────────────────────────────────────────────────
  summary_tbl <- qtl %>%
    left_join(all_results %>%
                filter(!is.na(Gene_Symbol)) %>%
                group_by(QTL_ID) %>%
                summarise(n_genes    = n_distinct(Gene_Symbol),
                          gene_names = paste(sort(unique(Gene_Symbol)), collapse = "; "),
                          .groups = "drop"),
              by = "QTL_ID") %>%
    mutate(n_genes = ifelse(is.na(n_genes), 0L, n_genes))

  shared_genes <- all_results %>%
    filter(!is.na(entrez_id)) %>%
    group_by(entrez_id, Gene_Symbol, CHR) %>%
    summarise(Traits = paste(sort(unique(Trait)), collapse = " & "),
              n_traits = n_distinct(Trait), .groups = "drop") %>%
    filter(n_traits > 1)

  ora_combined <- bind_rows(ora_results)
  split_ora <- function(db) {
    if (nrow(ora_combined) == 0) return(data.frame(Note = "No hits"))
    x <- filter(ora_combined, DB == db); if (nrow(x) > 0) x else data.frame(Note = "No hits")
  }

  cat("\n=== QTL SUMMARY ===\n"); print(summary_tbl %>% select(QTL_ID, VE_max, n_genes, gene_names))
  cat("\n=== SHARED GENES ===\n"); print(shared_genes)
  if (nrow(ora_combined) > 0) {
    cat("\n=== ORA ===\n")
    print(ora_combined %>% select(Trait, DB, ID, Description, pvalue, p.adjust, geneID))
  }

  out_sheets <- c(out_sheets[c("Settings")],
                  list(QTL_Summary = summary_tbl),
                  out_sheets[c("Sig_Windows", "SNPs_in_QTL")],
                  list(All_Genes    = all_results,
                       Shared_Genes = if (nrow(shared_genes)) shared_genes else data.frame(Note = "None"),
                       Background   = data.frame(entrez_id = bg_entrez_ids),
                       ORA_KEGG     = split_ora("KEGG"),
                       ORA_GO_BP    = split_ora("GO_BP"),
                       ORA_GO_MF    = split_ora("GO_MF"),
                       ORA_GO_CC    = split_ora("GO_CC")),
                  out_sheets["All_Windows"])
  write.table(all_results, file.path(out_dir, "QTL_genes_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ── 9. QTL TABLE (round 3 windows) ───────────────────────────────────────────
# One row per QTL called from the round-3 windows_variance, with the top SNPs
# from the lead window and (if enrichment ran) the genes in the interval,
# ordered by distance to the top SNP.
gene_by_qtl <- if (exists("all_results") && "Gene_Symbol" %in% names(all_results)) {
  all_results %>%
    filter(!is.na(Gene_Symbol)) %>%
    arrange(QTL_ID, Dist_to_TopP_SNP_kb) %>%
    group_by(QTL_ID) %>%
    summarise(n_genes      = n_distinct(Gene_Symbol),
              Nearest_gene = { g <- Gene_Symbol[!grepl("^LOC", Gene_Symbol)]
                               if (length(g)) g[1] else Gene_Symbol[1] },
              Genes        = paste(unique(Gene_Symbol), collapse = ", "),
              .groups = "drop")
} else NULL

qtl_table <- qtl %>%
  left_join(sig_windows %>% select(Trait, Window_ID, Top_VE_SNP, Max_SNP_VE,
                                   Top_Effect_SNP, Top_Effect),
            by = c("Trait", "Lead_Window" = "Window_ID")) %>%
  left_join(qtl_rounds_wide, by = "QTL_ID")
if (!is.null(gene_by_qtl)) {
  qtl_table <- qtl_table %>% left_join(gene_by_qtl, by = "QTL_ID") %>%
    mutate(n_genes = ifelse(is.na(n_genes), 0L, n_genes))
} else {
  qtl_table <- qtl_table %>% mutate(n_genes = NA_integer_, Nearest_gene = NA_character_,
                                    Genes = NA_character_)
}
qtl_table <- qtl_table %>% arrange(Trait, chr_num(CHR), Start)

# Per-round variance columns (R1, R2, ... R<CALL_ROUND>)
round_cols <- do.call(cbind, lapply(ROUNDS, function(r) {
  x <- data.frame(round(qtl_table[[sprintf("VE_sum_R%d", r)]], 2),
                  round(qtl_table[[sprintf("VE_max_R%d", r)]], 2),
                  ifelse(qtl_table[[sprintf("Above_R%d", r)]], "Yes", "No"))
  names(x) <- c(sprintf("Var. explained R%d (%%)", r),
                sprintf("Lead window var. R%d (%%)", r),
                sprintf("Window >= %s%% in R%d", VE_THRESHOLD, r))
  x
}))

qtl_table <- bind_cols(
  qtl_table %>% transmute(Trait,
            QTL                         = QTL_ID,
            Chr                         = CHR,
            `Start (bp)`                = Start,
            `End (bp)`                  = End,
            `Size (kb)`                 = round(Size_kb, 1),
            `Windows`                   = Windows,
            `No. windows`               = n_windows,
            `No. SNPs`                  = n_SNPs),
  round_cols,
  qtl_table %>% transmute(
            `Top SNP (P)`               = Top_P_SNP,
            `Top SNP position (bp)`     = Top_P_POS,
            `-log10(P)`                 = round(Max_neg_log10p, 2),
            `Top SNP (variance)`        = Top_VE_SNP,
            `Top SNP var. (%)`          = signif(Max_SNP_VE, 3),
            `Top SNP (effect)`          = Top_Effect_SNP,
            `SNP effect`                = signif(Top_Effect, 3),
            `No. genes`                 = n_genes,
            `Nearest gene to top SNP`   = Nearest_gene,
            `Genes (closest first)`     = Genes))

cat(sprintf("\n=== QTL TABLE (called in round %d; values in rounds %s) ===\n",
            CALL_ROUND, paste(ROUNDS, collapse = ", ")))
print(qtl_table[, c("Trait", "QTL", grep("^Var. explained R", names(qtl_table), value = TRUE),
                    "Top SNP (P)", "Nearest gene to top SNP")], row.names = FALSE)
write.csv(qtl_table, file.path(out_dir, sprintf("Table_QTL_round%d_VE%s.csv", CALL_ROUND, VE_THRESHOLD)),
          row.names = FALSE)
out_sheets <- c(setNames(list(qtl_table), sprintf("Table_QTL_round%d", CALL_ROUND)), out_sheets)

# ── 10. SUPPLEMENTARY FIGURES: SNP var-explained Manhattan, rounds as rows ─
# Same data as BLUPF90.GWAS.plots.R manhattan_var_explained: per-SNP
# var_Explained from round_<k>/chrsnpvar (variance of the window starting at
# each SNP). One figure per trait (Fig S1 = first trait in TRAITS, ...),
# row k = round k. With FIG_ANNOTATE = TRUE, adds the VE_THRESHOLD line,
# light grey shading over the final QTL (called in the last round) in EVERY
# round (padded 2 Mb each side so ~1-Mb QTL stay visible), and the nearest gene (or Chr:Mb) at each QTL's peak in the last-round
# row. FIG_ANNOTATE = FALSE gives rows identical to manhattan_var_explained.
read_snpvar <- function(path, trait, round) {
  x <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  data.frame(Trait = trait, Round = round, CHR = as.character(x$V5),
             POS = as.numeric(x$V6), VE = as.numeric(x$V3), stringsAsFactors = FALSE)
}
snpvar_all <- bind_rows(lapply(seq_len(nrow(TRAITS)), function(i)
  bind_rows(lapply(ROUNDS, function(r) {
    f <- file.path(TRAITS$Win_Dir[i], sprintf(SNPVAR_PATTERN, r))
    if (file.exists(f)) return(read_snpvar(f, TRAITS$Trait[i], r))
    if (r == CALL_ROUND) {   # fall back to the summary CSV (= round-3 chrsnpvar)
      cat(sprintf("  %s round %d: %s not found, using GWAS summary var_Explained\n",
                  TRAITS$Trait[i], r, f))
      g <- gwas[[TRAITS$Trait[i]]]
      return(data.frame(Trait = TRAITS$Trait[i], Round = r, CHR = g$CHR, POS = g$POS,
                        VE = g$var_Explained, stringsAsFactors = FALSE))
    }
    warning("Missing chrsnpvar (row skipped): ", f); NULL
  })))) %>% filter(!is.na(VE), !is.na(POS))

# Cumulative x-axis exactly as in BLUPF90.GWAS.plots.R
fig_chr <- snpvar_all %>% mutate(chr_n = chr_num(CHR)) %>%
  group_by(CHR, chr_n) %>% summarise(chr_len = max(POS), .groups = "drop") %>%
  arrange(chr_n) %>% mutate(cum_start = cumsum(as.numeric(chr_len)) - chr_len)
snp_plot_all <- snpvar_all %>%
  left_join(fig_chr %>% select(CHR, cum_start), by = "CHR") %>%
  mutate(bp_cum = POS + cum_start,
         sig = VE >= VE_THRESHOLD,
         CHR = factor(CHR, levels = fig_chr$CHR),
         Round_lab = factor(paste("Round", Round), levels = paste("Round", ROUNDS)))
# Flag SNPs inside the final (last-round) QTL so they are highlighted in every round
snp_plot_all$in_qtl <- FALSE
for (i in seq_len(nrow(qtl))) {
  q <- qtl[i, ]
  hit <- snp_plot_all$Trait == q$Trait & as.character(snp_plot_all$CHR) == q$CHR &
         snp_plot_all$POS >= q$Win_Start & snp_plot_all$POS <= q$Win_End
  snp_plot_all$in_qtl[hit] <- TRUE
}
fig_axis <- snp_plot_all %>% group_by(CHR) %>%
  summarise(center = (min(bp_cum) + max(bp_cum)) / 2, .groups = "drop")

qtl_bands <- qtl %>% left_join(fig_chr %>% select(CHR, cum_start), by = "CHR") %>%
  mutate(xmin = Win_Start + cum_start, xmax = Win_End + cum_start)

# Label at the highest last-round SNP inside each QTL
lab_df <- bind_rows(lapply(seq_len(nrow(qtl)), function(i) {
  q <- qtl[i, ]
  x <- snp_plot_all %>% filter(Trait == q$Trait, Round == CALL_ROUND,
                               as.character(CHR) == q$CHR,
                               POS >= q$Win_Start, POS <= q$Win_End)
  if (nrow(x) == 0) return(NULL)
  x[which.max(x$VE), ] %>% mutate(QTL_ID = q$QTL_ID)
})) 
if (nrow(lab_df)) lab_df <- lab_df %>%
  left_join(qtl_table %>% select(QTL, Nearest = `Nearest gene to top SNP`),
            by = c("QTL_ID" = "QTL")) %>%
  mutate(label = ifelse(is.na(Nearest), sprintf("Chr%s:%.1f Mb", CHR, POS / 1e6), Nearest))

# Merge labels of QTL on the same chromosome within LABEL_MERGE_MB into one
# label at the highest peak (avoids a fan of leader lines, e.g. chr16)
LABEL_MERGE_MB <- 15
if (nrow(lab_df)) lab_df <- lab_df %>%
  mutate(CHR_chr = as.character(CHR)) %>%
  arrange(Trait, chr_num(CHR_chr), POS) %>%
  group_by(Trait, CHR_chr) %>%
  mutate(grp = cumsum(c(TRUE, diff(POS) > LABEL_MERGE_MB * 1e6))) %>%
  group_by(Trait, CHR_chr, grp) %>%
  summarise(label = if (n() == 1) first(label)
                    else if (all(is.na(Nearest)))
                      sprintf("Chr%s: %s Mb", first(CHR_chr),
                              paste(sprintf("%.1f", POS / 1e6), collapse = ", "))
                    else paste(ifelse(is.na(Nearest), sprintf("%.1f Mb", POS / 1e6), Nearest),
                               collapse = ", "),
            bp_cum = bp_cum[which.max(VE)], VE = max(VE),
            Round = first(Round), Round_lab = first(Round_lab), CHR = first(CHR),
            .groups = "drop")

# Style matched to BLUPF90.GWAS.plots.R (manhattan_var_explained):
# trait-specific palettes cycled across chromosomes, accent colour for
# windows above the threshold, same theme / fonts / axis weights.
trait_palette <- function(trait_label) {
  tl <- tolower(trait_label)
  is_avg <- grepl("(^|[^a-z])(ave|avg|aver|average|mean)([^a-z]|$)", tl)
  is_max <- grepl("(^|[^a-z])(max|maximum|peak)([^a-z]|$)", tl)
  if (is_max) {
    list(chr = c("#6A040F", "#E85D04", "#FFBA08"), hi = "#0072B2", line = "#0072B2",
         name = "maximum temperature (maroon / vermillion / amber)")
  } else if (is_avg) {
    list(chr = c("#000080", "#00CED1", "#FF7F50"), hi = "#E6001A", line = "#E6001A",
         name = "average temperature (navy / teal / coral)")
  } else {
    list(chr = c("#1B262C", "#3B9EFF", "#FFBF00"), hi = "#E6001A", line = "#E6001A",
         name = "generic (charcoal / bright blue / amber)")
  }
}

TIFF_DPI  <- suppressWarnings(as.numeric(Sys.getenv("TIFF_DPI", unset = "600")))
if (is.na(TIFF_DPI) || TIFF_DPI <= 0) TIFF_DPI <- 600
# Sizes scaled for a 10 x 9 in, 3-row figure (the single-plot script uses
# 15-pt text / 1.1 lines on 13 x 5.2 in, which crowds a stacked layout)
base_font <- 11
axis_lwd  <- 0.6
tick_lwd  <- 0.6

gwas_theme <- theme_minimal(base_size = base_font) +
  theme(
    plot.title         = element_blank(),
    plot.subtitle      = element_blank(),
    plot.caption       = element_blank(),
    axis.title         = element_text(size = base_font + 1, face = "bold"),
    axis.title.x       = element_text(margin = margin(t = 6)),
    axis.title.y       = element_text(margin = margin(r = 6)),
    axis.text.x        = element_text(size = base_font - 2.5, color = "black"),
    axis.text.y        = element_text(size = base_font - 2, color = "black"),
    axis.line          = element_line(color = "black", linewidth = axis_lwd, lineend = "square"),
    axis.ticks         = element_line(color = "black", linewidth = tick_lwd),
    axis.ticks.length  = unit(3, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.border       = element_blank(),
    legend.position    = "none",
    strip.text         = element_text(size = base_font, face = "bold", hjust = 0,
                                      margin = margin(b = 2, t = 2)),
    strip.background   = element_blank(),
    panel.spacing.y    = unit(0.9, "lines"),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(8, 12, 8, 8)
  )

# Keep a tick for every chromosome but only label those with room
# (small chromosomes at the right end otherwise collide)
thin_labels <- function(ax, min_frac = 0.022) {
  gap <- min_frac * diff(range(ax$center)); last <- -Inf
  lab <- as.character(ax$CHR)
  for (i in seq_len(nrow(ax))) {
    if (ax$center[i] - last >= gap) last <- ax$center[i] else lab[i] <- ""
  }
  lab
}

save_plot <- function(plot, stem, width, height) {
  ok_pdf <- tryCatch({ ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
                              device = cairo_pdf); TRUE }, error = function(e) FALSE)
  if (!ok_pdf) ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ok <- tryCatch({ ggsave(paste0(stem, ".tiff"), plot, width = width, height = height,
                          device = "tiff", dpi = TIFF_DPI, compression = "lzw",
                          type = "cairo", bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok) ok <- tryCatch({ ggsave(paste0(stem, ".tiff"), plot, width = width, height = height,
                                   device = "tiff", dpi = TIFF_DPI, compression = "lzw",
                                   bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok) warning("Could not write TIFF: ", stem, ".tiff")
}

make_supp_fig <- function(tr) {
  pal <- trait_palette(tr)
  d   <- filter(snp_plot_all, Trait == tr)
  ld  <- if (nrow(lab_df)) filter(lab_df, Trait == tr) else lab_df
  qb  <- filter(qtl_bands, Trait == tr)
  chr_colors <- setNames(pal$chr[((seq_along(fig_chr$CHR) - 1) %% length(pal$chr)) + 1],
                         fig_chr$CHR)
  cat(sprintf("  %s palette: %s\n", tr, pal$name))

  p <- ggplot(d, aes(x = bp_cum, y = VE, color = CHR))
  if (FIG_ANNOTATE && nrow(qb))
    p <- p + geom_rect(data = qb, inherit.aes = FALSE,
                       aes(xmin = xmin - 2e6, xmax = xmax + 2e6, ymin = -Inf, ymax = Inf),
                       fill = "grey85", alpha = 0.6)
  p <- p + geom_point(alpha = 0.8, size = 0.7, shape = 16)
  if (FIG_ANNOTATE) {
    p <- p + geom_hline(yintercept = VE_THRESHOLD, linetype = "dashed",
                        color = pal$line, linewidth = 0.5)
    if (nrow(ld)) {
      # labels pushed up into the headroom of the last-round panel
      top <- max(filter(d, Round == CALL_ROUND)$VE, VE_THRESHOLD, na.rm = TRUE)
      p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_text_repel(data = ld, aes(label = label), size = 2.8,
                                 fontface = "italic", color = "black",
                                 ylim = c(top * 1.05, NA), direction = "x",
                                 segment.size = 0.3, segment.color = "grey40",
                                 min.segment.length = 0, box.padding = 0.25,
                                 max.overlaps = Inf, seed = 1)
      else geom_text(data = ld, aes(label = label), vjust = -0.9, size = 2.8,
                     fontface = "italic", color = "black")
    }
  }
  p + scale_color_manual(values = chr_colors) +
    scale_x_continuous(breaks = fig_axis$center, labels = thin_labels(fig_axis),
                       expand = expansion(mult = 0.005)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                       labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE),
                       expand = expansion(mult = c(0.02, if (FIG_ANNOTATE) 0.22 else 0.06))) +
    facet_wrap(~Round_lab, ncol = 1, scales = "free_y", strip.position = "top") +
    labs(x = "Chromosome", y = "% of Additive variance explained") +
    gwas_theme
}

fig_w <- 10
fig_h <- 1.2 + 2.6 * length(ROUNDS)  # 10 x 9 in for 3 rounds
for (k in seq_len(nrow(TRAITS))) {
  tr   <- TRAITS$Trait[k]
  stem <- file.path(out_dir, sprintf("FigS%d_WindowVE_Manhattan_rounds_%s", k, gsub(" ", "_", tr)))
  save_plot(make_supp_fig(tr), stem, fig_w, fig_h)
  cat(sprintf("Figure S%d (%s): %s.[pdf|png|tiff] (TIFF %g dpi)\n", k, tr, stem, TIFF_DPI))
}

# ── 11. EXPORT ───────────────────────────────────────────────────────────────
write_xlsx(out_sheets,
           path = file.path(out_dir, sprintf("QTL_windows_OarV4_NCBI_results_VE%s.xlsx", VE_THRESHOLD)))
write.csv(sig_windows, file.path(out_dir, sprintf("Sig_windows_VE%s.csv", VE_THRESHOLD)), row.names = FALSE)
cat(sprintf("\nAll outputs saved to: %s\nDone.\n", out_dir))