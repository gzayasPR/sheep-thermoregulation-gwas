#!/usr/bin/env Rscript
# =============================================================================
# popgen.plots.R
#
# Plots from popgen pipeline outputs. Reads env vars exported by
# run_popgen.pipe.sh: PREP_DIR, ADMIX_DIR, PCA_DIR, PLOT_DIR, POP_PREFIX,
# PHENO, model_csv, GROUP_COL, GROUP_COL2 (optional), TIFF_DPI (optional,
# default 600)
#
# GROUP_COL / GROUP_COL2 can each be a Name from the model csv (e.g. "BREED")
# or a header name in the phenotype file (e.g. "FARM").
#   GROUP_COL  -> PCA colour; outer split in ADMIXTURE bar plots
#   GROUP_COL2 -> PCA shape;  inner split in ADMIXTURE bar plots
#
# Output (each as .pdf + LZW .tiff) in PLOT_DIR. <tag> is built from the
# grouping labels, e.g. "_by_BREED", "_by_FARM", "_by_BREED_FARM":
#   pca_PC1_PC2<tag>, pca_PC3_PC4<tag>, pca_scree
#   admixture_CV_error
#   admixture_barplot_all_K<tag>   (all K stacked)
#   admixture_barplot_bestK<tag>   (K with lowest CV error)
#   With both labels set, bar plots are written for GROUP_COL alone,
#   GROUP_COL2 alone, and GROUP_COL x GROUP_COL2 (nested).
#   + admixture_Q_K<best>.csv       (IDs, group labels, ancestry proportions)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

env <- function(x, default = "") Sys.getenv(x, unset = default)
prep_dir   <- env("PREP_DIR");  admix_dir <- env("ADMIX_DIR")
pca_dir    <- env("PCA_DIR");   plot_dir  <- env("PLOT_DIR", ".")
prefix     <- env("POP_PREFIX", "popgen")
pheno_file <- env("PHENO");     model_csv <- env("model_csv")
group_col  <- env("GROUP_COL"); group_col2 <- env("GROUP_COL2")
tiff_dpi   <- suppressWarnings(as.numeric(env("TIFF_DPI", "600")))
if (is.na(tiff_dpi) || tiff_dpi <= 0) tiff_dpi <- 600
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Theme + save helper (same look as BLUPF90.GWAS.plots.R)
# ---------------------------------------------------------------------------
base_font <- 15
pop_theme <- theme_minimal(base_size = base_font) +
  theme(
    axis.title       = element_text(size = base_font + 2, face = "bold"),
    axis.text        = element_text(color = "black"),
    axis.line        = element_line(color = "black", linewidth = 1.1),
    axis.ticks       = element_line(color = "black", linewidth = 1.1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.title     = element_text(face = "bold"),
    plot.margin      = margin(12, 18, 12, 14)
  )

save_plot <- function(plot, basename, width, height) {
  ggsave(file.path(plot_dir, paste0(basename, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  tf <- file.path(plot_dir, paste0(basename, ".tiff"))
  ok <- tryCatch({
    ggsave(tf, plot, width = width, height = height, device = "tiff",
           dpi = tiff_dpi, compression = "lzw", type = "cairo", bg = "white"); TRUE
  }, error = function(e) tryCatch({
    ggsave(tf, plot, width = width, height = height, device = "tiff",
           dpi = tiff_dpi, compression = "lzw", bg = "white"); TRUE
  }, error = function(e2) FALSE))
  if (!ok) warning("Could not write TIFF: ", tf)
}

# Colorblind-friendly qualitative palette (Okabe-Ito + extras), up to 12
qual_pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
              "#D55E00", "#F0E442", "#000000", "#999999", "#882255",
              "#44AA99", "#117733")
# Distinguishable point shapes (filled first, then open, then line-based)
shape_pal <- c(16, 17, 15, 18, 1, 2, 0, 5, 6, 3, 4, 8, 7, 9, 10, 12, 13, 14)

safe_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ---------------------------------------------------------------------------
# 1. Sample IDs + optional group labels from the phenotype file
# ---------------------------------------------------------------------------
fam <- fread(file.path(prep_dir, paste0(prefix, ".pruned.fam")), header = FALSE)
ids <- data.table(FID = as.character(fam$V1), IID = as.character(fam$V2))

pheno <- meta <- NULL
if (file.exists(pheno_file) && file.exists(model_csv)) {
  pheno <- fread(pheno_file)
  meta  <- fread(model_csv)
  setnames(meta, 1, "ColN")                         # guard against BOM in header
}

# Returns a character vector of labels aligned with ids$IID, or NULL
get_group <- function(col) {
  if (!nzchar(col)) return(NULL)
  if (is.null(pheno)) {
    warning("PHENO / model_csv not found; ignoring grouping '", col, "'")
    return(NULL)
  }
  id_col  <- meta$ColN[meta$What == "GenoID"]
  # col may be a Name from the model csv (e.g. "BREED" -> ColN 5)
  # or a header name in the phenotype file
  grp_idx <- meta$ColN[meta$Name == col]
  grp_idx <- if (length(grp_idx) == 1) grp_idx else match(col, names(pheno))
  if (length(id_col) != 1 || is.na(grp_idx)) {
    warning("Group column '", col, "' or GenoID column not found; ignoring it")
    return(NULL)
  }
  lut <- unique(data.table(IID   = as.character(pheno[[id_col]]),
                           Group = as.character(pheno[[grp_idx]])), by = "IID")
  g <- lut$Group[match(ids$IID, lut$IID)]
  cat(sprintf("Grouping by '%s': %d of %d animals labeled, %d levels\n",
              col, sum(!is.na(g)), nrow(ids), length(unique(na.omit(g)))))
  g[is.na(g) | g == ""] <- "Unknown"
  g
}

g1 <- get_group(group_col)
g2 <- get_group(group_col2)
# If only GROUP_COL2 worked, promote it to the primary grouping
if (is.null(g1) && !is.null(g2)) { g1 <- g2; group_col <- group_col2; g2 <- NULL }
ids$G1 <- if (!is.null(g1)) g1 else NA_character_
ids$G2 <- if (!is.null(g2)) g2 else NA_character_
has_g1 <- !is.null(g1); has_g2 <- !is.null(g2)

# ---------------------------------------------------------------------------
# 2. ADMIXTURE: CV error + Q matrices
# ---------------------------------------------------------------------------
cv_file <- file.path(admix_dir, "CV_error.txt")
cv <- fread(cv_file)
best_k <- cv$K[which.min(cv$CV_error)]
cat(sprintf("Best K by CV error: %d\n", best_k))

p_cv <- ggplot(cv, aes(K, CV_error)) +
  geom_line(color = "#0072B2", linewidth = 1) +
  geom_point(color = "#0072B2", size = 3) +
  geom_point(data = cv[K == best_k], color = "#D55E00", size = 4.5) +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K", y = "Cross-validation error") +
  pop_theme
save_plot(p_cv, "admixture_CV_error", 7, 5)

read_q <- function(k) {
  q <- fread(file.path(admix_dir, sprintf("%s.pruned.%d.Q", prefix, k)), header = FALSE)
  setnames(q, paste0("Q", seq_len(k)))
  cbind(ids, q)
}

# Dominant cluster at best K, used to sort animals within each group
q_best <- read_q(best_k)
qcols  <- paste0("Q", seq_len(best_k))
q_best$dom      <- max.col(as.matrix(q_best[, ..qcols]), ties.method = "first")
q_best$dom_prop <- apply(as.matrix(q_best[, ..qcols]), 1, max)

# Q table for best K, with group columns under their real names
q_out <- q_best[order(G1, G2, dom, -dom_prop)]
q_out[, dom_prop := NULL]
if (has_g1) setnames(q_out, "G1", group_col) else q_out[, G1 := NULL]
if (has_g2) setnames(q_out, "G2", group_col2) else q_out[, G2 := NULL]
fwrite(q_out, file.path(plot_dir, sprintf("admixture_Q_K%d.csv", best_k)))

all_q <- rbindlist(lapply(cv$K, function(k) {
  m <- melt(read_q(k), id.vars = c("FID", "IID", "G1", "G2"),
            variable.name = "Cluster", value.name = "Proportion")
  m$K <- factor(paste0("K = ", k), levels = paste0("K = ", cv$K))
  m
}))

# x positions for every animal: sorted by outer group, inner group, then
# dominant cluster; a wide gap between outer groups, a narrow one between
# inner groups. Returns per-animal positions plus label positions.
bar_layout <- function(outer, inner) {
  d <- data.table(IID = q_best$IID, outer = outer, inner = inner,
                  dom = q_best$dom, dom_prop = q_best$dom_prop)
  if (is.null(inner)) d$inner <- ""
  if (is.null(outer)) d$outer <- ""
  d <- d[order(outer, inner, dom, -dom_prop)]
  n <- nrow(d)
  gap_big   <- max(1.5, 0.015 * n)
  gap_small <- max(0.6, 0.005 * n)
  step <- c(0, ifelse(d$outer[-1] != d$outer[-n], 1 + gap_big,
                ifelse(d$inner[-1] != d$inner[-n], 1 + gap_small, 1)))
  d$pos <- cumsum(step)
  list(animals = d,
       outer   = d[, .(mid = (min(pos) + max(pos)) / 2), by = outer],
       inner   = d[, .(mid = (min(pos) + max(pos)) / 2), by = .(outer, inner)])
}

# outer / inner: label vectors aligned with q_best (NULL = not used)
admix_bar <- function(df, outer = NULL, inner = NULL) {
  lay <- bar_layout(outer, inner)
  df  <- merge(df, lay$animals[, .(IID, pos)], by = "IID")
  pal <- rep(qual_pal, length.out = length(unique(df$Cluster)))

  # Bottom axis: inner labels if nested, else outer labels.
  # Top axis (nested only): outer labels over each block of inner groups.
  if (!is.null(inner)) {
    bottom <- lay$inner; bot_lab <- bottom$inner
    top_axis <- sec_axis(~ ., breaks = lay$outer$mid, labels = lay$outer$outer)
  } else if (!is.null(outer)) {
    bottom <- lay$outer; bot_lab <- bottom$outer
    top_axis <- waiver()
  } else {
    bottom <- NULL
  }

  p <- ggplot(df, aes(pos, Proportion, fill = Cluster)) +
    geom_col(width = 1, color = NA) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    facet_grid(K ~ .) +
    labs(x = NULL, y = "Ancestry proportion") +
    pop_theme +
    theme(axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_text(size = base_font - 5),
          legend.position = "none",
          panel.spacing.y = unit(0.8, "lines"),
          strip.text = element_text(face = "bold"),
          plot.margin = margin(12, 18, 12, 40))

  if (is.null(bottom)) {
    return(p + scale_x_continuous(expand = c(0.005, 0)) +
             theme(axis.text.x = element_blank()))
  }
  # Labels tilted 45 deg; each ends at the centre of its group
  # Nested: extra right margin so the last tilted top label isn't clipped
  r_marg <- if (!is.null(inner)) 110 else 18
  p + scale_x_continuous(expand = c(0.005, 0), breaks = bottom$mid,
                         labels = bot_lab, sec.axis = top_axis) +
    theme(plot.margin = margin(12, r_marg, 12, 40),
          axis.text.x.bottom = element_text(size = base_font - if (is.null(inner)) 5 else 6,
                                            face = if (is.null(inner)) "bold" else "plain",
                                            angle = 45, hjust = 1, vjust = 1),
          axis.text.x.top    = element_text(size = base_font - 4, face = "bold",
                                            angle = 45, hjust = 0, vjust = 0))
}

# Which groupings to draw: each label alone, plus nested when both are set
bar_sets <- list()
if (!has_g1) bar_sets[["none"]] <- list(o = NULL, i = NULL, tag = "")
if (has_g1)  bar_sets[["g1"]] <- list(o = ids$G1, i = NULL,
                                      tag = paste0("_by_", safe_tag(group_col)))
if (has_g2)  bar_sets[["g2"]] <- list(o = ids$G2, i = NULL,
                                      tag = paste0("_by_", safe_tag(group_col2)))
if (has_g2)  bar_sets[["g12"]] <- list(o = ids$G1, i = ids$G2,
                                       tag = paste0("_by_", safe_tag(group_col), "_",
                                                    safe_tag(group_col2)))
# ids rows and q_best rows share the .fam order, so these vectors line up

for (s in bar_sets) {
  extra_h <- if (!is.null(s$i)) 1.5 else 0      # room for the top label row
  save_plot(admix_bar(all_q, s$o, s$i),
            paste0("admixture_barplot_all_K", s$tag),
            width = 14, height = 1.4 * length(cv$K) + 1.5 + extra_h)
  save_plot(admix_bar(all_q[K == paste0("K = ", best_k)], s$o, s$i) +
              theme(legend.position = "right"),
            paste0("admixture_barplot_bestK", s$tag),
            width = 14, height = 4 + extra_h)
}

# ---------------------------------------------------------------------------
# 3. PCA  (colour = GROUP_COL, shape = GROUP_COL2)
# ---------------------------------------------------------------------------
evec <- fread(file.path(pca_dir, paste0(prefix, ".pca.eigenvec")))
eval <- scan(file.path(pca_dir, paste0(prefix, ".pca.eigenval")), quiet = TRUE)
# % variance relative to the sum of the computed eigenvalues (top N PCs only)
pve <- 100 * eval / sum(eval)

evec[, IID := as.character(IID)]
evec <- merge(evec, q_best[, .(IID, G1, G2, dom)], by = "IID", all.x = TRUE)
evec$Color <- if (has_g1) evec$G1 else paste0("Q", evec$dom)
color_title <- if (has_g1) group_col else sprintf("Cluster (K=%d)", best_k)
n_col <- length(unique(evec$Color))
if (n_col > length(qual_pal))
  warning(n_col, " colour groups but only ", length(qual_pal), " colours; some repeat")

pca_tag <- if (has_g1) paste0("_by_", safe_tag(group_col)) else ""
if (has_g2) {
  pca_tag <- paste0(pca_tag, "_", safe_tag(group_col2))
  n_shp <- length(unique(evec$G2))
  if (n_shp > length(shape_pal))
    warning(n_shp, " '", group_col2, "' levels but only ", length(shape_pal),
            " shapes; some repeat")
}

pca_plot <- function(a, b) {
  pts <- if (has_g2) {
    geom_point(aes(color = Color, shape = G2), size = 2.4, alpha = 0.85)
  } else {
    geom_point(aes(color = Color), size = 2.2, alpha = 0.85, shape = 16)
  }
  p <- ggplot(evec, aes(.data[[paste0("PC", a)]], .data[[paste0("PC", b)]])) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    pts +
    scale_color_manual(values = rep(qual_pal, length.out = n_col), name = color_title) +
    labs(x = sprintf("PC%d (%.1f%%)", a, pve[a]),
         y = sprintf("PC%d (%.1f%%)", b, pve[b])) +
    pop_theme
  if (has_g2) {
    p <- p + scale_shape_manual(values = rep(shape_pal, length.out = n_shp),
                                name = group_col2) +
      guides(color = guide_legend(order = 1, override.aes = list(size = 3.5)),
             shape = guide_legend(order = 2, override.aes = list(size = 3.5)))
  }
  p
}
pca_w <- if (has_g2) 9.5 else 8
save_plot(pca_plot(1, 2), paste0("pca_PC1_PC2", pca_tag), pca_w, 6.5)
if (length(eval) >= 4) save_plot(pca_plot(3, 4), paste0("pca_PC3_PC4", pca_tag), pca_w, 6.5)

scree <- data.table(PC = seq_along(eval), PVE = pve)
p_scree <- ggplot(scree, aes(PC, PVE)) +
  geom_col(fill = "#0072B2", width = 0.7) +
  geom_line(color = "#D55E00", linewidth = 0.9) +
  geom_point(color = "#D55E00", size = 2) +
  scale_x_continuous(breaks = scree$PC) +
  labs(x = "Principal component", y = "% variance (of top PCs)") +
  pop_theme
save_plot(p_scree, "pca_scree", 8, 5)

cat("\n============================================================\n")
cat(sprintf(" Plots written to: %s\n", plot_dir))
cat(sprintf(" Animals        : %d\n", nrow(ids)))
cat(sprintf(" Groups         : %s%s\n", if (has_g1) group_col else "(none)",
            if (has_g2) paste0(" x ", group_col2) else ""))
cat(sprintf(" Best K (CV)    : %d\n", best_k))
cat(sprintf(" PC1 / PC2      : %.1f%% / %.1f%% (of top %d PCs)\n", pve[1], pve[2], length(eval)))
cat("============================================================\n")