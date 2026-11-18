#!/usr/bin/env Rscript
# =============================================================================
# BLUPF90.GWAS.plots.R
#
# Generates four plots from a GWAS_summary.*.csv produced by
# BLUPF90.GWAS.summarize.sh:
#
#   1. Manhattan plot  — SNP effect (snp_effect)
#   2. Manhattan plot  — Variance explained (var_Explained), with dashed
#                        1% threshold line (override with VAR_THRESH env var)
#   3. Manhattan plot  — -log10(p-value) (neg_log_pval)
#   4. QQ plot         — observed vs expected -log10(p-value)
#
# Usage (command line):
#   Rscript BLUPF90.GWAS.plots.R  <summary_csv>  <out_dir>  [trait_label]
#
# Usage (from shell wrapper):
#   Called automatically by BLUPF90.GWAS.summarize.sh when sourced from
#   run_blup90.pipe.sh; env vars SUMMARY_CSV, PLOT_DIR, and TRAIT are used.
#
# Output (each plot written as BOTH .pdf and 600-dpi LZW-compressed .tiff):
#   <out_dir>/manhattan_snp_effect.{pdf,tiff}
#   <out_dir>/manhattan_var_explained.{pdf,tiff}
#   <out_dir>/manhattan_neg_log_pval.{pdf,tiff}
#   <out_dir>/qqplot_pval.{pdf,tiff}
#   <out_dir>/GWAS_plots_combined.{pdf,tiff}   (all 4 on one page)
#
# Color scheme is chosen automatically from the trait label:
#   average temperature traits -> blue/cyan palette
#   maximum temperature traits -> red/orange palette
#   anything else              -> neutral high-contrast palette
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# 0. Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

summary_csv  <- if (length(args) >= 1) args[1] else Sys.getenv("SUMMARY_CSV",  unset = NA)
out_dir      <- if (length(args) >= 2) args[2] else Sys.getenv("PLOT_DIR",     unset = ".")
trait_label  <- if (length(args) >= 3) args[3] else Sys.getenv("TRAIT",        unset = "Trait")

# Resolution for TIFF output (dpi); override with TIFF_DPI env var
tiff_dpi <- suppressWarnings(as.numeric(Sys.getenv("TIFF_DPI", unset = "600")))
if (is.na(tiff_dpi) || tiff_dpi <= 0) tiff_dpi <- 600

if (is.na(summary_csv) || !file.exists(summary_csv)) {
  stop("ERROR: summary CSV not found. Pass as first argument or set SUMMARY_CSV env var.\n",
       "  Usage: Rscript BLUPF90.GWAS.plots.R <summary_csv> [out_dir] [trait_label]")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Input  : %s\n", summary_csv))
cat(sprintf("Output : %s\n", out_dir))
cat(sprintf("Trait  : %s\n", trait_label))
cat(sprintf("TIFF   : %g dpi (LZW)\n", tiff_dpi))

# ---------------------------------------------------------------------------
# 1. Load and prepare data
# ---------------------------------------------------------------------------
gwas <- read.csv(summary_csv, stringsAsFactors = FALSE)

required_cols <- c("SNPID", "CHR", "POS", "snp_effect", "var_Explained", "neg_log_pval")
missing_cols  <- setdiff(required_cols, colnames(gwas))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing columns in summary CSV: ", paste(missing_cols, collapse = ", "))
}

# Remove rows missing all three metrics (uninformative)
gwas <- gwas[!(is.na(gwas$snp_effect) & is.na(gwas$var_Explained) & is.na(gwas$neg_log_pval)), ]

# Sort by chromosome then position
gwas <- gwas[order(gwas$CHR, gwas$POS), ]

# ---------------------------------------------------------------------------
# 2. Build cumulative x-axis (bp offset per chromosome)
#    Standard Manhattan layout: chromosomes laid end-to-end with a small gap
# ---------------------------------------------------------------------------
chr_info <- gwas %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(cum_start = cumsum(as.numeric(chr_len)) - chr_len)

gwas <- gwas %>%
  left_join(chr_info[, c("CHR", "cum_start")], by = "CHR") %>%
  mutate(bp_cum = POS + cum_start)

# Chromosome label positions (midpoint of each chromosome)
chr_labels <- gwas %>%
  group_by(CHR) %>%
  summarise(center = (min(bp_cum) + max(bp_cum)) / 2, .groups = "drop")

# ---------------------------------------------------------------------------
# 2b. Trait-specific, high-contrast color palettes
#     - alternating chromosome colors are dark/bright pairs so odd vs even
#       chromosomes separate cleanly even in grayscale print
#     - average- and maximum-temperature traits get distinct hue families
# ---------------------------------------------------------------------------
tl <- tolower(trait_label)

is_avg <- grepl("(^|[^a-z])(ave|avg|aver|average|mean)([^a-z]|$)", tl)
is_max <- grepl("(^|[^a-z])(max|maximum|peak)([^a-z]|$)", tl)

if (is_max) {
  # Maroon / Vermillion / Amber — a warm "heat" ramp, deliberately opposite the
  # cool ave-temp palette so the two traits are distinguishable at a glance.
  # Dark -> mid -> light lightness ramp keeps chromosomes separable in
  # grayscale and under deuteranopia/protanopia.
  trait_class    <- "maximum temperature (maroon / vermillion / amber)"
  chr_pal        <- c("#6A040F", "#E85D04", "#FFBA08")
  point_hi_color <- "#0072B2"                 # blue accent for significant SNPs
  qq_point_color <- "#6A040F"
  qq_line_color  <- "#0072B2"
  sugg_color     <- "#6A040F"
  gw_color       <- "#0072B2"
} else if (is_avg) {
  # Navy / Teal / Coral — cycled across chromosomes (colorblind-safe)
  # (teal lightened from #008080 so the three differ in lightness, not just hue)
  trait_class    <- "average temperature (navy / teal / coral)"
  chr_pal        <- c("#000080", "#00CED1", "#FF7F50")
  point_hi_color <- "#E6001A"                 # red accent for significant SNPs
  qq_point_color <- "#000080"
  qq_line_color  <- "#E6001A"
  sugg_color     <- "#000080"
  gw_color       <- "#E6001A"
} else {
  # Neutral fallback when the trait label matches neither pattern
  trait_class    <- "generic (charcoal / bright blue / amber)"
  chr_pal        <- c("#1B262C", "#3B9EFF", "#FFBF00")
  point_hi_color <- "#E6001A"
  qq_point_color <- "#36454F"
  qq_line_color  <- "#E6001A"
  sugg_color     <- "#36454F"
  gw_color       <- "#E6001A"
}
cat(sprintf("Palette: %s\n", trait_class))

# Cycle the three palette colors across chromosomes in order:
# chr 1 = color 1, chr 2 = color 2, chr 3 = color 3, chr 4 = color 1, ...
chr_sorted <- sort(unique(gwas$CHR))
chr_colors <- setNames(
  chr_pal[((seq_along(chr_sorted) - 1) %% length(chr_pal)) + 1],
  chr_sorted
)

# ---------------------------------------------------------------------------
# 3. Theme and helpers
# ---------------------------------------------------------------------------
# Global size controls — bump these to scale all text / line weights at once
base_font   <- 15    # base font size (pt)
axis_lwd    <- 1.1   # thickness of the x and y axis lines
tick_lwd    <- 1.1   # thickness of the axis tick marks

gwas_theme <- theme_minimal(base_size = base_font) +
  theme(
    plot.title        = element_blank(),
    plot.subtitle     = element_blank(),
    plot.caption      = element_blank(),
    plot.tag          = element_blank(),
    axis.title        = element_text(size = base_font + 2, face = "bold"),
    axis.title.x      = element_text(margin = margin(t = 8)),
    axis.title.y      = element_text(margin = margin(r = 8)),
    axis.text.x       = element_text(size = base_font - 3, angle = 45, hjust = 1, color = "black"),
    axis.text.y       = element_text(size = base_font - 1, color = "black"),
    axis.line         = element_line(color = "black", linewidth = axis_lwd, lineend = "square"),
    axis.ticks        = element_line(color = "black", linewidth = tick_lwd),
    axis.ticks.length = unit(4, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    panel.border      = element_blank(),
    legend.position   = "none",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.margin       = margin(12, 18, 12, 14)
  )

# --- Save helper: writes a PDF and a high-resolution TIFF of the same plot ---
save_plot <- function(plot, basename, width, height) {
  pdf_file  <- file.path(out_dir, paste0(basename, ".pdf"))
  tiff_file <- file.path(out_dir, paste0(basename, ".tiff"))

  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)

  ok <- tryCatch({
    ggsave(tiff_file, plot, width = width, height = height,
           device = "tiff", dpi = tiff_dpi,
           compression = "lzw", type = "cairo", bg = "white")
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    # Fallback for builds without cairo TIFF support
    ok <- tryCatch({
      ggsave(tiff_file, plot, width = width, height = height,
             device = "tiff", dpi = tiff_dpi, compression = "lzw", bg = "white")
      TRUE
    }, error = function(e) FALSE)
  }
  if (!ok) warning("Could not write TIFF: ", tiff_file)

  invisible(c(pdf_file, tiff_file))
}

# Significance thresholds for -log(p) plot
suggestive_thresh  <- -log10(1e-5)
genome_wide_thresh <- -log10(5e-8)

manhattan_plot <- function(df, y_col, y_label, title = NULL, subtitle = NULL,
                           highlight_thresh = NULL,
                           point_alpha = 0.85, point_size = 1.2) {

  df <- df[!is.na(df[[y_col]]), ]

  p <- ggplot(df, aes(x = bp_cum, y = .data[[y_col]], color = factor(CHR))) +
    geom_point(alpha = point_alpha, size = point_size, shape = 16) +
    scale_color_manual(values = chr_colors) +
    scale_x_continuous(
      breaks = chr_labels$center,
      labels = chr_labels$CHR,
      expand = expansion(mult = 0.01)
    ) +
    labs(
      x        = "Chromosome",
      y        = y_label
    ) +
    gwas_theme

  # Optional threshold lines
  if (!is.null(highlight_thresh)) {
    for (i in seq_along(highlight_thresh)) {
      p <- p + geom_hline(
        yintercept = highlight_thresh[[i]]$val,
        linetype   = highlight_thresh[[i]]$lty,
        color      = highlight_thresh[[i]]$col,
        linewidth  = 0.9
      )
    }
    # Highlight SNPs above the top threshold
    top_thresh <- max(sapply(highlight_thresh, `[[`, "val"))
    top_snps   <- df[df[[y_col]] >= top_thresh, ]
    if (nrow(top_snps) > 0) {
      p <- p + geom_point(
        data  = top_snps,
        size  = 2.1,
        shape = 16,
        alpha = 1
      )
    }
  }
  p
}

# ---------------------------------------------------------------------------
# 4. Plot 1 — SNP effect
# ---------------------------------------------------------------------------
cat("Generating plot 1: Manhattan — SNP effect...\n")

p1 <- manhattan_plot(
  df       = gwas,
  y_col    = "snp_effect",
  y_label  = "SNP effect"
)

save_plot(p1, "manhattan_snp_effect", width = 13, height = 5.2)

# ---------------------------------------------------------------------------
# 5. Plot 2 — Variance explained
# ---------------------------------------------------------------------------
cat("Generating plot 2: Manhattan — Variance explained...\n")

# Reviewer request: dashed 1% line + y-axis label that states the values are
# 1-Mb window statistics, not single-SNP values.
# Threshold (in %) can be overridden with the VAR_THRESH env var.
var_thresh <- suppressWarnings(as.numeric(Sys.getenv("VAR_THRESH", unset = "1")))
if (is.na(var_thresh) || var_thresh <= 0) var_thresh <- 1

p2 <- manhattan_plot(
  df               = gwas,
  y_col            = "var_Explained",
  y_label          = "% of additive genetic variance\nexplained by 1-Mb window",
  highlight_thresh = list(
    list(val = var_thresh, lty = "dashed", col = gw_color)
  )

)

save_plot(p2, "manhattan_var_explained", width = 13, height = 5.2)

# ---------------------------------------------------------------------------
# 6. Plot 3 — -log10(p-value) with significance lines
# ---------------------------------------------------------------------------
cat("Generating plot 3: Manhattan — -log10(p-value)...\n")

thresh_lines <- list(
  list(val = suggestive_thresh,  lty = "dashed", col = sugg_color),   # suggestive
  list(val = genome_wide_thresh, lty = "solid",  col = gw_color)      # genome-wide
)

p3 <- manhattan_plot(
  df                = gwas,
  y_col             = "neg_log_pval",
  y_label           = expression(-log[10](italic(p))),
  highlight_thresh  = thresh_lines
)

save_plot(p3, "manhattan_neg_log_pval", width = 13, height = 5.2)

# ---------------------------------------------------------------------------
# 7. Plot 4 — QQ plot
# ---------------------------------------------------------------------------
cat("Generating plot 4: QQ plot...\n")

pval_obs <- gwas$neg_log_pval[!is.na(gwas$neg_log_pval)]
n        <- length(pval_obs)

qq_df <- data.frame(
  expected = -log10(ppoints(n)),              # uniform distribution under H0
  observed = sort(pval_obs, decreasing = TRUE)
) %>% arrange(expected)

# Genomic inflation factor lambda
# Back-convert -log10(p) to chi-sq (1 df) to compute lambda
chi2_obs <- qchisq(10^(-pval_obs), df = 1, lower.tail = FALSE)
lambda   <- median(chi2_obs, na.rm = TRUE) / qchisq(0.5, df = 1)

p4 <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, color = qq_line_color,
              linewidth = 1.1, linetype = "dashed") +
  geom_point(color = qq_point_color, alpha = 0.7, size = 1.3, shape = 16) +
  labs(
    x        = expression("Expected " * -log[10](italic(p))),
    y        = expression("Observed " * -log[10](italic(p)))
  ) +
  gwas_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3)
  )

save_plot(p4, "qqplot_pval", width = 7, height = 7)

# ---------------------------------------------------------------------------
# 8. Combined 4-panel figure (PDF + high-resolution TIFF)
# ---------------------------------------------------------------------------
cat("Generating combined 4-panel figure...\n")

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
has_gridExtra <- requireNamespace("gridExtra", quietly = TRUE)

draw_combined <- function() {
  if (has_patchwork) {
    library(patchwork)
    combined <- (p1 / p2 / p3 / p4)
    print(combined)
  } else if (has_gridExtra) {
    library(gridExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 1)
  } else {
    print(p1); print(p2); print(p3); print(p4)
    message("TIP: Install 'patchwork' for a polished combined layout: install.packages('patchwork')")
  }
}

combined_pdf <- file.path(out_dir, "GWAS_plots_combined.pdf")
pdf(combined_pdf, width = 15, height = 18, onefile = TRUE)
draw_combined()
dev.off()

combined_tiff <- file.path(out_dir, "GWAS_plots_combined.tiff")
ok_tiff <- tryCatch({
  tiff(combined_tiff, width = 15, height = 18, units = "in",
       res = tiff_dpi, compression = "lzw", type = "cairo", bg = "white")
  draw_combined()
  dev.off()
  TRUE
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  tryCatch({
    tiff(combined_tiff, width = 15, height = 18, units = "in",
         res = tiff_dpi, compression = "lzw", bg = "white")
    draw_combined()
    dev.off()
    TRUE
  }, error = function(e2) { try(dev.off(), silent = TRUE); FALSE })
})
if (!ok_tiff) warning("Could not write TIFF: ", combined_tiff)

# ---------------------------------------------------------------------------
# 9. Summary
# ---------------------------------------------------------------------------
cat("\n============================================================\n")
cat(sprintf(" GWAS plots written to: %s\n", out_dir))
cat("   manhattan_snp_effect.pdf     / .tiff\n")
cat("   manhattan_var_explained.pdf  / .tiff\n")
cat("   manhattan_neg_log_pval.pdf   / .tiff\n")
cat("   qqplot_pval.pdf              / .tiff\n")
cat("   GWAS_plots_combined.pdf      / .tiff\n")
cat(sprintf(" TIFF resolution: %g dpi (LZW)\n", tiff_dpi))
cat(sprintf(" Palette        : %s\n", trait_class))
cat(sprintf(" SNPs plotted   : %d\n", nrow(gwas)))
cat(sprintf(" Chromosomes    : %s\n", paste(sort(unique(gwas$CHR)), collapse = ", ")))
cat(sprintf(" Lambda (GIF)   : %.4f\n", lambda))
cat("============================================================\n")