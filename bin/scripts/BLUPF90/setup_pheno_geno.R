library(data.table)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop(paste(
    "Incorrect number of arguments.\n",
    "Usage:\n",
    "  Rscript msetup_pheno_geno.R",
    " <pheno> <metadata> <model> <plinkfam> <trait> <out_dir>\n",
    sep = ""
  ), call. = FALSE)
}

pheno <- args[1]
meta <- args[2]
model <- args[3]
fam <- args[4]
trait <- args[5]
out_dir <- args[6]

# library(data.table)
# pheno <- "/blue/mateescu/gzayas97/For_UF/Sheep_Project/data/Sheep_Master_File.csv"
# meta <- "/blue/mateescu/gzayas97/For_UF/Sheep_Project/bin/scripts/BLUPF90/sheep_model.csv"
# model <- "~ REP"
# fam <- "/blue/mateescu/gzayas97/For_UF/Sheep_Project/data/UFL_Sheep_plink.fam"
# trait <- "AVE"
# out_dir <- "/blue/mateescu/gzayas97/For_UF/Sheep_Project/results/BLUPF90/Setup/AVE"


pheno_dt <- fread(pheno)
meta_dt <- fread(meta)
fam_file <- fread(fam)[,1:2]
geno_ids <- fam_file[[2]]

geno_id_col <- meta_dt$ColN[meta_dt$What == "GenoID"]
pheno_id_col <- meta_dt$ColN[meta_dt$Name == trait]

# Parse model string to extract term names
model_terms <- attr(terms(as.formula(model)), "term.labels")

# Look up ColN for each term
model_cols <- meta_dt$ColN[meta_dt$Name %in% model_terms]
cols_to_check <- c(pheno_id_col, model_cols)
complete_rows <- complete.cases(pheno_dt[, ..cols_to_check])
pheno_dt_clean <- pheno_dt[complete_rows]
pheno_ids <- na.omit(pheno_dt_clean [[ meta_dt$ColN[meta_dt$What == "GenoID"]]])

intersection_ids <- intersect(geno_ids,pheno_ids )

# Get intersection individuals and save keep file
keep_dt <- fam_file[fam_file[[1]] %in% intersection_ids, ]  # FamID and IID
fwrite(keep_dt, file.path(out_dir, "keep.txt"), sep = " ", col.names = FALSE)

# Filter pheno to intersection only
pheno_dt_final <- pheno_dt_clean[pheno_dt_clean[[geno_id_col]] %in% intersection_ids]

# Build BLUPF90 phenotype file: GenoID, trait, fixed effects
blup_cols <- c(geno_id_col, pheno_id_col, model_cols)
pheno_blup <- pheno_dt_final[, ..blup_cols]

fwrite(pheno_blup, file.path(out_dir, "phenotype.txt"), sep = " ", col.names = FALSE)
