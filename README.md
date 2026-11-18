# Sheep Thermoregulation Genomics

Analysis pipelines for the genomic study of thermoregulation in Florida sheep: population structure (PCA, ADMIXTURE) and weighted single-step GBLUP genome-wide association (BLUPF90 suite), including QTL window annotation and functional enrichment on the Oar_v4.0 assembly.

This code accompanies:

> Cuellar et al. *(title TBD)*. *Journal of Animal Breeding and Genetics* (under revision).

---

## Repository structure

```
bin/scripts/
├── PopGen/                         # Population structure
│   ├── run_popgen.pipe.sh          # Driver: runs the full PopGen pipeline
│   ├── popgen.prep.sh              # Genotype QC and LD pruning (PLINK)
│   ├── PCA.sh                      # Principal component analysis (PLINK)
│   ├── ADMIXTURE.sh                # ADMIXTURE across K = 2–10 with cross-validation
│   └── popgen.plots.R              # PCA, scree, CV-error and admixture bar plots
│
└── BLUPF90/                        # Variance components and GWAS
    ├── run_blup90.pipe.sh          # Driver: runs the full BLUPF90 pipeline
    ├── setup_pheno_geno.R          # Builds phenotype file and keep list per trait
    ├── sheep_model.csv             # Model specification (traits, fixed/random effects)
    ├── BLUPF90.setup.geno.sh       # Subsets genotypes per trait
    ├── plink_2_blupf90.sh          # Converts PLINK files to BLUPF90 genotype format
    ├── AIREML.parameter            # Parameter template for AIREMLF90
    ├── AIREML.SNP.sh               # Variance component estimation (AIREMLF90)
    ├── GBLUP.parameter             # Parameter template for BLUPF90+/postGSf90
    ├── BLUPF90.GWAS.geno.add.sh    # Iterative weighted ssGBLUP GWAS (3 rounds)
    ├── BLUPF90.GWAS.summarize.sh   # Collects SNP effects, p-values and window variance
    ├── BLUPF90.GWAS.plots.R        # Manhattan and QQ plots
    └── find_QTL_genes_OarV4_fromWindows.R  # Gene annotation of significant windows + GO/KEGG ORA
```

---

## Software requirements

| Software | Version | Used for |
|---|---|---|
| PLINK | 1.9 | Genotype QC, LD pruning, PCA |
| ADMIXTURE | 1.3.0 | Ancestry estimation |
| BLUPF90 suite (renumf90, airemlf90, blupf90+, preGSf90, postGSf90) | 1.65  | Variance components, ssGBLUP GWAS |
| R | 4.4 | Data preparation, plotting, annotation |


Analyses were run on the University of Florida HiPerGator cluster (SLURM).

---

## Input data

Raw genotype and phenotype data are **not** included in this repository.

The pipelines expect:

- **Genotypes:** PLINK binary files (`.bed/.bim/.fam`) after QC, coordinates on Oar_v4.0
- **Phenotypes:** a CSV with animal ID, trait records, and covariates (e.g. farm, breed)
- **Annotation:** Oar_v4.0 GFF (NCBI) for gene annotation

Data availability: *(Phenotypes and Genotypes are private, can be requested)*

---

## Usage

Paths are set in `project_env.sh` (not included); edit the paths at the top of each driver script to match your system.

### 1. Population structure

```bash
bash bin/scripts/PopGen/run_popgen.pipe.sh
```

Steps: QC and LD pruning → PCA → ADMIXTURE (K = 2–10, CV error) → plots.

### 2. Variance components and GWAS

```bash
bash bin/scripts/BLUPF90/run_blup90.pipe.sh
```

Steps, per trait:

1. Set up phenotypes and genotypes, convert to BLUPF90 format
2. Estimate variance components with AIREMLF90
3. Run weighted ssGBLUP GWAS for three iterations, updating SNP weights each round
4. Summarize SNP effects, p-values, and variance explained by genomic windows
5. Produce Manhattan and QQ plots
6. Annotate significant windows with Oar_v4.0 genes and run GO/KEGG over-representation analysis

---

## Outputs

Each run writes to `results/`:

- `results/PopGen/`: PCA eigenvectors, ADMIXTURE Q/P matrices, CV error, figures
- `results/BLUPF90/aireml/<TRAIT>/`: variance components and heritability
- `results/BLUPF90/GBLUP_GWAS/<TRAIT>/`: per-round BLUPF90 output, GWAS summaries, significant windows, plots, QTL gene tables, enrichment results

---

## Citation

If you use this code, please cite the paper above. *(DOI to be added on publication.)*

## Contact

Gabriel A. Zayas Santiago, Ph.D.: zayas@wisc.edu | gazspr@gmail.com
