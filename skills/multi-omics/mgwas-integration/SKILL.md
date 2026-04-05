---
name: bio-multi-omics-mgwas-integration
description: Metabolite-gene association analysis via metabolite genome-wide association studies (mGWAS). Covers PLINK-based GWAS with metabolite phenotypes, colocalization with coloc, Mendelian randomization with TwoSampleMR, and mQTL annotation. Use when linking genetic variants to metabolite levels or inferring causal metabolite-disease relationships.
tool_type: cli
primary_tool: plink
---

## Version Compatibility

Reference examples tested with: PLINK 1.9 (clumping), PLINK 2.0+ (association), coloc 5.2+, TwoSampleMR 0.5+, qqman 0.1.9+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `plink --version` or `plink2 --version`
- R: `packageVersion("<pkg>")` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolite-Gene Association (mGWAS)

**"Find genetic variants associated with metabolite levels"** -> Run genome-wide association with metabolite concentrations as quantitative phenotypes, then colocalize hits with disease loci and test causal relationships via Mendelian randomization.
- CLI: `plink2 --glm` for association testing
- R: `coloc::coloc.abf()` for colocalization, `TwoSampleMR` for causal inference

## Installation

```bash
# PLINK 2.0 (association testing)
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20231018.zip
unzip plink2_linux_x86_64_20231018.zip
chmod +x plink2
sudo mv plink2 /usr/local/bin/

# PLINK 1.9 (required for --clump, which is not available in PLINK 2.0)
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip
unzip plink_linux_x86_64_20231018.zip
chmod +x plink
sudo mv plink /usr/local/bin/

# R packages
Rscript -e '
install.packages(c("coloc", "remotes"), repos = "https://cran.r-project.org")
remotes::install_github("MRCIEU/TwoSampleMR")
'
```

## Prepare Metabolite Phenotype File

**Goal:** Format metabolite concentrations as PLINK-compatible phenotype files.

**Approach:** Rank-based inverse normal transform metabolite levels to ensure normality, then write a tab-delimited phenotype file with FID/IID columns.

```r
library(dplyr)

# Load metabolite quantification (samples x metabolites)
metab <- read.csv("metabolite_levels.csv")

# Inverse normal transformation for GWAS phenotypes
inv_normal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# Transform each metabolite column
metab_transformed <- metab %>%
  mutate(across(starts_with("M_"), inv_normal))

# Format for PLINK: FID, IID, then phenotype columns
pheno_file <- metab_transformed %>%
  transmute(
    FID = sample_id,
    IID = sample_id,
    glucose = M_glucose,
    lactate = M_lactate,
    alanine = M_alanine
  )

write.table(pheno_file, "metabolite_pheno.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
```

## Run mGWAS with PLINK 2.0

**Goal:** Test association between each SNP and metabolite levels genome-wide.

**Approach:** Run linear regression with covariates (age, sex, PCs) for each metabolite phenotype, then filter by genome-wide significance.

```bash
# Quality control on genotype data
plink2 \
  --bfile genotypes \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.02 \
  --mind 0.05 \
  --make-bed \
  --out genotypes_qc

# Run association for each metabolite phenotype
plink2 \
  --bfile genotypes_qc \
  --pheno metabolite_pheno.txt \
  --pheno-name glucose \
  --covar covariates.txt \
  --covar-name age,sex,PC1,PC2,PC3,PC4,PC5 \
  --glm hide-covar cols=+a1freq \
  --out mgwas_glucose

# Loop over multiple metabolites
for METAB in glucose lactate alanine; do
  plink2 \
    --bfile genotypes_qc \
    --pheno metabolite_pheno.txt \
    --pheno-name "$METAB" \
    --covar covariates.txt \
    --covar-name age,sex,PC1,PC2,PC3,PC4,PC5 \
    --glm hide-covar cols=+a1freq \
    --out "mgwas_${METAB}"
done

# Extract genome-wide significant mQTLs (p < 5e-8)
# Use header-based column selection (P column position varies by PLINK2 version)
head -1 mgwas_glucose.glucose.glm.linear > mgwas_glucose_significant.txt
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="P") pcol=i; next} $pcol+0 < 5e-8' \
  mgwas_glucose.glucose.glm.linear >> mgwas_glucose_significant.txt
```

## mQTL Identification and Annotation

**Goal:** Identify metabolic quantitative trait loci and annotate with biological context.

**Approach:** Clump significant hits into independent loci, then query HMDB and KEGG for metabolite pathway context.

```bash
# Clump associated SNPs into independent loci
# NOTE: --clump is only available in PLINK 1.9, not PLINK 2.0
plink \
  --bfile genotypes_qc \
  --clump mgwas_glucose.glucose.glm.linear \
  --clump-snp-field ID \
  --clump-p1 5e-8 \
  --clump-p2 1e-5 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --out mgwas_glucose_clumped
```

```r
library(httr)
library(jsonlite)

# Read clumped mQTL results
mqtls <- read.table("mgwas_glucose_clumped.clumped", header = TRUE)

# Query HMDB for metabolite information
query_hmdb <- function(metabolite_name) {
  url <- paste0("https://hmdb.ca/unearth/q?query=", URLencode(metabolite_name),
                "&searcher=metabolites&button=")
  resp <- GET(url)
  content(resp, "text", encoding = "UTF-8")
}

# Query KEGG for pathway context
query_kegg_compound <- function(compound_id) {
  url <- paste0("https://rest.kegg.jp/get/", compound_id)
  resp <- GET(url)
  content(resp, "text", encoding = "UTF-8")
}

# Example: get KEGG pathways for glucose (C00031)
kegg_info <- query_kegg_compound("cpd:C00031")

# Map mQTL SNPs to nearest genes using KEGG pathway links
kegg_pathways <- GET("https://rest.kegg.jp/link/pathway/hsa")
```

## Manhattan and QQ Plots

**Goal:** Visualize mGWAS results with publication-ready Manhattan and QQ plots.

**Approach:** Use qqman R package to generate standard GWAS visualizations with significance thresholds annotated.

```r
library(qqman)

# Read PLINK results
results <- read.table("mgwas_glucose.glucose.glm.linear",
                       header = TRUE, sep = "\t")

# Rename columns for qqman compatibility
gwas_data <- data.frame(
  SNP = results$ID,
  CHR = results$`#CHROM`,
  BP = results$POS,
  P = results$P
)

# Manhattan plot
png("manhattan_glucose.png", width = 1200, height = 600, res = 150)
manhattan(gwas_data,
          main = "mGWAS: Glucose Levels",
          ylim = c(0, 30),
          cex = 0.6,
          col = c("steelblue", "coral"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          annotatePval = 5e-8,
          annotateTop = TRUE)
dev.off()

# QQ plot
png("qq_glucose.png", width = 600, height = 600, res = 150)
qq(gwas_data$P, main = "QQ Plot: Glucose mGWAS")
dev.off()

# Calculate genomic inflation factor (lambda)
chisq <- qchisq(1 - gwas_data$P, 1)
lambda_gc <- median(chisq) / qchisq(0.5, 1)
message("Genomic inflation factor: ", round(lambda_gc, 3))
```

## Colocalization Analysis with coloc

**Goal:** Test whether mQTL signals share a causal variant with disease GWAS loci.

**Approach:** Extract regional summary statistics from both the mGWAS and a disease GWAS, then run Bayesian colocalization to compute posterior probabilities for shared vs distinct causal variants.

```r
library(coloc)

# Load mGWAS summary stats for a locus (chr2:100-101 Mb)
mqtl_region <- read.table("mgwas_glucose.glucose.glm.linear", header = TRUE) %>%
  filter(`#CHROM` == 2, POS >= 100e6, POS <= 101e6)

# Load disease GWAS summary stats for the same region
disease_region <- read.table("disease_gwas_chr2.txt", header = TRUE) %>%
  filter(CHR == 2, BP >= 100e6, BP <= 101e6)

# Prepare coloc datasets
dataset_mqtl <- list(
  beta = mqtl_region$BETA,
  varbeta = mqtl_region$SE^2,
  snp = mqtl_region$ID,
  position = mqtl_region$POS,
  type = "quant",
  N = 5000,
  MAF = mqtl_region$A1_FREQ,
  sdY = 1  # sdY = 1 when phenotype is inverse-normal transformed
)

dataset_disease <- list(
  beta = disease_region$BETA,
  varbeta = disease_region$SE^2,
  snp = disease_region$SNP,
  position = disease_region$BP,
  type = "cc",
  N = 50000,
  s = 0.3
)

# Run colocalization
# H4 posterior > 0.8 supports shared causal variant
result <- coloc.abf(dataset1 = dataset_mqtl, dataset2 = dataset_disease)
print(result$summary)
```

## Mendelian Randomization with TwoSampleMR

**Goal:** Test causal effect of metabolite levels on disease risk using mQTL instruments.

**Approach:** Select independent mQTLs as genetic instruments, extract their effects on a disease outcome from a separate GWAS, then apply MR methods (IVW, weighted median, MR-Egger) to estimate causal effects.

```r
library(TwoSampleMR)

# Format mQTL results as exposure data
exposure <- read_exposure_data(
  filename = "mgwas_glucose_significant.txt",
  sep = "\t",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "REF",  # Works when A1 == ALT; if A1 can equal REF, use --glm cols=+ax and "AX"
  pval_col = "P",
  eaf_col = "A1_FREQ"
)
exposure$exposure <- "Glucose"

# Clump instruments to ensure independence
exposure_clumped <- clump_data(exposure, clump_r2 = 0.001)

# Extract outcome data from IEU GWAS database
# Example: type 2 diabetes (ieu-a-26)
outcome <- extract_outcome_data(
  snps = exposure_clumped$SNP,
  outcomes = "ieu-a-26"
)

# Harmonize exposure and outcome
dat <- harmonise_data(exposure_clumped, outcome)

# Run MR analysis
mr_results <- mr(dat, method_list = c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode"
))
print(mr_results)

# Sensitivity analyses
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

# Visualization
mr_scatter_plot(mr_results, dat)
mr_forest_plot(mr_singlesnp(dat))
mr_funnel_plot(mr_singlesnp(dat))
```

## Key Databases

| Database | URL | Use |
|----------|-----|-----|
| mGWAS Atlas | https://metabolomics.helmholtz-muenchen.de/gwas/ | Published mGWAS results |
| GWAS Catalog | https://www.ebi.ac.uk/gwas/ | Metabolite trait associations |
| HMDB | https://hmdb.ca | Metabolite annotation |
| KEGG | https://www.kegg.jp | Pathway context |
| IEU OpenGWAS | https://gwas.mrcieu.ac.uk | MR outcome datasets |

## Best Practices

- Always inverse-normal transform metabolite phenotypes before association testing
- Include population structure covariates (PCs) to avoid confounding
- Apply genomic control if lambda > 1.05
- Use LD clumping (r2 < 0.1, 1000 kb window) before downstream analyses
- For MR: require F-statistic > 10 for each instrument to avoid weak instrument bias
- Run multiple MR methods; consistent results across methods strengthen causal claims
- Check for horizontal pleiotropy with MR-Egger intercept test

## Related Skills

- multi-omics/mixomics-analysis - Multi-omics integration
- multi-omics/data-harmonization - Harmonize omics datasets
- metabolomics-analysis/statistical-analysis - Differential metabolite analysis
