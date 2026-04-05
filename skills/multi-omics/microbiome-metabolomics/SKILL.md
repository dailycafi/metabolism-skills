---
name: bio-multi-omics-microbiome-metabolomics
description: Host-microbiome metabolite interaction analysis integrating 16S/shotgun metagenomics with untargeted and targeted metabolomics. Covers SCFA quantification, bile acid profiling, tryptophan pathway analysis, and microbe-metabolite co-occurrence modeling. Use when linking microbial community composition to metabolic phenotypes.
tool_type: mixed
primary_tool: mmvec
---

## Version Compatibility

Reference examples tested with: scikit-bio 0.6+, pandas 2.2+, scipy 1.12+, matplotlib 3.8+, seaborn 0.13+, mmvec 1.0+, mixOmics 6.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Microbiome Metabolomics

**"Link my gut microbiome data with metabolomics profiles"** → Integrate microbial abundance tables with metabolite measurements to identify host-microbiome metabolic interactions.
- Python: scikit-bio for diversity, mmvec for co-occurrence, pandas for integration
- R: mixOmics for sPLS correlation, MIMOSA2 for community metabolic modeling

## Installation

```bash
# Python environment
pip install scikit-bio pandas scipy matplotlib seaborn networkx

# mmvec (neural network microbe-metabolite co-occurrence)
pip install mmvec

# R packages
Rscript -e '
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("mixOmics")
  install.packages(c("pheatmap", "igraph", "ggplot2"))
'
```

## SCFA Quantification from Targeted LC-MS/GC-MS

**Goal:** Quantify short-chain fatty acids (acetate, propionate, butyrate) from targeted assay data using calibration curves and internal standard normalization.

**Approach:** Load standard curve and sample peak areas, fit weighted regression per SCFA, back-calculate concentrations with IS correction.

```python
import pandas as pd
import numpy as np
from scipy.stats import linregress

scfa_names = ['acetate', 'propionate', 'butyrate']

standards = pd.DataFrame({
    'concentration_uM': [0, 10, 50, 100, 250, 500, 1000],
    'acetate_area': [200, 15000, 72000, 148000, 365000, 740000, 1480000],
    'propionate_area': [150, 12000, 58000, 120000, 295000, 600000, 1190000],
    'butyrate_area': [180, 13500, 65000, 135000, 330000, 670000, 1340000],
    'istd_area': [500000, 498000, 502000, 497000, 501000, 499000, 503000],
})

calibrations = {}
for scfa in scfa_names:
    ratio = standards[f'{scfa}_area'] / standards['istd_area']
    conc = standards['concentration_uM'].values
    mask = conc > 0
    slope, intercept, r_value, _, _ = linregress(conc[mask], ratio.values[mask])
    calibrations[scfa] = {'slope': slope, 'intercept': intercept, 'r2': r_value**2}

samples = pd.read_csv('scfa_peak_areas.csv')
for scfa in scfa_names:
    ratio = samples[f'{scfa}_area'] / samples['istd_area']
    cal = calibrations[scfa]
    samples[f'{scfa}_uM'] = np.maximum((ratio - cal['intercept']) / cal['slope'], 0)

samples[['sample_id'] + [f'{s}_uM' for s in scfa_names]].to_csv('scfa_concentrations.csv', index=False)
```

## Bile Acid Profiling

**Goal:** Classify bile acids into primary vs secondary and conjugated vs unconjugated, then summarize profiles per sample.

**Approach:** Map measured bile acids to biological categories, compute group totals and ratios for downstream analysis.

```python
import pandas as pd

bile_acid_classes = {
    'CA': {'type': 'primary', 'conjugation': 'unconjugated'},
    'CDCA': {'type': 'primary', 'conjugation': 'unconjugated'},
    'GCA': {'type': 'primary', 'conjugation': 'glycine'},
    'TCA': {'type': 'primary', 'conjugation': 'taurine'},
    'GCDCA': {'type': 'primary', 'conjugation': 'glycine'},
    'TCDCA': {'type': 'primary', 'conjugation': 'taurine'},
    'DCA': {'type': 'secondary', 'conjugation': 'unconjugated'},
    'LCA': {'type': 'secondary', 'conjugation': 'unconjugated'},
    'UDCA': {'type': 'secondary', 'conjugation': 'unconjugated'},
    'GDCA': {'type': 'secondary', 'conjugation': 'glycine'},
    'TDCA': {'type': 'secondary', 'conjugation': 'taurine'},
    'GLCA': {'type': 'secondary', 'conjugation': 'glycine'},
    'TLCA': {'type': 'secondary', 'conjugation': 'taurine'},
}

ba_data = pd.read_csv('bile_acid_concentrations.csv', index_col='sample_id')
class_df = pd.DataFrame(bile_acid_classes).T

profile = pd.DataFrame(index=ba_data.index)
measured = [col for col in ba_data.columns if col in class_df.index]

for category in ['primary', 'secondary']:
    members = class_df[class_df['type'] == category].index
    cols = [m for m in members if m in measured]
    profile[f'total_{category}'] = ba_data[cols].sum(axis=1)

for conj in ['unconjugated', 'glycine', 'taurine']:
    members = class_df[class_df['conjugation'] == conj].index
    cols = [m for m in members if m in measured]
    profile[f'total_{conj}'] = ba_data[cols].sum(axis=1)

profile['primary_secondary_ratio'] = profile['total_primary'] / profile['total_secondary'].replace(0, np.nan)
profile['conjugated_unconjugated_ratio'] = (
    (profile['total_glycine'] + profile['total_taurine']) /
    profile['total_unconjugated'].replace(0, np.nan)
)
profile.to_csv('bile_acid_profiles.csv')
```

## Tryptophan Metabolite Pathway Analysis

**Goal:** Map measured tryptophan metabolites to the indole and kynurenine pathway branches and compute pathway activity scores.

**Approach:** Assign metabolites to pathway branches, sum normalized abundances per branch, and compute the kynurenine-to-tryptophan ratio (KTR) as an IDO/TDO activity proxy.

```python
import pandas as pd
import numpy as np

trp_pathways = {
    'indole': ['indole', 'indole_3_acetate', 'indole_3_propionate',
               'indole_3_aldehyde', 'indoxyl_sulfate', 'tryptamine'],
    'kynurenine': ['kynurenine', 'kynurenic_acid', '3_hydroxykynurenine',
                   'anthranilic_acid', 'quinolinic_acid', 'picolinic_acid'],
    'serotonin': ['serotonin', '5_hydroxyindoleacetic_acid', 'melatonin'],
}

trp_data = pd.read_csv('tryptophan_metabolites.csv', index_col='sample_id')

trp_data_norm = trp_data.div(trp_data.sum(axis=1), axis=0)

pathway_scores = pd.DataFrame(index=trp_data.index)
for pathway, metabolites in trp_pathways.items():
    cols = [m for m in metabolites if m in trp_data_norm.columns]
    pathway_scores[f'{pathway}_score'] = trp_data_norm[cols].sum(axis=1)

if 'kynurenine' in trp_data.columns and 'tryptophan' in trp_data.columns:
    pathway_scores['KTR'] = trp_data['kynurenine'] / trp_data['tryptophan'].replace(0, np.nan)

pathway_scores.to_csv('tryptophan_pathway_scores.csv')
```

## 16S/Metagenomics Data Loading and Diversity

**Goal:** Load microbial abundance tables and compute alpha/beta diversity metrics for integration with metabolomics.

**Approach:** Read OTU/ASV tables, compute Shannon diversity and Bray-Curtis dissimilarity using scikit-bio.

```python
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity, beta_diversity
from skbio import DistanceMatrix

otu_table = pd.read_csv('otu_table.csv', index_col='sample_id')

alpha_div = pd.DataFrame({
    'shannon': alpha_diversity('shannon', otu_table.values, otu_table.index),
    'observed_otus': alpha_diversity('observed_otus', otu_table.values, otu_table.index),
    'chao1': alpha_diversity('chao1', otu_table.values, otu_table.index),
}, index=otu_table.index)

bc_dm = beta_diversity('braycurtis', otu_table.values, otu_table.index)

from skbio.stats.ordination import pcoa
pcoa_results = pcoa(bc_dm)
pcoa_coords = pcoa_results.samples[['PC1', 'PC2', 'PC3']]
pcoa_coords.index = otu_table.index
```

## Microbiome-Metabolome Data Integration

**Goal:** Align microbial abundance and metabolomics matrices on shared samples for joint analysis.

**Approach:** Intersect sample IDs, apply CLR transform to compositional microbiome data, log-transform metabolomics, and compute Spearman correlations.

```python
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.special import log1p

otu_table = pd.read_csv('otu_table.csv', index_col='sample_id')
metab_table = pd.read_csv('metabolomics.csv', index_col='sample_id')

common_samples = otu_table.index.intersection(metab_table.index)
otu_aligned = otu_table.loc[common_samples]
metab_aligned = metab_table.loc[common_samples]

def clr_transform(df):
    pseudocount = 0.5
    log_data = np.log(df + pseudocount)
    geometric_mean = log_data.mean(axis=1)
    return log_data.sub(geometric_mean, axis=0)

otu_clr = clr_transform(otu_aligned)
metab_log = np.log1p(metab_aligned)

n_taxa = otu_clr.shape[1]
n_metab = metab_log.shape[1]
corr_matrix = np.zeros((n_taxa, n_metab))
pval_matrix = np.zeros((n_taxa, n_metab))

for i in range(n_taxa):
    for j in range(n_metab):
        rho, pval = spearmanr(otu_clr.iloc[:, i], metab_log.iloc[:, j])
        corr_matrix[i, j] = rho
        pval_matrix[i, j] = pval

corr_df = pd.DataFrame(corr_matrix, index=otu_clr.columns, columns=metab_log.columns)
pval_df = pd.DataFrame(pval_matrix, index=otu_clr.columns, columns=metab_log.columns)
```

## mmvec: Neural Network Co-Occurrence

**Goal:** Model microbe-metabolite co-occurrence probabilities using mmvec neural network embeddings.

**Approach:** Prepare BIOM-format inputs, run mmvec to learn conditional probabilities, then extract co-occurrence ranks.

```bash
# Convert tables to BIOM format
biom convert -i otu_table.tsv -o microbes.biom --table-type="OTU table" --to-hdf5
biom convert -i metabolite_table.tsv -o metabolites.biom --table-type="OTU table" --to-hdf5

# Run mmvec
mmvec paired-omics \
    --microbe-file microbes.biom \
    --metabolite-file metabolites.biom \
    --epochs 1000 \
    --batch-size 64 \
    --latent-dim 3 \
    --input-prior 0.01 \
    --output-prior 0.01 \
    --learning-rate 1e-3 \
    --summary-dir mmvec_output
```

```python
import pandas as pd

conditionals = pd.read_csv('mmvec_output/conditionals.tsv', sep='\t', index_col=0)

top_associations = []
for microbe in conditionals.index:
    ranked = conditionals.loc[microbe].sort_values(ascending=False)
    for metabolite in ranked.head(5).index:
        top_associations.append({
            'microbe': microbe,
            'metabolite': metabolite,
            'conditional_prob': ranked[metabolite],
        })

top_df = pd.DataFrame(top_associations)
top_df.to_csv('top_microbe_metabolite_associations.csv', index=False)
```

## sPLS Microbiome-Metabolome Correlation (R/mixOmics)

**Goal:** Identify correlated features between microbiome and metabolome using sparse PLS in R.

**Approach:** Load CLR-transformed microbiome and log-transformed metabolomics matrices, tune sPLS, and extract selected feature pairs.

```r
library(mixOmics)

X_microbiome <- as.matrix(read.csv('otu_clr.csv', row.names = 1))
X_metabolome <- as.matrix(read.csv('metabolomics_log.csv', row.names = 1))

common <- intersect(rownames(X_microbiome), rownames(X_metabolome))
X_microbiome <- X_microbiome[common, ]
X_metabolome <- X_metabolome[common, ]

tune_spls <- perf(spls(X_microbiome, X_metabolome, ncomp = 5),
                  validation = 'Mfold', folds = 5, nrepeat = 10)

spls_result <- spls(X_microbiome, X_metabolome, ncomp = 3,
                    keepX = c(30, 30, 30), keepY = c(20, 20, 20))

plotIndiv(spls_result, comp = c(1, 2), rep.space = 'XY-variate')
plotVar(spls_result, comp = c(1, 2), var.names = TRUE, cex = c(3, 3))
cim(spls_result, comp = 1, xlab = 'Metabolites', ylab = 'Microbes',
    margins = c(8, 12))

selected_microbes <- selectVar(spls_result, comp = 1)$X$name
selected_metabolites <- selectVar(spls_result, comp = 1)$Y$name
```

## Visualization: Paired Heatmaps and Correlation Networks

**Goal:** Generate publication-quality paired heatmaps and network plots of microbe-metabolite associations.

**Approach:** Build clustered heatmaps with significance masking and networkx correlation graphs with community detection.

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import false_discovery_control

from statsmodels.stats.multitest import multipletests

_, pval_corrected, _, _ = multipletests(pval_df.values.flatten(), method='fdr_bh')
pval_fdr = pd.DataFrame(
    pval_corrected.reshape(pval_df.shape),
    index=pval_df.index,
    columns=pval_df.columns
)

sig_mask = pval_fdr < 0.05
top_taxa = corr_df[sig_mask].abs().max(axis=1).nlargest(25).index
top_metab = corr_df[sig_mask].abs().max(axis=0).nlargest(25).index

fig, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(
    corr_df.loc[top_taxa, top_metab],
    cmap='RdBu_r', center=0, vmin=-0.8, vmax=0.8,
    xticklabels=True, yticklabels=True,
    linewidths=0.5, ax=ax,
)
ax.set_xlabel('Metabolites')
ax.set_ylabel('Microbial Taxa')
ax.set_title('Microbiome-Metabolome Spearman Correlations (FDR < 0.05)')
plt.tight_layout()
plt.savefig('microbiome_metabolome_heatmap.png', dpi=300)

G = nx.Graph()
for taxon in top_taxa:
    G.add_node(taxon, node_type='microbe')
for metab in top_metab:
    G.add_node(metab, node_type='metabolite')

for taxon in top_taxa:
    for metab in top_metab:
        if pval_fdr.loc[taxon, metab] < 0.05 and abs(corr_df.loc[taxon, metab]) > 0.4:
            G.add_edge(taxon, metab, weight=corr_df.loc[taxon, metab])

pos = nx.spring_layout(G, k=2, seed=42)
node_colors = ['#2196F3' if G.nodes[n]['node_type'] == 'microbe' else '#FF9800' for n in G.nodes()]
edge_colors = ['#d32f2f' if G.edges[e]['weight'] < 0 else '#388e3c' for e in G.edges()]

fig, ax = plt.subplots(figsize=(14, 12))
nx.draw_networkx(G, pos, node_color=node_colors, edge_color=edge_colors,
                 node_size=600, font_size=7, width=1.5, ax=ax)
ax.set_title('Microbe-Metabolite Correlation Network')
plt.tight_layout()
plt.savefig('correlation_network.png', dpi=300)
```

## Best Practices

- **Compositionality:** Always apply CLR or other compositional transform to microbiome relative abundance data before correlation analysis; raw relative abundances produce spurious correlations.
- **Multiple testing:** Apply FDR correction (Benjamini-Hochberg) when computing many pairwise correlations between taxa and metabolites.
- **Sample matching:** Verify that microbiome and metabolomics samples are from the same biological specimens; mismatched sample IDs silently corrupt results.
- **SCFA specificity:** Use deuterated internal standards (d3-acetate, d5-propionate, d7-butyrate) for SCFA quantification; matrix effects in fecal samples are substantial.
- **Bile acid isomers:** Confirm chromatographic separation of alpha/beta-MCA and UDCA/HDCA isomers before quantification; co-elution is common on C18 columns.
- **mmvec training:** Monitor training loss convergence; underfitting (too few epochs) yields random co-occurrence rankings.
- **Batch effects:** If samples were collected across sites or sequenced in different runs, include batch as a covariate or apply ComBat correction before integration.

## Related Skills

- mixomics-analysis - Supervised multi-omics integration with DIABLO
- targeted-analysis - Calibration curves and absolute quantification
- normalization-qc - QC-based normalization for metabolomics
- data-harmonization - Preprocess before multi-omics integration
