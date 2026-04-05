# metabolism-skills

Curated AI agent skills for metabolism and metabolomics research. Built for [Claude Code](https://claude.com/claude-code).

## Why

Metabolism-related skills are scattered across multiple repositories and represent only 1.6% of the scientific AI agent ecosystem ([Claw4Science, 2026](https://claw4science.org)). This project aggregates the best existing skills and fills gaps with original contributions, providing a single entry point for metabolomics researchers.

## Installation

```bash
# Add all skills
npx skills add dailycafi/metabolism-skills

# Or copy specific skills
cp -r skills/metabolomics-analysis/xcms-preprocessing ~/.claude/skills/
```

## Skill Catalog

### Mass Spectrometry Data Processing (`ms-data-processing/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [pyopenms](skills/ms-data-processing/pyopenms/) | pyOpenMS | LC-MS/MS data processing, untargeted metabolomics pipeline |
| [matchms](skills/ms-data-processing/matchms/) | matchms | Mass spectral similarity matching and compound identification |
| [gcms-processing](skills/ms-data-processing/gcms-processing/) | XCMS/matchms | GC-MS deconvolution, NIST matching, retention index calculation |
| [nmr-metabolomics](skills/ms-data-processing/nmr-metabolomics/) | nmrglue/speaq | 1D/2D NMR processing, binning, metabolite identification |
| [format-conversion](skills/ms-data-processing/format-conversion/) | msConvert | RAW/WIFF/mzML/mzXML format conversion |
| [spatial-metabolomics](skills/ms-data-processing/spatial-metabolomics/) | pyimzML | MALDI-MSI/DESI-MSI ion imaging and spatial analysis |

### Metabolic Modeling (`metabolic-modeling/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [cobrapy](skills/metabolic-modeling/cobrapy/) | COBRApy | Constraint-based modeling: FBA, FVA, gene knockouts, flux sampling |

### Metabolomics Analysis (`metabolomics-analysis/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [xcms-preprocessing](skills/metabolomics-analysis/xcms-preprocessing/) | XCMS (R) | Peak detection, alignment, gap filling |
| [msdial-preprocessing](skills/metabolomics-analysis/msdial-preprocessing/) | MS-DIAL | Untargeted metabolomics preprocessing |
| [metabolite-annotation](skills/metabolomics-analysis/metabolite-annotation/) | Various | Metabolite identification and annotation |
| [normalization-qc](skills/metabolomics-analysis/normalization-qc/) | R | QC, normalization, batch correction |
| [statistical-analysis](skills/metabolomics-analysis/statistical-analysis/) | R | Univariate/multivariate statistical analysis |
| [pathway-mapping](skills/metabolomics-analysis/pathway-mapping/) | MetaboAnalystR | KEGG/Reactome pathway enrichment |
| [targeted-analysis](skills/metabolomics-analysis/targeted-analysis/) | Various | Targeted metabolomics quantification |
| [lipidomics](skills/metabolomics-analysis/lipidomics/) | lipidr (R) | Lipidomics data analysis |
| [pharmacometabolomics](skills/metabolomics-analysis/pharmacometabolomics/) | RDKit | Drug metabolite profiling, CYP450 prediction, ADME |
| [clinical-metabolomics](skills/metabolomics-analysis/clinical-metabolomics/) | Python | IEM screening, newborn screening, clinical reporting |

### Pathway Analysis (`pathway-analysis/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [bioservices](skills/pathway-analysis/bioservices/) | bioservices (Python) | KEGG, Reactome, UniProt, STRING integration |

### Databases (`databases/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [hmdb](skills/databases/hmdb/) | REST API | Human Metabolome Database (220k+ metabolites) |
| [metabolights](skills/databases/metabolights/) | REST API | EMBL-EBI MetaboLights repository |
| [metabolomics-workbench](skills/databases/metabolomics-workbench/) | REST API | NIH Metabolomics Workbench (4,200+ studies) |
| [metabolomics-workbench-api](skills/databases/metabolomics-workbench-api/) | REST API | Metabolomics Workbench endpoint reference |
| [kegg-api](skills/databases/kegg-api/) | REST API | KEGG pathway/compound database |
| [reactome-api](skills/databases/reactome-api/) | REST API | Reactome pathway database |
| [string-api](skills/databases/string-api/) | REST API | STRING protein interaction network |

### Multi-Omics Integration (`multi-omics/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [mofa-integration](skills/multi-omics/mofa-integration/) | MOFA2 (R) | Multi-Omics Factor Analysis |
| [mixomics-analysis](skills/multi-omics/mixomics-analysis/) | mixOmics (R) | Multivariate integration (sPLS, DIABLO) |
| [data-harmonization](skills/multi-omics/data-harmonization/) | R | Cross-platform data harmonization |
| [similarity-network](skills/multi-omics/similarity-network/) | SNFtool (R) | Similarity Network Fusion |
| [mgwas-integration](skills/multi-omics/mgwas-integration/) | PLINK/coloc | Metabolite-gene association and Mendelian randomization |
| [microbiome-metabolomics](skills/multi-omics/microbiome-metabolomics/) | mmvec/mixOmics | Host-microbiome metabolite interactions |

### Systems Biology (`systems-biology/`)

| Skill | Tool | Description |
|-------|------|-------------|
| [flux-balance-analysis](skills/systems-biology/flux-balance-analysis/) | COBRApy | FBA for genome-scale models |
| [metabolic-reconstruction](skills/systems-biology/metabolic-reconstruction/) | CarveMe/gapseq | Genome-scale model reconstruction |
| [context-specific-models](skills/systems-biology/context-specific-models/) | COBRApy | Tissue/condition-specific models (GIMME/iMAT) |
| [gene-essentiality](skills/systems-biology/gene-essentiality/) | COBRApy | Gene essentiality prediction |
| [model-curation](skills/systems-biology/model-curation/) | memote | Model quality assessment and curation |
| [isotope-flux-analysis](skills/systems-biology/isotope-flux-analysis/) | IsoCor/COBRApy | 13C metabolic flux analysis and isotope correction |
| [network-visualization](skills/systems-biology/network-visualization/) | Escher/Cytoscape | Metabolic network maps with flux overlay |

## Credits & Attribution

This project curates skills from the open-source community. Every aggregated skill has an `upstream` field in its SKILL.md frontmatter linking to the original source.

| Source Repository | License | Skills |
|-------------------|---------|--------|
| [K-Dense-AI/claude-scientific-skills](https://github.com/K-Dense-AI/claude-scientific-skills) | MIT | cobrapy, pyopenms, bioservices, database API references |
| [wu-yc/LabClaw](https://github.com/wu-yc/LabClaw) | MIT | matchms, metabolomics-workbench |
| [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills) | MIT | 17 skills (metabolomics, multi-omics, systems-biology) |

See [LICENSES/NOTICE.md](LICENSES/NOTICE.md) for full per-skill attribution.

## Roadmap

See [docs/ROADMAP.md](docs/ROADMAP.md) for future planned skills and improvement ideas.

## Contributing

1. Fork the repo
2. Create a skill directory under the appropriate domain
3. Write a `SKILL.md` with real, runnable code blocks
4. Submit a pull request

## License

Original contributions are MIT licensed. Aggregated skills retain their original licenses. See [LICENSES/](LICENSES/) for details.
