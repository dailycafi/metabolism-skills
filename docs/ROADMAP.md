# MetaboClaw Roadmap

## Current Coverage (v0.1.0)

26 skills across 6 domains, aggregated from 3 open-source repositories.

| Domain | Skills | Source |
|--------|--------|--------|
| ms-data-processing | pyopenms, matchms | K-Dense, LabClaw |
| metabolic-modeling | cobrapy | K-Dense |
| metabolomics-analysis | xcms-preprocessing, msdial-preprocessing, metabolite-annotation, normalization-qc, statistical-analysis, pathway-mapping, targeted-analysis, lipidomics | bioSkills |
| pathway-analysis | bioservices (KEGG/Reactome/STRING) | K-Dense |
| databases | metabolomics-workbench (x2), kegg-api, reactome-api, string-api | K-Dense, LabClaw |
| multi-omics | mofa-integration, mixomics-analysis, data-harmonization, similarity-network | bioSkills |
| systems-biology | flux-balance-analysis, metabolic-reconstruction, context-specific-models, gene-essentiality, model-curation | bioSkills |

## Gaps to Fill (Original Contributions)

### High Priority

- [ ] **HMDB database lookup** -- Human Metabolome Database (220k+ metabolites), REST API access, metabolite annotation by HMDB ID, structure search
- [ ] **GC-MS data processing** -- Deconvolution, NIST/AMDIS library matching, retention index calculation, specific to volatile/semi-volatile metabolites
- [ ] **NMR metabolomics** -- 1D/2D NMR processing (Chenomx, MetaboMiner, rNMR), spectral binning, metabolite identification from NMR spectra
- [ ] **MetaboLights database** -- EMBL-EBI metabolomics repository API, study retrieval, ISA-Tab format parsing

### Medium Priority

- [ ] **13C metabolic flux analysis (13C-MFA)** -- Isotope labeling experiments, INCA/FluxML/OpenFLUX integration, mass isotopomer distribution analysis
- [ ] **Metabolomics data format conversion** -- mzML/mzXML/RAW/WIFF interconversion, msConvert wrapper, vendor format handling
- [ ] **Metabolite-gene association (mGWAS)** -- Metabolic QTL analysis, GWAS integration with metabolomics, Mendelian randomization for metabolites
- [ ] **Spatial metabolomics** -- MALDI-MSI, DESI-MSI data processing, ion image generation, molecular histology

### Low Priority

- [ ] **Metabolic network visualization** -- Escher/MetExplore metabolic maps, pathway-level flux visualization, publication-ready diagrams
- [ ] **Pharmacometabolomics** -- Drug metabolite profiling, CYP450 substrate prediction, ADME property prediction
- [ ] **Microbiome metabolomics** -- Host-microbiome metabolite interactions, SCFA analysis, bile acid profiling
- [ ] **Clinical metabolomics reporting** -- IEM (inborn errors of metabolism) screening, reference ranges, clinical decision support

## Contributing

To add a new skill:

1. Create a directory under the appropriate domain in `skills/`
2. Write a `SKILL.md` with this frontmatter:
   ```yaml
   ---
   name: your-skill-name
   description: "One-line description"
   tool_type: python | r | cli
   primary_tool: main-package-name
   ---
   ```
3. Include real, runnable code blocks (not pseudocode)
4. Submit a pull request
