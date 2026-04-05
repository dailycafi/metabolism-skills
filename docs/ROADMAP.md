# metabolism-skills Roadmap

## Current Coverage (v1.0.0)

38 skills across 7 domains: 22 aggregated from open-source repositories + 12 original contributions.

| Domain | Skills | Count |
|--------|--------|-------|
| ms-data-processing | pyopenms, matchms, gcms-processing, nmr-metabolomics, format-conversion, spatial-metabolomics | 6 |
| metabolic-modeling | cobrapy | 1 |
| metabolomics-analysis | xcms-preprocessing, msdial-preprocessing, metabolite-annotation, normalization-qc, statistical-analysis, pathway-mapping, targeted-analysis, lipidomics, pharmacometabolomics, clinical-metabolomics | 10 |
| pathway-analysis | bioservices (KEGG/Reactome/STRING) | 1 |
| databases | hmdb, metabolights, metabolomics-workbench (x2), kegg-api, reactome-api, string-api | 7 |
| multi-omics | mofa-integration, mixomics-analysis, data-harmonization, similarity-network, mgwas-integration, microbiome-metabolomics | 6 |
| systems-biology | flux-balance-analysis, metabolic-reconstruction, context-specific-models, gene-essentiality, model-curation, isotope-flux-analysis, network-visualization | 7 |

## Completed (v1.0.0)

All originally planned gaps have been filled:

- [x] HMDB database lookup
- [x] GC-MS data processing
- [x] NMR metabolomics
- [x] MetaboLights database
- [x] 13C metabolic flux analysis
- [x] Metabolomics data format conversion
- [x] Metabolite-gene association (mGWAS)
- [x] Spatial metabolomics
- [x] Metabolic network visualization
- [x] Pharmacometabolomics
- [x] Microbiome metabolomics
- [x] Clinical metabolomics reporting

## Future Ideas

- [ ] **Stable isotope resolved metabolomics (SIRM)** -- Multi-tracer experiments, isotopologue analysis
- [ ] **Single-cell metabolomics** -- SpaceM, mass cytometry metabolite panels
- [ ] **Metabolomics data repositories** -- MassIVE, GNPS molecular networking
- [ ] **Exposomics** -- Environmental metabolite profiling, biomonitoring
- [ ] **Fluxomics visualization** -- Sankey diagrams for metabolic flux distributions
- [ ] **Metabolite structure elucidation** -- SIRIUS/CSI:FingerID, molecular formula prediction from MS/MS

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
