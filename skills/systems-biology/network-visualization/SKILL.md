---
name: bio-systems-biology-network-visualization
description: Metabolic network visualization using Escher, py4cytoscape, KEGG API, and NetworkX. Overlay flux data from FBA onto pathway maps, automate Cytoscape sessions, color KEGG pathways, and generate publication-ready figures. Use when visualizing metabolic pathways, flux distributions, or network topology.
tool_type: python
primary_tool: escher
---

## Version Compatibility

Reference examples tested with: escher 1.7+, py4cytoscape 1.9+, networkx 3.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolic Network Visualization

**"Visualize metabolic fluxes on pathway maps"** -> Render genome-scale flux distributions on interactive pathway maps with Escher, automate Cytoscape network layouts, color KEGG pathways by expression or flux, and build custom metabolic networks with NetworkX.
- Python: `escher.Builder()` for pathway maps, `py4cytoscape` for Cytoscape automation, `networkx` for graph analysis

## Installation

```bash
pip install escher py4cytoscape networkx matplotlib requests pandas

# For Jupyter integration
pip install jupyterlab
jupyter labextension install @jupyter-widgets/jupyterlab-manager

# Cytoscape must be running locally for py4cytoscape
# Download from https://cytoscape.org/download.html
```

## Escher: Interactive Pathway Maps with Flux Overlay

**Goal:** Display FBA flux results on a metabolic pathway map in Jupyter.

**Approach:** Load a COBRApy model, run FBA, then pass the flux dictionary to Escher's Builder to render an interactive SVG map with reaction arrows scaled by flux magnitude.

```python
import escher
import cobra

# Load model and run FBA
model = cobra.io.load_model("textbook")
solution = model.optimize()

# Build flux dictionary from solution
flux_data = {rxn.id: solution.fluxes[rxn.id] for rxn in model.reactions}

# Create Escher map with flux overlay
builder = escher.Builder(
    map_name="e_coli_core.Core metabolism",
    model=model,
    reaction_data=flux_data,
    reaction_scale=[
        {"type": "min", "color": "#2166ac", "size": 4},
        {"type": "value", "value": 0, "color": "#f7f7f7", "size": 8},
        {"type": "max", "color": "#b2182b", "size": 20},
    ],
    reaction_styles=["color", "size", "text"],
    identifiers_on_map="bigg_id",
)

# Display in Jupyter
builder.display_in_notebook()
```

## Escher: Metabolite Concentration Overlay

**Goal:** Map metabolomics measurements onto pathway nodes.

**Approach:** Provide a metabolite-level dictionary to Escher's Builder to color and size metabolite nodes by measured concentration.

```python
import escher

# Metabolite data (BiGG IDs -> measured concentrations in mM)
metabolite_data = {
    "glc__D_c": 5.2,
    "g6p_c": 0.8,
    "f6p_c": 0.3,
    "fdp_c": 1.5,
    "pyr_c": 0.9,
    "lac__D_c": 2.1,
    "atp_c": 4.8,
    "adp_c": 1.2,
}

builder = escher.Builder(
    map_name="e_coli_core.Core metabolism",
    metabolite_data=metabolite_data,
    metabolite_styles=["color", "size", "text"],
    metabolite_scale=[
        {"type": "min", "color": "#ffffcc", "size": 10},
        {"type": "median", "color": "#fd8d3c", "size": 20},
        {"type": "max", "color": "#800026", "size": 30},
    ],
)

builder.display_in_notebook()
```

## Escher: Export Publication Figures

**Goal:** Save Escher maps as static images for manuscripts.

**Approach:** Use the Builder's save methods to export HTML (for interactive supplementary) or embed in a figure workflow.

```python
import escher

builder = escher.Builder(
    map_name="e_coli_core.Core metabolism",
    reaction_data=flux_data,
    reaction_styles=["color", "size"],
)

# Save as standalone HTML
builder.save_html("flux_map.html")

# Save as embedded HTML for supplementary materials
builder.save_html("flux_map_embedded.html", js_source="local")

# For static PNG/SVG: open the HTML in a browser and use
# the built-in Escher export menu (Export > SVG or PNG)
```

## py4cytoscape: Automated Network Analysis

**Goal:** Build and style metabolic networks in Cytoscape programmatically.

**Approach:** Create a network from reaction-metabolite edges, apply layout algorithms, and map data to visual properties via py4cytoscape.

```python
import py4cytoscape as p4c
import pandas as pd

# Verify Cytoscape is running
p4c.cytoscape_version_info()

# Build edge table from metabolic reactions
edges = pd.DataFrame([
    {"source": "Glucose", "target": "G6P", "reaction": "HK", "flux": 8.2},
    {"source": "G6P", "target": "F6P", "reaction": "PGI", "flux": 4.9},
    {"source": "F6P", "target": "FBP", "reaction": "PFK", "flux": 7.5},
    {"source": "FBP", "target": "G3P", "reaction": "FBA", "flux": 7.5},
    {"source": "G3P", "target": "Pyruvate", "reaction": "Lower glycolysis", "flux": 15.0},
    {"source": "Pyruvate", "target": "Acetyl-CoA", "reaction": "PDH", "flux": 9.3},
    {"source": "Pyruvate", "target": "Lactate", "reaction": "LDH", "flux": 5.1},
    {"source": "Acetyl-CoA", "target": "Citrate", "reaction": "CS", "flux": 6.0},
])

# Create network in Cytoscape
p4c.create_network_from_data_frames(
    edges=edges,
    title="Glycolysis_TCA",
    collection="Metabolic Networks"
)

# Apply layout
p4c.layout_network("force-directed")

# Style: map flux to edge width
style_name = "MetabolicFlux"
p4c.create_visual_style(style_name)
p4c.set_visual_style(style_name)

p4c.set_edge_line_width_mapping(
    table_column="flux",
    table_column_values=[0, 5, 15],
    widths=[1, 4, 10],
    mapping_type="c",
    style_name=style_name
)

p4c.set_edge_color_mapping(
    table_column="flux",
    table_column_values=[0, 7.5, 15],
    colors=["#2166ac", "#f7f7f7", "#b2182b"],
    mapping_type="c",
    style_name=style_name
)

# Set node shape and color
p4c.set_node_shape_default("ELLIPSE", style_name=style_name)
p4c.set_node_color_default("#66c2a5", style_name=style_name)
p4c.set_node_label_mapping("name", style_name=style_name)

# Export
p4c.export_image("metabolic_network.png", resolution=300)
p4c.export_image("metabolic_network.svg", type="SVG")
p4c.save_session("metabolic_network.cys")
```

## KEGG Pathway Coloring

**Goal:** Color KEGG pathway diagrams by measured metabolite levels or gene expression.

**Approach:** Use the KEGG REST API to retrieve pathway images with user-specified coloring of enzymes and compounds.

```python
import requests
from urllib.parse import urlencode

def color_kegg_pathway(pathway_id, gene_colors=None, compound_colors=None):
    """Color a KEGG pathway map via the KEGG API."""
    color_entries = []
    for src in [gene_colors, compound_colors]:
        if src:
            color_entries.extend(f"{k} {v}" for k, v in src.items())

    params = {"map": pathway_id, "multi_query": "\n".join(color_entries)}
    response = requests.get("https://www.kegg.jp/kegg-bin/show_pathway", params=params)
    return response.url


# Color glycolysis pathway by fold change
gene_colors = {
    "hsa:2645": "#ff0000",   # GCK (upregulated, red)
    "hsa:5213": "#ff6666",   # PFKM (moderately up)
    "hsa:5315": "#6666ff",   # PKM (downregulated, blue)
}

compound_colors = {
    "C00031": "#ff9900",   # Glucose (elevated)
    "C00022": "#0099ff",   # Pyruvate (decreased)
}

colored_url = color_kegg_pathway("hsa00010", gene_colors, compound_colors)
print(f"Colored pathway: {colored_url}")
```

## KEGG Pathway Retrieval

**Goal:** Programmatically fetch KEGG pathway membership for metabolites and genes.

**Approach:** Query the KEGG REST API to retrieve pathway-compound and pathway-gene linkages.

```python
import requests

def get_kegg_pathways_for_compound(compound_id):
    """Get all pathways containing a KEGG compound."""
    url = f"https://rest.kegg.jp/link/pathway/{compound_id}"
    resp = requests.get(url)
    pathways = []
    for line in resp.text.strip().split("\n"):
        if line:
            parts = line.split("\t")
            pathways.append(parts[1].replace("path:", ""))
    return pathways

def get_pathway_compounds(pathway_id):
    """Get all compounds in a KEGG pathway."""
    url = f"https://rest.kegg.jp/link/compound/{pathway_id}"
    resp = requests.get(url)
    compounds = []
    for line in resp.text.strip().split("\n"):
        if line:
            parts = line.split("\t")
            compounds.append(parts[1].replace("cpd:", ""))
    return compounds


# Example: find pathways containing pyruvate
pyruvate_pathways = get_kegg_pathways_for_compound("cpd:C00022")
print(f"Pyruvate in {len(pyruvate_pathways)} pathways: {pyruvate_pathways[:5]}")

# Get compounds in glycolysis
glycolysis_compounds = get_pathway_compounds("path:hsa00010")
print(f"Glycolysis compounds: {glycolysis_compounds}")
```

## NetworkX: Custom Metabolic Networks

**Goal:** Build and analyze metabolic network topology (degree distribution, hubs, shortest paths).

**Approach:** Construct a bipartite graph of reactions and metabolites from a COBRA model, compute centrality metrics, and render with matplotlib.

```python
import networkx as nx
import cobra
import matplotlib.pyplot as plt
import numpy as np

def build_metabolic_graph(model):
    """Build a metabolite-centric directed graph from a COBRA model."""
    graph = nx.DiGraph()
    for rxn in model.reactions:
        substrates = [m.id for m in rxn.reactants]
        products = [m.id for m in rxn.products]
        for s in substrates:
            for p in products:
                graph.add_edge(s, p, reaction=rxn.id)
            if rxn.reversibility:
                for p in products:
                    graph.add_edge(p, s, reaction=rxn.id)
    return graph


model = cobra.io.load_model("textbook")
graph = build_metabolic_graph(model)

# Remove currency metabolites for cleaner visualization
currency = ["h_c", "h2o_c", "atp_c", "adp_c", "nad_c", "nadh_c",
            "pi_c", "h_e", "h2o_e", "coa_c"]
graph.remove_nodes_from([n for n in currency if n in graph])

# Centrality analysis
degree_cent = nx.degree_centrality(graph)
betweenness = nx.betweenness_centrality(graph)
top_hubs = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:10]
for met, score in top_hubs:
    print(f"  {met}: degree={score:.3f}, betweenness={betweenness.get(met, 0):.3f}")

# Visualize with flux overlay
solution = model.optimize()
flux_dict = {rxn.id: solution.fluxes[rxn.id] for rxn in model.reactions}
edge_colors = [flux_dict.get(d["reaction"], 0) for _, _, d in graph.edges(data=True)]

pos = nx.spring_layout(graph, k=2, seed=42)
node_sizes = [300 * degree_cent.get(n, 0.01) + 50 for n in graph.nodes()]

fig, ax = plt.subplots(figsize=(16, 12))
nx.draw_networkx_nodes(graph, pos, node_size=node_sizes, node_color="#66c2a5",
                        alpha=0.8, ax=ax)
nx.draw_networkx_edges(graph, pos, edge_color=edge_colors,
                        edge_cmap=plt.cm.RdBu_r, width=1.5, alpha=0.6, ax=ax)
nx.draw_networkx_labels(graph, pos, font_size=6, ax=ax)
sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r,
                            norm=plt.Normalize(vmin=min(edge_colors),
                                               vmax=max(edge_colors)))
sm.set_array([])
plt.colorbar(sm, ax=ax, label="Flux (mmol/gDW/h)")
ax.set_title("E. coli Core Metabolic Network with FBA Fluxes")
ax.axis("off")
plt.savefig("metabolic_network_flux.png", dpi=300, bbox_inches="tight")
plt.close()
```

## Differential Flux Visualization

**Goal:** Compare flux distributions between two conditions and highlight changes on the pathway map.

**Approach:** Compute flux differences (condition B minus condition A), then pass the delta dictionary to Escher for diverging color scale visualization.

```python
import escher
import cobra

model = cobra.io.load_model("textbook")

# Condition A: aerobic (default)
with model:
    sol_aerobic = model.optimize()
    flux_aerobic = {r.id: sol_aerobic.fluxes[r.id] for r in model.reactions}

# Condition B: anaerobic
with model:
    model.reactions.get_by_id("EX_o2_e").lower_bound = 0
    sol_anaerobic = model.optimize()
    flux_anaerobic = {r.id: sol_anaerobic.fluxes[r.id] for r in model.reactions}

# Compute flux change
flux_diff = {rxn_id: flux_anaerobic.get(rxn_id, 0) - flux_aerobic.get(rxn_id, 0)
             for rxn_id in set(list(flux_aerobic.keys()) + list(flux_anaerobic.keys()))}

builder = escher.Builder(
    map_name="e_coli_core.Core metabolism",
    model=model,
    reaction_data=flux_diff,
    reaction_scale=[
        {"type": "min", "color": "#2166ac", "size": 15},
        {"type": "value", "value": 0, "color": "#f7f7f7", "size": 5},
        {"type": "max", "color": "#b2182b", "size": 15},
    ],
    reaction_styles=["color", "size", "text"],
)

builder.save_html("flux_diff_aerobic_vs_anaerobic.html")
```

## Best Practices

- Remove currency metabolites (H2O, H+, ATP, NAD+) from network graphs to reduce visual clutter
- Use diverging color scales (blue-white-red) for flux data centered on zero
- For large networks (>500 nodes), use Cytoscape instead of matplotlib for interactive exploration
- Always include a color scale legend in exported figures
- Use Escher's built-in maps (available at https://escher.github.io) before creating custom maps
- For publication: export as SVG then edit in Inkscape/Illustrator for final polish
- When overlaying multiple data types, use separate visual channels (color for flux, size for confidence)

## Related Skills

- systems-biology/flux-balance-analysis - Generate flux data for overlay
- metabolomics-analysis/pathway-mapping - Map metabolites to pathways
- pathway-analysis/bioservices - Query pathway databases
