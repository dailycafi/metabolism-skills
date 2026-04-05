---
name: hmdb-database
description: Access the Human Metabolome Database (HMDB) via REST API (220,000+ metabolites). Retrieve metabolite details, pathways, disease associations, spectra references, and cross-references to KEGG, PubChem, and ChEBI.
tool_type: python
primary_tool: requests
metadata:
    skill-author: Fei Wen
---

# HMDB Database Lookup

## Overview

The Human Metabolome Database (HMDB) is the most comprehensive freely available database of human metabolites, containing detailed information on over 220,000 metabolite entries. Each entry includes chemical, clinical, enzymatic, and molecular biology data. HMDB provides programmatic access through XML-based REST endpoints and a CSV search API for name, mass, and formula queries.

## When to Use This Skill

Use this skill when you need to look up metabolite identifiers, retrieve pathway and disease associations, obtain spectral references, cross-reference metabolites across databases (KEGG, PubChem, ChEBI), or perform batch queries against the human metabolome.

## Installation

```bash
uv pip install requests
```

Verify:

```python
import requests
response = requests.get("https://hmdb.ca/metabolites/HMDB0000122.xml", timeout=30)
print(f"Status: {response.status_code}")
```

## Core Capabilities

### 1. Retrieve Metabolite Details by Accession

Each HMDB metabolite has a unique accession (HMDB#######) and returns XML with comprehensive data.

```python
import requests
import xml.etree.ElementTree as ET

def fetch_metabolite(accession: str) -> ET.Element:
    """Fetch a single metabolite record from HMDB."""
    url = f"https://hmdb.ca/metabolites/{accession}.xml"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    return ET.fromstring(response.content)

# Example: look up glucose (HMDB0000122)
root = fetch_metabolite("HMDB0000122")

ns = "{http://www.hmdb.ca}"
name = root.findtext(f"{ns}name")
formula = root.findtext(f"{ns}chemical_formula")
mono_mass = root.findtext(f"{ns}monisotopic_molecular_weight")
avg_mass = root.findtext(f"{ns}average_molecular_weight")
iupac = root.findtext(f"{ns}iupac_name")
smiles = root.findtext(f"{ns}smiles")
inchikey = root.findtext(f"{ns}inchikey")

print(f"Name: {name}")
print(f"Formula: {formula}")
print(f"Monoisotopic mass: {mono_mass}")
print(f"Average mass: {avg_mass}")
print(f"IUPAC: {iupac}")
print(f"SMILES: {smiles}")
print(f"InChIKey: {inchikey}")
```

### 2. Extract Cross-Reference Identifiers

Pull KEGG, PubChem, ChEBI, and other database IDs from a metabolite record.

```python
def extract_cross_references(root: ET.Element) -> dict:
    """Extract external database identifiers from an HMDB XML record."""
    ns = "{http://www.hmdb.ca}"
    return {
        "hmdb_id": root.findtext(f"{ns}accession"),
        "kegg_id": root.findtext(f"{ns}kegg_id"),
        "pubchem_cid": root.findtext(f"{ns}pubchem_compound_id"),
        "chebi_id": root.findtext(f"{ns}chebi_id"),
        "drugbank_id": root.findtext(f"{ns}drugbank_id"),
        "chemspider_id": root.findtext(f"{ns}chemspider_id"),
        "metlin_id": root.findtext(f"{ns}metlin_id"),
        "cas_registry": root.findtext(f"{ns}cas_registry_number"),
    }

root = fetch_metabolite("HMDB0000122")
xrefs = extract_cross_references(root)
for db, identifier in xrefs.items():
    print(f"{db}: {identifier}")
```

### 3. Retrieve Pathway and Disease Associations

```python
def extract_pathways(root: ET.Element) -> list:
    """Extract biological pathways linked to a metabolite."""
    ns = "{http://www.hmdb.ca}"
    pathways = []
    for pw in root.findall(f".//{ns}pathway"):
        pathways.append({
            "name": pw.findtext(f"{ns}name"),
            "smpdb_id": pw.findtext(f"{ns}smpdb_id"),
            "kegg_map_id": pw.findtext(f"{ns}kegg_map_id"),
        })
    return pathways

def extract_diseases(root: ET.Element) -> list:
    """Extract disease associations for a metabolite."""
    ns = "{http://www.hmdb.ca}"
    diseases = []
    for disease in root.findall(f".//{ns}disease"):
        diseases.append({
            "name": disease.findtext(f"{ns}name"),
            "omim_id": disease.findtext(f"{ns}omim_id"),
        })
    return diseases

root = fetch_metabolite("HMDB0000122")
for pw in extract_pathways(root):
    print(f"Pathway: {pw['name']} (SMPDB: {pw['smpdb_id']})")

for disease in extract_diseases(root):
    print(f"Disease: {disease['name']} (OMIM: {disease['omim_id']})")
```

### 4. Search by Name, Mass, or Formula

HMDB provides a search interface that can be queried programmatically.

```python
def search_by_name(query: str) -> list:
    """Search HMDB metabolites by name."""
    url = "https://hmdb.ca/unearth/q"
    params = {
        "query": query,
        "searcher": "metabolites",
        "button": "",
    }
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    # Returns HTML; parse accession links
    import re
    accessions = re.findall(r'href="/metabolites/(HMDB\d+)"', response.text)
    return list(dict.fromkeys(accessions))  # deduplicate, preserve order

results = search_by_name("tryptophan")
print(f"Found {len(results)} results: {results[:5]}")
```

```python
import time

def search_by_mass(mass: float, tolerance: float = 0.05) -> list:
    """Search HMDB by monoisotopic mass within a tolerance window."""
    url = "https://hmdb.ca/spectra/ms/search"
    params = {
        "utf8": "true",
        "query_masses": str(mass),
        "tolerance": str(tolerance),
        "mode": "positive",
        "adduct_type": "M+H",
    }
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    import re
    accessions = re.findall(r'href="/metabolites/(HMDB\d+)"', response.text)
    return list(dict.fromkeys(accessions))

# Search for compounds near glucose monoisotopic mass (180.0634 Da)
hits = search_by_mass(180.0634, tolerance=0.05)
print(f"Mass search hits: {hits[:5]}")
```

### 5. Batch Queries with Rate Limiting

```python
import time

def batch_fetch_metabolites(accessions: list, delay: float = 1.0) -> list:
    """Fetch multiple metabolite records with rate limiting.

    Args:
        accessions: List of HMDB accession IDs.
        delay: Seconds to wait between requests (respect server limits).

    Returns:
        List of parsed XML root elements.
    """
    results = []
    for i, accession in enumerate(accessions):
        try:
            root = fetch_metabolite(accession)
            ns = "{http://www.hmdb.ca}"
            results.append({
                "accession": accession,
                "name": root.findtext(f"{ns}name"),
                "formula": root.findtext(f"{ns}chemical_formula"),
                "mass": root.findtext(f"{ns}monisotopic_molecular_weight"),
            })
        except requests.exceptions.RequestException as e:
            results.append({"accession": accession, "error": str(e)})
        if i < len(accessions) - 1:
            time.sleep(delay)
    return results

targets = ["HMDB0000122", "HMDB0000161", "HMDB0000064", "HMDB0000148"]
records = batch_fetch_metabolites(targets, delay=1.5)
for rec in records:
    print(rec)
```

### 6. Extract Spectra References

```python
def extract_spectra_references(root: ET.Element) -> dict:
    """Extract NMR and MS spectra references from an HMDB record."""
    ns = "{http://www.hmdb.ca}"
    spectra = {"nmr": [], "ms": []}

    for spec in root.findall(f".//{ns}spectrum"):
        spec_type = spec.findtext(f"{ns}type")
        entry = {
            "type": spec_type,
            "spectrum_id": spec.findtext(f"{ns}spectrum_id"),
        }
        if spec_type and "NMR" in spec_type.upper():
            spectra["nmr"].append(entry)
        elif spec_type and "MS" in spec_type.upper():
            spectra["ms"].append(entry)
    return spectra

root = fetch_metabolite("HMDB0000122")
spectra = extract_spectra_references(root)
print(f"NMR spectra: {len(spectra['nmr'])}")
print(f"MS spectra: {len(spectra['ms'])}")
```

## Common Workflows

### Workflow 1: Identify Unknown Metabolite from MS Data

1. Search by observed monoisotopic mass:
   ```python
   candidates = search_by_mass(180.0634, tolerance=0.02)
   ```

2. Retrieve full records for top candidates:
   ```python
   records = batch_fetch_metabolites(candidates[:3], delay=1.0)
   ```

3. Cross-reference with KEGG and PubChem for validation:
   ```python
   for acc in candidates[:3]:
       root = fetch_metabolite(acc)
       xrefs = extract_cross_references(root)
       print(f"{acc}: KEGG={xrefs['kegg_id']}, PubChem={xrefs['pubchem_cid']}")
   ```

### Workflow 2: Build a Pathway-Metabolite Map

1. Start with a list of metabolites of interest:
   ```python
   metabolites = ["HMDB0000122", "HMDB0000161", "HMDB0000190"]
   ```

2. Extract pathways for each:
   ```python
   pathway_map = {}
   for acc in metabolites:
       root = fetch_metabolite(acc)
       ns = "{http://www.hmdb.ca}"
       name = root.findtext(f"{ns}name")
       pathways = extract_pathways(root)
       pathway_map[name] = [pw["name"] for pw in pathways]
       time.sleep(1.0)
   ```

3. Find shared pathways:
   ```python
   from collections import Counter
   all_pathways = [pw for pws in pathway_map.values() for pw in pws]
   shared = [pw for pw, count in Counter(all_pathways).items() if count > 1]
   print(f"Shared pathways: {shared}")
   ```

## Best Practices

1. **Rate limiting**: Always include at least a 1-second delay between requests. HMDB does not publish rate limits, but aggressive querying can result in IP blocks.

2. **Use the namespace**: HMDB XML uses the `http://www.hmdb.ca` namespace. All `findtext` and `findall` calls must use the `{http://www.hmdb.ca}` prefix.

3. **Handle missing fields**: Not all metabolites have complete records. Always check for `None` return values from `findtext`.

4. **Cache responses locally**: For batch analysis, save XML files to disk to avoid repeated downloads.

5. **Prefer accession lookups over search**: Direct accession queries (`/metabolites/HMDBXXXXXXX.xml`) are the most reliable and fastest endpoint.

6. **Normalize accession format**: HMDB accessions should be zero-padded to 7 digits after the prefix (e.g., `HMDB0000122`, not `HMDB00122`).

## Resources

- **HMDB home**: https://hmdb.ca
- **HMDB downloads**: https://hmdb.ca/downloads (full XML dumps, SDF files)
- **HMDB API documentation**: https://hmdb.ca/about#citing
- **SMPDB (linked pathways)**: https://smpdb.ca
