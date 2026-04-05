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

The Human Metabolome Database (HMDB) is the most comprehensive freely available database of human metabolites, containing detailed information on over 220,000 metabolite entries. Each entry includes chemical, clinical, enzymatic, and molecular biology data.

HMDB is organized into five interconnected sub-databases:

1. **Metabolites** -- small molecule metabolites found in the human body
2. **Proteins** -- enzymes and transporters linked to metabolites
3. **Pathways** -- metabolic and signaling pathways (linked to SMPDB)
4. **Diseases** -- metabolite-disease associations with biomarker information
5. **Spectra** -- experimental NMR, MS/MS, and GC-MS reference spectra

HMDB provides programmatic access through XML-based REST endpoints at `https://hmdb.ca`. Note that HMDB uses Cloudflare protection, so automated requests require a proper `User-Agent` header and may occasionally be challenged. For large-scale data access, use the bulk download files.

## When to Use This Skill

Use this skill when you need to look up metabolite identifiers, retrieve pathway and disease associations, obtain spectral references, cross-reference metabolites across databases (KEGG, PubChem, ChEBI), or perform batch queries against the human metabolome.

## Installation

```bash
uv pip install requests
```

Verify:

```python
import requests

response = requests.get(
    "https://hmdb.ca/metabolites/HMDB0000122.xml",
    headers={"User-Agent": "MetabolismSkill/1.0 (research use)"},
    timeout=30,
)
print(f"Status: {response.status_code}")
# Note: HMDB uses Cloudflare; a 403 response means the request was
# challenged.  If this happens consistently, fall back to the bulk
# XML download (see "Bulk Data Access" section below).
```

## Core Capabilities

### 1. Retrieve Metabolite Details by Accession

Each HMDB metabolite has a unique accession in the format HMDB####### (7 digits after the prefix, e.g., HMDB0000122 for D-glucose). The XML endpoint returns comprehensive data.

```python
import requests
import xml.etree.ElementTree as ET

HEADERS = {"User-Agent": "MetabolismSkill/1.0 (research use)"}

def fetch_metabolite(accession: str) -> ET.Element:
    """Fetch a single metabolite record from HMDB.

    Args:
        accession: HMDB accession ID (e.g. "HMDB0000122").
                   Must use the 7-digit format.

    Returns:
        Parsed XML root element.

    Raises:
        requests.exceptions.HTTPError: On non-2xx response (including
            403 from Cloudflare challenge).
    """
    url = f"https://hmdb.ca/metabolites/{accession}.xml"
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()
    return ET.fromstring(response.content)

# Example: look up D-glucose (HMDB0000122)
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

HMDB provides a search interface that can be queried programmatically. These endpoints return HTML, so results are parsed with regex.

```python
import re

def search_by_name(query: str) -> list:
    """Search HMDB metabolites by name.

    Returns a deduplicated list of HMDB accession IDs found on the
    search results page.
    """
    url = "https://hmdb.ca/unearth/q"
    params = {
        "query": query,
        "searcher": "metabolites",
        "button": "",
    }
    response = requests.get(url, params=params, headers=HEADERS, timeout=30)
    response.raise_for_status()
    # Returns HTML; parse accession links
    accessions = re.findall(r'href="/metabolites/(HMDB\d+)"', response.text)
    return list(dict.fromkeys(accessions))  # deduplicate, preserve order

results = search_by_name("tryptophan")
print(f"Found {len(results)} results: {results[:5]}")
```

```python
import time

def search_by_mass(mass: float, tolerance: float = 0.05) -> list:
    """Search HMDB by monoisotopic mass within a tolerance window.

    Returns a deduplicated list of HMDB accession IDs.
    """
    url = "https://hmdb.ca/spectra/ms/search"
    params = {
        "utf8": "true",
        "query_masses": str(mass),
        "tolerance": str(tolerance),
        "mode": "positive",
        "adduct_type": "M+H",
    }
    response = requests.get(url, params=params, headers=HEADERS, timeout=30)
    response.raise_for_status()
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
        List of dicts with basic metabolite info or error messages.
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

# L-Alanine, L-Histidine, Creatine
targets = ["HMDB0000161", "HMDB0000177", "HMDB0000064"]
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

### 7. Bulk Data Access (Recommended for Large-Scale Analysis)

For analysis involving more than a few dozen metabolites, download the full XML dump instead of making individual API requests. This avoids Cloudflare rate limiting and is significantly faster.

```python
import os
import zipfile

def download_hmdb_bulk(output_dir: str) -> str:
    """Download the full HMDB metabolites XML archive.

    The file is ~7 GB compressed and contains one XML record per
    metabolite.  This is the recommended approach for large-scale
    analysis.
    """
    url = "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "hmdb_metabolites.zip")

    response = requests.get(
        url,
        headers=HEADERS,
        stream=True,
        timeout=600,
    )
    response.raise_for_status()

    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    return output_path

def parse_hmdb_bulk_xml(zip_path: str):
    """Iterate over metabolite records in the bulk XML download.

    Yields one ET.Element per metabolite using iterparse to keep
    memory usage manageable.
    """
    with zipfile.ZipFile(zip_path, "r") as zf:
        xml_name = [n for n in zf.namelist() if n.endswith(".xml")][0]
        with zf.open(xml_name) as xml_file:
            ns = "{http://www.hmdb.ca}"
            for event, elem in ET.iterparse(xml_file, events=("end",)):
                if elem.tag == f"{ns}metabolite":
                    yield elem
                    elem.clear()
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
       time.sleep(1.0)
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

1. **Include a User-Agent header**: HMDB is behind Cloudflare. Always send a descriptive `User-Agent` header (e.g. `"MetabolismSkill/1.0 (research use)"`) to reduce the chance of being blocked.

2. **Rate limiting**: Always include at least a 1-second delay between requests. HMDB does not publish rate limits, but aggressive querying can result in IP blocks.

3. **Use bulk downloads for large-scale work**: For more than ~50 metabolites, download the full XML dump from `https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip` rather than hitting the API per-metabolite.

4. **Use the namespace**: HMDB XML uses the `http://www.hmdb.ca` namespace. All `findtext` and `findall` calls must use the `{http://www.hmdb.ca}` prefix.

5. **Handle missing fields**: Not all metabolites have complete records. Always check for `None` return values from `findtext`.

6. **Cache responses locally**: For batch analysis, save XML files to disk to avoid repeated downloads.

7. **Prefer accession lookups over search**: Direct accession queries (`/metabolites/HMDBXXXXXXX.xml`) are the most reliable and fastest endpoint.

8. **Normalize accession format**: HMDB accessions must be zero-padded to 7 digits after the prefix (e.g., `HMDB0000122`, not `HMDB00122`). Older literature may use a 5-digit format; always convert to 7 digits.

## Resources

- **HMDB home**: https://hmdb.ca
- **HMDB downloads**: https://hmdb.ca/downloads (full XML dumps, SDF files, CSV)
- **HMDB citing/about**: https://hmdb.ca/about#citing
- **SMPDB (linked pathways)**: https://smpdb.ca
