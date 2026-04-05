---
name: metabolights-database
description: Access EMBL-EBI MetaboLights via REST API. Search and retrieve metabolomics study metadata, ISA-Tab files, assay data, and study files programmatically. Supports organism, technology, and metabolite-based queries across 2,800+ public studies.
tool_type: python
primary_tool: requests
metadata:
    skill-author: Fei Wen
---

# MetaboLights Database

## Overview

MetaboLights is the primary open-access repository for metabolomics studies, hosted by EMBL-EBI. It stores over 2,800 public studies with raw and processed data in the ISA-Tab metadata standard. The REST API provides programmatic access to study metadata, sample information, assay details, and downloadable data files.

MetaboLights is built on the **ISA (Investigation-Study-Assay) framework**, a community standard for describing multi-omics experiments. Each study contains:

- **i_Investigation.txt** -- top-level study design, contacts, publications, ontology references
- **s_*.txt** -- sample tables with organism, tissue, and factor information
- **a_*.txt** -- assay tables linking samples to data files and analytical methods
- **m_*.tsv** -- MAF (Metabolite Assignment File) tables with identified metabolites, chemical identifiers, and quantification data

The MAF format is MetaboLights-specific and provides curated metabolite annotations with columns for database identifiers (HMDB, ChEBI, KEGG, PubChem), chemical formula, SMILES, and InChI.

## When to Use This Skill

Use this skill when you need to search for public metabolomics studies by organism or technology, retrieve ISA-Tab metadata (investigation, study, assay tables), download study data files, or access curated metabolite annotations from deposited studies.

## Installation

```bash
uv pip install requests pandas
```

Verify:

```python
import requests

response = requests.get(
    "https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS1",
    headers={"Accept": "application/json"},
    timeout=30,
)
print(f"Status: {response.status_code}")
```

## Core Capabilities

### 1. Retrieve Study Metadata by Accession

Each MetaboLights study has an accession in the format MTBLS followed by a number. The API response contains two main sections: `mtblsStudy` (status, permissions, URLs) and `isaInvestigation` (scientific metadata following the ISA framework).

```python
import requests

BASE_URL = "https://www.ebi.ac.uk/metabolights/ws"
HEADERS = {"Accept": "application/json"}

def get_study(accession: str) -> dict:
    """Retrieve study metadata from MetaboLights.

    Returns a dict with keys: mtblsStudy, isaInvestigation, validation.
    """
    url = f"{BASE_URL}/studies/{accession}"
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()
    return response.json()

# Example: retrieve MTBLS1 (the first MetaboLights study)
data = get_study("MTBLS1")
isa = data["isaInvestigation"]
mtbls = data["mtblsStudy"]

print(f"Title: {isa.get('title')}")
print(f"Submission date: {isa.get('submissionDate')}")
print(f"Status: {mtbls.get('studyStatus')}")

# Study-level details are nested inside isaInvestigation > studies[]
if isa.get("studies"):
    study_detail = isa["studies"][0]
    print(f"Description: {study_detail.get('description', '')[:200]}")
    print(f"Factors: {[f.get('factorName') for f in study_detail.get('factors', [])]}")
```

### 2. List Study Contacts and Protocols

```python
def get_study_contacts(accession: str) -> list:
    """Retrieve contact information for a study."""
    url = f"{BASE_URL}/studies/{accession}/contacts"
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()
    return response.json().get("contacts", [])

def get_study_protocols(accession: str) -> list:
    """Retrieve protocols used in a study."""
    url = f"{BASE_URL}/studies/{accession}/protocols"
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()
    return response.json().get("protocols", [])

contacts = get_study_contacts("MTBLS1")
for c in contacts:
    print(f"Contact: {c.get('firstName')} {c.get('lastName')} ({c.get('email', 'N/A')})")

protocols = get_study_protocols("MTBLS1")
for p in protocols:
    print(f"Protocol: {p.get('name')}")
```

### 3. Access Assay, Sample, and MAF Tables

MetaboLights uses ISA-Tab format with specific file naming conventions:

- `i_Investigation.txt` -- investigation metadata (one per study)
- `s_<study>.txt` -- sample table
- `a_<study>_<assay_type>.txt` -- assay table
- `m_<study>_<assay_type>_maf.tsv` -- Metabolite Assignment File (MAF)

```python
import pandas as pd
from io import StringIO

def get_study_files_list(accession: str) -> list:
    """List all files in a MetaboLights study.

    Returns a list of dicts with keys: file, type, status, directory,
    createdAt, timestamp.
    """
    url = f"{BASE_URL}/studies/{accession}/files"
    response = requests.get(url, headers=HEADERS, timeout=60)
    response.raise_for_status()
    return response.json().get("study", [])

def download_isa_table(accession: str, filename: str) -> pd.DataFrame:
    """Download and parse an ISA-Tab table file (s_*.txt, a_*.txt, or m_*.tsv).

    Uses the /download?file= endpoint which returns raw file content.
    """
    url = f"{BASE_URL}/studies/{accession}/download"
    params = {"file": filename}
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()
    return pd.read_csv(StringIO(response.text), sep="\t")

# List files in MTBLS1
files = get_study_files_list("MTBLS1")
for f in files:
    file_name = f.get("file", "")
    file_type = f.get("type", "")
    if file_name.startswith(("s_", "a_", "i_", "m_")):
        print(f"ISA file: {file_name}  (type: {file_type})")
```

### 4. Parse ISA-Tab Investigation File

The investigation file contains top-level study design information including ontology references, study design descriptors, and publication details.

```python
def parse_investigation_file(accession: str) -> dict:
    """Download and parse the i_Investigation.txt file.

    Returns a dict mapping section names to lists of (key, value) tuples.
    Sections include ONTOLOGY SOURCE REFERENCE, INVESTIGATION,
    STUDY, STUDY DESIGN DESCRIPTORS, STUDY PROTOCOLS, etc.
    """
    url = f"{BASE_URL}/studies/{accession}/download"
    params = {"file": "i_Investigation.txt"}
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()

    sections = {}
    current_section = None
    for line in response.text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current_section = line.strip("[]")
            sections[current_section] = []
        elif current_section is not None:
            parts = line.split("\t", 1)
            key = parts[0].strip('"')
            value = parts[1].strip('"') if len(parts) > 1 else ""
            sections[current_section].append((key, value))
    return sections

investigation = parse_investigation_file("MTBLS1")
for section, fields in investigation.items():
    print(f"\n[{section}]")
    for key, value in fields[:3]:
        print(f"  {key}: {value}")
```

### 5. Get Assay Filenames

```python
def get_study_assays(accession: str) -> list:
    """Retrieve assay filenames for a study.

    Returns a list of dicts with at least a 'filename' key.
    """
    url = f"{BASE_URL}/studies/{accession}/assays"
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()
    data = response.json()
    return data.get("data", {}).get("assays", [])

assays = get_study_assays("MTBLS1")
for a in assays:
    print(f"Assay file: {a.get('filename')}")
```

### 6. List Public Studies

The `/studies` endpoint returns all public study accessions and a total count.

```python
def get_public_studies() -> dict:
    """Retrieve public study accessions and total count.

    Returns a dict with 'content' (list of accession strings) and
    'studies' (total count integer).
    """
    url = f"{BASE_URL}/studies"
    response = requests.get(url, headers=HEADERS, timeout=60)
    response.raise_for_status()
    return response.json()

result = get_public_studies()
print(f"Total public studies: {result.get('studies')}")
```

### 7. Download Study Data Files

For individual file downloads, use the `/download?file=` endpoint. For downloading entire studies, use the FTP/HTTP mirror.

```python
import os

def download_study_file(accession: str, filename: str, output_dir: str) -> str:
    """Download a specific file from a MetaboLights study.

    Uses the REST API download endpoint for metadata files. For large
    raw data files (mzML, RAW), prefer the FTP/HTTP mirror instead.
    """
    url = f"{BASE_URL}/studies/{accession}/download"
    params = {"file": filename}
    response = requests.get(url, params=params, stream=True, timeout=120)
    response.raise_for_status()

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    return output_path

# Download the sample table from MTBLS1
path = download_study_file("MTBLS1", "s_MTBLS1.txt", "/tmp/mtbls1")
print(f"Downloaded to: {path}")
```

### 8. Download via FTP/HTTP Mirror (Bulk Access)

For downloading entire studies or large raw data files, use the EBI public HTTP/FTP mirror.

```python
def download_from_ftp_mirror(
    accession: str,
    filename: str,
    output_dir: str,
) -> str:
    """Download a file from the EBI public HTTP mirror of MetaboLights.

    This is faster and more reliable for large files than the REST API.
    """
    url = (
        f"https://ftp.ebi.ac.uk/pub/databases/metabolights/"
        f"studies/public/{accession}/{filename}"
    )
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    return output_path

# Example: download the sample table via HTTP mirror
path = download_from_ftp_mirror("MTBLS1", "s_MTBLS1.txt", "/tmp/mtbls1")
print(f"Downloaded to: {path}")
```

## Common Workflows

### Workflow 1: Explore a Study End-to-End

1. Get study metadata:
   ```python
   data = get_study("MTBLS1")
   isa = data["isaInvestigation"]
   print(f"Title: {isa.get('title')}")
   ```

2. List available files:
   ```python
   files = get_study_files_list("MTBLS1")
   isa_files = [
       f.get("file") for f in files
       if f.get("file", "").startswith(("s_", "a_", "m_"))
   ]
   print(f"ISA tables: {isa_files}")
   ```

3. Download and inspect the sample table:
   ```python
   sample_df = download_isa_table("MTBLS1", "s_MTBLS1.txt")
   print(f"Samples: {len(sample_df)}")
   print(sample_df.columns.tolist())
   ```

4. Download and inspect the MAF (metabolite annotations):
   ```python
   maf_df = download_isa_table(
       "MTBLS1",
       "m_MTBLS1_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv",
   )
   print(f"Identified metabolites: {len(maf_df)}")
   print(maf_df.columns.tolist())
   # Typical MAF columns: database_identifier, chemical_formula,
   # smiles, inchi, metabolite_identification, etc.
   ```

### Workflow 2: Find Studies for a Specific Organism

1. Get study descriptors to check for organism information:
   ```python
   import time

   data = get_study("MTBLS1")
   study_detail = data["isaInvestigation"]["studies"][0]
   descriptors = study_detail.get("studyDesignDescriptors", [])
   for d in descriptors:
       print(f"Descriptor: {d.get('annotationValue')}")
   ```

2. Retrieve and filter multiple studies:
   ```python
   # Get a set of known accessions and check each
   accessions = ["MTBLS1", "MTBLS2", "MTBLS3"]
   for acc in accessions:
       try:
           data = get_study(acc)
           isa = data["isaInvestigation"]
           print(f"{acc}: {isa.get('title', '')[:80]}")
       except requests.exceptions.RequestException as e:
           print(f"{acc}: error - {e}")
       time.sleep(0.5)
   ```

### Workflow 3: Extract Metabolite Annotations from MAF Files

1. Get the MAF file for a study:
   ```python
   files = get_study_files_list("MTBLS1")
   maf_files = [
       f.get("file") for f in files
       if f.get("file", "").startswith("m_") and f.get("file", "").endswith(".tsv")
   ]
   print(f"MAF files: {maf_files}")
   ```

2. Download and inspect metabolite identifiers:
   ```python
   if maf_files:
       maf_df = download_isa_table("MTBLS1", maf_files[0])
       id_col = "database_identifier"
       if id_col in maf_df.columns:
           identifiers = maf_df[id_col].dropna().unique()
           print(f"Unique metabolite IDs: {len(identifiers)}")
           print(f"Sample IDs: {list(identifiers[:5])}")
   ```

## Best Practices

1. **Use Accept headers**: Always include `"Accept": "application/json"` for JSON responses. Without it, some endpoints return HTML.

2. **Rate limiting**: Add 0.5-1 second delays between API calls to be respectful of EMBL-EBI servers.

3. **Handle large files carefully**: Some studies contain hundreds of raw data files (mzML, RAW). Use streaming downloads (`stream=True`) and check file sizes before downloading.

4. **ISA-Tab column names vary**: Column names in study and assay tables are not standardized across studies. Always inspect `df.columns` before accessing specific fields.

5. **Check study status**: Not all studies are public. Use the `studyStatus` field in `mtblsStudy` to verify a study is available before downloading files.

6. **Cache metadata locally**: For repeated queries, cache study metadata and file lists to reduce API calls.

7. **Use the HTTP/FTP mirror for bulk downloads**: For downloading entire studies or large raw files, the EBI HTTP mirror (`https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/`) is more appropriate and reliable than the REST API.

8. **Use the correct download endpoint**: File downloads use `GET /studies/{accession}/download?file={filename}`, not a URL path with the filename appended.

9. **Understand the ISA-MAF relationship**: The MAF file is the key output for metabolite identification results. It links assay measurements to curated chemical identifiers (HMDB, ChEBI, KEGG, PubChem) and should be your starting point for extracting metabolite lists from a study.

## Resources

- **MetaboLights home**: https://www.ebi.ac.uk/metabolights/
- **REST API docs**: https://www.ebi.ac.uk/metabolights/ws/api-docs
- **ISA-Tab format specification**: https://isa-specs.readthedocs.io/
- **HTTP bulk access**: https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/
- **MetaboLights GitHub**: https://github.com/EBI-Metabolights/MetaboLights
