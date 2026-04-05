---
name: metabolights-database
description: Access EMBL-EBI MetaboLights via REST API. Search and retrieve metabolomics study metadata, ISA-Tab files, assay data, and study files programmatically. Supports organism, technology, and metabolite-based queries across 3,000+ public studies.
tool_type: python
primary_tool: requests
metadata:
    skill-author: Fei Wen
---

# MetaboLights Database

## Overview

MetaboLights is the primary open-access repository for metabolomics studies, hosted by EMBL-EBI. It stores over 3,000 public studies with raw and processed data in the ISA-Tab metadata standard. The REST API provides programmatic access to study metadata, sample information, assay details, and downloadable data files.

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

Each MetaboLights study has an accession in the format MTBLS followed by a number.

```python
import requests

BASE_URL = "https://www.ebi.ac.uk/metabolights/ws"

def get_study(accession: str) -> dict:
    """Retrieve study metadata from MetaboLights."""
    url = f"{BASE_URL}/studies/{accession}"
    response = requests.get(
        url,
        headers={"Accept": "application/json"},
        timeout=30,
    )
    response.raise_for_status()
    return response.json()

# Example: retrieve MTBLS1 (the first MetaboLights study)
study = get_study("MTBLS1")
print(f"Title: {study.get('title')}")
print(f"Description: {study.get('description', '')[:200]}")
print(f"Submission date: {study.get('submissionDate')}")
print(f"Status: {study.get('studyStatus')}")
```

### 2. List Study Contacts and Protocols

```python
def get_study_contacts(accession: str) -> list:
    """Retrieve contact information for a study."""
    url = f"{BASE_URL}/studies/{accession}/contacts"
    response = requests.get(
        url,
        headers={"Accept": "application/json"},
        timeout=30,
    )
    response.raise_for_status()
    return response.json().get("contacts", [])

def get_study_protocols(accession: str) -> list:
    """Retrieve protocols used in a study."""
    url = f"{BASE_URL}/studies/{accession}/protocols"
    response = requests.get(
        url,
        headers={"Accept": "application/json"},
        timeout=30,
    )
    response.raise_for_status()
    return response.json().get("protocols", [])

contacts = get_study_contacts("MTBLS1")
for c in contacts:
    print(f"Contact: {c.get('firstName')} {c.get('lastName')} ({c.get('email', 'N/A')})")

protocols = get_study_protocols("MTBLS1")
for p in protocols:
    print(f"Protocol: {p.get('name')}")
```

### 3. Access Assay and Sample Tables

MetaboLights uses ISA-Tab format with three table types: investigation (i_), study (s_), and assay (a_) files.

```python
import pandas as pd
from io import StringIO

def get_study_files_list(accession: str) -> list:
    """List all files in a MetaboLights study."""
    url = f"{BASE_URL}/studies/{accession}/files"
    response = requests.get(
        url,
        headers={"Accept": "application/json"},
        timeout=60,
    )
    response.raise_for_status()
    return response.json().get("study", [])

def download_isa_table(accession: str, filename: str) -> pd.DataFrame:
    """Download and parse an ISA-Tab table file (s_*.txt or a_*.txt)."""
    url = (
        f"https://www.ebi.ac.uk/metabolights/{accession}/files/"
        f"{filename}"
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return pd.read_csv(StringIO(response.text), sep="\t")

# List files in MTBLS1
files = get_study_files_list("MTBLS1")
for f in files:
    file_name = f.get("file", "")
    if file_name.startswith(("s_", "a_", "i_", "m_")):
        print(f"ISA file: {file_name}")
```

### 4. Parse ISA-Tab Investigation File

The investigation file contains top-level study design information.

```python
def parse_investigation_file(accession: str) -> dict:
    """Download and parse the i_Investigation.txt file."""
    url = (
        f"https://www.ebi.ac.uk/metabolights/{accession}/files/"
        f"i_Investigation.txt"
    )
    response = requests.get(url, timeout=60)
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

### 5. Search Studies by Organism, Technology, or Metabolite

```python
def search_studies(query: str) -> list:
    """Search MetaboLights studies by keyword."""
    url = f"{BASE_URL}/search/studies"
    params = {"query": query}
    response = requests.get(
        url,
        params=params,
        headers={"Accept": "application/json"},
        timeout=30,
    )
    response.raise_for_status()
    return response.json().get("content", [])

# Search for human plasma studies
results = search_studies("human plasma")
for study in results[:5]:
    print(f"{study.get('accession')}: {study.get('title', '')[:80]}")
```

```python
def get_public_studies() -> list:
    """Retrieve a list of all public MetaboLights study accessions."""
    url = f"{BASE_URL}/studies"
    response = requests.get(
        url,
        headers={"Accept": "application/json"},
        timeout=60,
    )
    response.raise_for_status()
    return response.json().get("content", [])

public_studies = get_public_studies()
print(f"Total public studies: {len(public_studies)}")
```

### 6. Download Study Data Files

```python
import os

def download_study_file(accession: str, filename: str, output_dir: str) -> str:
    """Download a specific file from a MetaboLights study."""
    url = (
        f"https://www.ebi.ac.uk/metabolights/{accession}/files/"
        f"{filename}"
    )
    response = requests.get(url, stream=True, timeout=120)
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

## Common Workflows

### Workflow 1: Explore a Study End-to-End

1. Get study metadata:
   ```python
   study = get_study("MTBLS1")
   print(f"Title: {study.get('title')}")
   print(f"Organism: {study.get('organism', [])}")
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

4. Download and inspect the assay table:
   ```python
   assay_df = download_isa_table("MTBLS1", "a_MTBLS1_metabolite_profiling_mass_spectrometry.txt")
   print(f"Assays: {len(assay_df)}")
   ```

### Workflow 2: Find Studies for a Specific Organism

1. Search for studies:
   ```python
   results = search_studies("Arabidopsis thaliana")
   accessions = [s.get("accession") for s in results]
   print(f"Found {len(accessions)} studies")
   ```

2. Retrieve metadata for each:
   ```python
   import time
   for acc in accessions[:5]:
       study = get_study(acc)
       print(f"{acc}: {study.get('title', '')[:80]}")
       time.sleep(0.5)
   ```

3. Collect sample counts:
   ```python
   for acc in accessions[:5]:
       try:
           files = get_study_files_list(acc)
           sample_files = [
               f.get("file") for f in files
               if f.get("file", "").startswith("s_")
           ]
           print(f"{acc}: {len(sample_files)} sample table(s)")
       except requests.exceptions.RequestException as e:
           print(f"{acc}: error - {e}")
       time.sleep(0.5)
   ```

### Workflow 3: Compare Sample Metadata Across Studies

1. Define studies of interest:
   ```python
   study_ids = ["MTBLS1", "MTBLS2", "MTBLS3"]
   ```

2. Download sample tables:
   ```python
   import time
   sample_tables = {}
   for sid in study_ids:
       try:
           files = get_study_files_list(sid)
           s_file = next(
               (f.get("file") for f in files if f.get("file", "").startswith("s_")),
               None,
           )
           if s_file:
               sample_tables[sid] = download_isa_table(sid, s_file)
       except Exception as e:
           print(f"Skipping {sid}: {e}")
       time.sleep(0.5)
   ```

3. Compare column schemas:
   ```python
   for sid, df in sample_tables.items():
       print(f"\n{sid} columns ({len(df)} samples):")
       print(f"  {df.columns.tolist()[:8]}")
   ```

## Best Practices

1. **Use Accept headers**: Always include `"Accept": "application/json"` for JSON responses. Without it, some endpoints return HTML.

2. **Rate limiting**: Add 0.5-1 second delays between API calls to be respectful of EMBL-EBI servers.

3. **Handle large files carefully**: Some studies contain hundreds of raw data files (mzML, RAW). Use streaming downloads (`stream=True`) and check file sizes before downloading.

4. **ISA-Tab column names vary**: Column names in study and assay tables are not standardized across studies. Always inspect `df.columns` before accessing specific fields.

5. **Check study status**: Not all studies are public. Use the `studyStatus` field to verify a study is available before downloading files.

6. **Cache metadata locally**: For repeated queries, cache study metadata and file lists to reduce API calls.

7. **Use the FTP mirror for bulk downloads**: For downloading entire studies, the EMBL-EBI FTP server (`ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/`) is more appropriate than the REST API.

## Resources

- **MetaboLights home**: https://www.ebi.ac.uk/metabolights/
- **REST API docs**: https://www.ebi.ac.uk/metabolights/ws/api-docs
- **ISA-Tab format specification**: https://isa-specs.readthedocs.io/
- **FTP bulk access**: ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/
- **MetaboLights GitHub**: https://github.com/EBIBioStudies/metabolern
