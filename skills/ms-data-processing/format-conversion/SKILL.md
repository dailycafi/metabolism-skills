---
name: ms-format-conversion
description: "Convert mass spectrometry vendor files to open formats using msConvert and pyopenms. Use when: user needs to convert RAW files to mzML, convert WIFF/Agilent .d/Waters .raw to open format, batch convert instrument files, or troubleshoot format issues. Triggers: convert RAW file, mzML conversion, mzXML, vendor format, msConvert, ProteoWizard, Thermo RAW, AB SCIEX WIFF, Agilent .d, Waters .raw, Bruker .d, centroid mode, profile mode, batch conversion, format conversion."
tool_type: cli
primary_tool: msconvert
metadata:
    skill-author: Fei Wen
---

# Metabolomics Data Format Conversion

## Overview

Mass spectrometry instruments produce data in proprietary vendor formats that must be converted to open standards (mzML, mzXML) for downstream analysis. This skill covers msConvert (ProteoWizard) for vendor-to-open conversion, pyopenms for mzML/mzXML reading and interconversion, batch processing scripts, and common pitfalls around centroiding, encoding precision, and vendor-specific quirks.

## When to Use This Skill

Use this skill when converting raw instrument files to mzML or mzXML, interconverting between open formats, performing batch conversion across many files, or troubleshooting format-related issues in metabolomics pipelines.

## Vendor Format Reference

| Vendor | Format | Extension | Platform | Notes |
|--------|--------|-----------|----------|-------|
| Thermo Fisher | RAW | `.raw` | Windows | Most common in proteomics/metabolomics |
| AB SCIEX | WIFF | `.wiff`, `.wiff2` | Windows | `.wiff` requires `.wiff.scan` sidecar; `.wiff2` requires `.wiff2.scan` |
| Agilent | MassHunter | `.d` (directory) | Windows | Directory containing `AcqData/` |
| Waters | MassLynx | `.raw` (directory) | Windows | **Directory**, not a single file (unlike Thermo `.raw`) |
| Bruker | BAF/TDF | `.d` (directory) | Windows/Linux | Contains `ser`/`fid` files; timsTOF uses `.tdf`/`.tdf_bin` inside `.d` |

## Installation

### msConvert (ProteoWizard)

msConvert is the standard tool for vendor format conversion. Vendor format reading requires Windows or Windows libraries via Wine/Docker.

**Windows (native):**

Download from https://proteowizard.sourceforge.io/download.html and install. After installation:

```bash
# Verify installation (Windows cmd or PowerShell)
msconvert --help
```

**Linux/macOS (Docker):**

```bash
# Pull the ProteoWizard Docker image with Wine for vendor support
docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest

# Verify
docker run --rm chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest wine msconvert --help
```

### pyopenms (for open format operations)

```bash
uv pip install pyopenms
```

## Core Capabilities

### 1. Basic Vendor-to-mzML Conversion with msConvert

```bash
# Convert a single Thermo RAW file to mzML
msconvert sample.raw --mzML --outdir /data/converted/

# Convert with centroiding (recommended for most workflows)
msconvert sample.raw --mzML --filter "peakPicking vendor msLevel=1-2" --outdir /data/converted/

# Convert to mzXML instead
msconvert sample.raw --mzXML --filter "peakPicking vendor msLevel=1-2" --outdir /data/converted/
```

### 2. Centroid vs Profile Mode

Profile mode stores the full peak shape; centroid mode stores only peak centroids (m/z, intensity pairs). Most analysis tools expect centroid data.

```bash
# Vendor-native centroiding (highest quality, requires vendor libraries)
msconvert sample.raw --mzML \
    --filter "peakPicking vendor msLevel=1-2"

# Algorithm-based centroiding (works without vendor libraries)
msconvert sample.mzML --mzML \
    --filter "peakPicking cwt snr=1.0 peakSpace=0.1 msLevel=1-2"

# Check if a file is centroid or profile using pyopenms
```

```python
import pyopenms as ms

exp = ms.MSExperiment()
ms.MzMLFile().load("sample.mzML", exp)

for i, spec in enumerate(exp):
    if i >= 3:
        break
    data_type = "centroid" if spec.getType() == 1 else "profile"
    print(f"Spectrum {i}: MS{spec.getMSLevel()}, {data_type}, {len(spec.get_peaks()[0])} peaks")
```

### 3. Encoding Precision (32-bit vs 64-bit)

64-bit encoding preserves full precision but doubles file size. 32-bit is sufficient for most metabolomics workflows.

```bash
# Force 32-bit encoding (smaller files)
msconvert sample.raw --mzML \
    --filter "peakPicking vendor msLevel=1-2" \
    --32

# Force 64-bit encoding (maximum precision)
msconvert sample.raw --mzML \
    --filter "peakPicking vendor msLevel=1-2" \
    --64

# Mixed: 64-bit m/z, 32-bit intensity (good balance)
msconvert sample.raw --mzML \
    --filter "peakPicking vendor msLevel=1-2" \
    --mz64 --inten32
```

### 4. Filtering During Conversion

msConvert filters reduce file size by selecting only needed data.

```bash
# Keep only MS1 spectra
msconvert sample.raw --mzML \
    --filter "msLevel 1"

# Keep only a retention time window (seconds)
msconvert sample.raw --mzML \
    --filter "scanTime [60,600]"

# Keep only a specific m/z range
msconvert sample.raw --mzML \
    --filter "mzWindow [100,1000]"

# Combine multiple filters (applied in order)
msconvert sample.raw --mzML \
    --filter "peakPicking vendor msLevel=1-2" \
    --filter "msLevel 1-2" \
    --filter "scanTime [60,1200]" \
    --filter "mzWindow [50,1500]"
```

### 5. Batch Conversion Scripts

**Bash (Linux/macOS with Docker):**

```bash
#!/usr/bin/env bash
# batch_convert.sh - Convert all RAW files in a directory
INPUT_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"

for raw_file in "$INPUT_DIR"/*.raw; do
    [ -f "$raw_file" ] || continue
    echo "Converting: $(basename "$raw_file")"
    docker run --rm \
        -v "$INPUT_DIR":/input \
        -v "$OUTPUT_DIR":/output \
        chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest \
        wine msconvert "/input/$(basename "$raw_file")" \
        --mzML \
        --filter "peakPicking vendor msLevel=1-2" \
        --outdir /output
done
echo "Conversion complete."
```

**Python (cross-platform):**

```python
import subprocess
import os
from pathlib import Path

def batch_convert(
    input_dir: str,
    output_dir: str,
    output_format: str = "mzML",
    centroid: bool = True,
) -> list:
    """Batch convert vendor files using msconvert."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    extensions = {".raw", ".wiff", ".d"}
    input_files = [
        f for f in input_path.iterdir()
        if f.suffix.lower() in extensions or (f.is_dir() and f.suffix.lower() in extensions)
    ]

    converted = []
    for raw_file in input_files:
        cmd = [
            "msconvert", str(raw_file),
            f"--{output_format}",
            "--outdir", str(output_path),
        ]
        if centroid:
            cmd.extend(["--filter", "peakPicking vendor msLevel=1-2"])

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            converted.append(raw_file.name)
            print(f"Converted: {raw_file.name}")
        except subprocess.CalledProcessError as e:
            print(f"Failed: {raw_file.name} - {e.stderr}")
    return converted

# Convert all RAW files in a directory
results = batch_convert("/data/raw", "/data/mzml", output_format="mzML", centroid=True)
print(f"Successfully converted {len(results)} files")
```

### 6. mzML/mzXML Interconversion with pyopenms

```python
import pyopenms as ms

def mzml_to_mzxml(input_path: str, output_path: str) -> None:
    """Convert mzML to mzXML using pyopenms."""
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_path, exp)
    ms.MzXMLFile().store(output_path, exp)
    print(f"Converted {input_path} -> {output_path} ({exp.getNrSpectra()} spectra)")

def mzxml_to_mzml(input_path: str, output_path: str) -> None:
    """Convert mzXML to mzML using pyopenms."""
    exp = ms.MSExperiment()
    ms.MzXMLFile().load(input_path, exp)
    ms.MzMLFile().store(output_path, exp)
    print(f"Converted {input_path} -> {output_path} ({exp.getNrSpectra()} spectra)")

# Convert mzML to mzXML
mzml_to_mzxml("sample.mzML", "sample.mzXML")

# Convert mzXML to mzML
mzxml_to_mzml("legacy_data.mzXML", "legacy_data.mzML")
```

### 7. Validate Converted Files

```python
import pyopenms as ms

def validate_mzml(filepath: str) -> dict:
    """Validate an mzML file and report basic statistics."""
    exp = ms.MSExperiment()
    ms.MzMLFile().load(filepath, exp)

    ms1_count = 0
    ms2_count = 0
    centroid_count = 0
    profile_count = 0

    for spec in exp:
        level = spec.getMSLevel()
        if level == 1:
            ms1_count += 1
        elif level == 2:
            ms2_count += 1
        if spec.getType() == 1:
            centroid_count += 1
        else:
            profile_count += 1

    rt_range = (
        exp.getSpectrum(0).getRT(),
        exp.getSpectrum(exp.getNrSpectra() - 1).getRT(),
    )

    return {
        "file": filepath,
        "total_spectra": exp.getNrSpectra(),
        "ms1_spectra": ms1_count,
        "ms2_spectra": ms2_count,
        "centroid_spectra": centroid_count,
        "profile_spectra": profile_count,
        "chromatograms": exp.getNrChromatograms(),
        "rt_range_seconds": rt_range,
    }

stats = validate_mzml("converted_sample.mzML")
for key, value in stats.items():
    print(f"{key}: {value}")
```

## Common Workflows

### Workflow 1: Convert a Batch of Thermo RAW Files for XCMS

1. Convert with centroiding:
   ```bash
   msconvert /data/raw/*.raw --mzML \
       --filter "peakPicking vendor msLevel=1-2" \
       --outdir /data/mzml/
   ```

2. Validate output:
   ```python
   from pathlib import Path
   for mzml in Path("/data/mzml").glob("*.mzML"):
       stats = validate_mzml(str(mzml))
       print(f"{mzml.name}: {stats['total_spectra']} spectra, "
             f"{stats['centroid_spectra']} centroid")
   ```

### Workflow 2: Convert Legacy mzXML Archive to mzML

1. Find all mzXML files:
   ```python
   from pathlib import Path
   mzxml_files = list(Path("/data/legacy").glob("**/*.mzXML"))
   print(f"Found {len(mzxml_files)} mzXML files")
   ```

2. Batch convert:
   ```python
   output_dir = Path("/data/converted")
   output_dir.mkdir(parents=True, exist_ok=True)
   for mzxml in mzxml_files:
       out_path = output_dir / mzxml.with_suffix(".mzML").name
       mzxml_to_mzml(str(mzxml), str(out_path))
   ```

## Best Practices

1. **Always centroid during conversion**: Apply vendor-native peak picking (`peakPicking vendor`) during the initial conversion step. Re-centroiding downstream is lower quality.

2. **Use vendor libraries when possible**: Vendor-native centroiding (via msConvert on Windows or Docker) produces superior results compared to algorithm-based alternatives.

3. **Prefer mzML over mzXML**: mzML is the current HUPO-PSI standard. mzXML is a legacy format that lacks support for newer features (ion mobility, chromatogram metadata).

4. **Keep original vendor files**: Always retain the original vendor files. Conversion is lossy for some metadata (tune parameters, calibration info).

5. **Check centroid/profile status after conversion**: Use the validation function to confirm spectra are in the expected mode before running analysis pipelines.

6. **Match encoding to downstream tools**: Use 64-bit m/z encoding if your pipeline performs high-precision mass calculations. Use 32-bit intensity to save space without meaningful precision loss.

7. **WIFF files need sidecar**: AB SCIEX `.wiff` files require the accompanying `.wiff.scan` file in the same directory. Conversion will fail without it.

8. **Agilent and Waters use directories**: These vendor formats are directories (`.d/`, `.raw/`), not single files. Point msConvert at the directory path.

## Resources

- **ProteoWizard / msConvert**: https://proteowizard.sourceforge.io/
- **ProteoWizard Docker image**: https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
- **mzML specification**: https://www.psidev.info/mzML
- **pyopenms documentation**: https://pyopenms.readthedocs.io/
- **HUPO-PSI standards**: https://www.psidev.info/
