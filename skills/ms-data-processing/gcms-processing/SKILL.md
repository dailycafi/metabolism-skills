---
name: gcms-processing
description: "Process GC-MS data for metabolomics and volatile compound analysis including peak detection, deconvolution, and NIST library matching. Use when: user has GC-MS data, needs retention index calculation, wants to identify volatiles, match EI spectra against NIST, or process derivatized metabolites. Triggers: GC-MS, gas chromatography, volatile analysis, NIST library, retention index, Kovats index, EI spectrum, electron ionization, deconvolution, AMDIS, TMS derivatives, volatile profiling, headspace analysis, SPME."
tool_type: python
primary_tool: matchms
---

# GC-MS Data Processing

## Overview

GC-MS metabolomics requires specialized workflows distinct from LC-MS: electron ionization (EI) fragmentation produces reproducible spectra suitable for library matching, but introduces challenges like derivatization artifacts, column bleed interference, and co-eluting peaks requiring deconvolution. This skill covers the full pipeline from raw data through metabolite identification.

## Version Compatibility

Reference examples tested with: matchms 0.18-0.26, pyopenms 3.x, xcms 4.x

The matchms API for `Scores` objects and `CosineGreedy.pair()` return types changed
across versions. If code throws `TypeError` or `KeyError`, introspect the return value
(e.g., `print(type(result), result)`) and adapt the unpacking accordingly.

## Installation

### Python Tools

```bash
uv pip install matchms pyopenms numpy pandas scipy
```

### R Tools (XCMS for GC-MS)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("xcms", "MSnbase"))
```

### Verify Installation

```python
import matchms
import pyopenms as ms
print(f"matchms: {matchms.__version__}")
print(f"pyopenms: {ms.__version__}")
```

## Core Capabilities

### 1. Loading GC-MS Data

Read raw GC-MS files (mzML, mzXML, CDF) and inspect TIC vs EIC.

```python
import pyopenms as ms
import numpy as np

exp = ms.MSExperiment()
ms.MzMLFile().load("gcms_sample.mzML", exp)

# Build total ion chromatogram (TIC)
retention_times = []
tic_intensities = []
for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        retention_times.append(spectrum.getRT())
        mz, intensity = spectrum.get_peaks()
        tic_intensities.append(np.sum(intensity))

print(f"Spectra: {len(retention_times)}")
print(f"RT range: {min(retention_times):.1f} - {max(retention_times):.1f} s")
```

Extract an extracted ion chromatogram (EIC) for a target m/z:

```python
def extract_eic(experiment, target_mz, tolerance=0.5):
    rts = []
    intensities = []
    for spectrum in experiment:
        if spectrum.getMSLevel() == 1:
            rts.append(spectrum.getRT())
            mz_array, int_array = spectrum.get_peaks()
            mask = np.abs(mz_array - target_mz) <= tolerance
            intensities.append(np.sum(int_array[mask]) if np.any(mask) else 0.0)
    return np.array(rts), np.array(intensities)

# EIC for m/z 73 (common TMS fragment)
rt_arr, eic_arr = extract_eic(exp, target_mz=73.0, tolerance=0.5)
```

### 2. Peak Detection with XCMS (R, GC-MS Optimized)

CentWave parameters tuned for GC-MS: narrower peak widths, wider ppm tolerance for unit-mass-resolution instruments.

```r
library(xcms)
library(MSnbase)

raw_files <- list.files("gcms_data", pattern = "\\.mzML$", full.names = TRUE)
raw_data <- readMSData(raw_files, mode = "onDisk")

# CentWave parameters optimized for GC-MS
# - peakwidth: GC peaks are typically 2-10 seconds (narrower than LC)
# - ppm: unit-mass GC-MS needs wider tolerance (~30-40 ppm);
#         for HR-GC-MS (e.g., GC-Orbitrap), use 5-10 ppm instead
# - noise: adjust to instrument baseline
cwp <- CentWaveParam(
    peakwidth = c(2, 10),
    ppm = 40,
    snthresh = 5,
    prefilter = c(3, 500),
    mzdiff = 0.5,
    noise = 500,
    integrate = 2
)

xdata <- findChromPeaks(raw_data, param = cwp)
cat("Peaks detected:", nrow(chromPeaks(xdata)), "\n")
```

### 3. Retention Index Calculation

Retention indices normalize retention times against an alkane ladder, enabling cross-instrument comparison.

**Kovats Retention Index** (isothermal GC, uses logarithmic interpolation):

```python
import numpy as np

def kovats_ri(rt_compound, alkane_rts, alkane_carbons):
    """Calculate Kovats retention index for isothermal GC.

    Uses the original Kovats logarithmic formula:
    RI = 100 * (n + (log(tR(x)) - log(tR(n))) / (log(tR(n+1)) - log(tR(n))))

    Args:
        rt_compound: retention time of the target compound (seconds)
        alkane_rts: list of alkane retention times (seconds), sorted ascending
        alkane_carbons: list of carbon numbers matching alkane_rts
    """
    alkane_rts = np.array(alkane_rts, dtype=np.float64)
    alkane_carbons = np.array(alkane_carbons)

    idx = np.searchsorted(alkane_rts, rt_compound) - 1
    idx = max(0, min(idx, len(alkane_rts) - 2))

    n_z = alkane_carbons[idx]
    rt_z = alkane_rts[idx]
    rt_z1 = alkane_rts[idx + 1]

    ri = 100 * (n_z + (np.log(rt_compound) - np.log(rt_z)) / (np.log(rt_z1) - np.log(rt_z)))
    return ri
```

**Linear Retention Index (LTPRI)** (temperature-programmed GC, most common in metabolomics):

```python
def linear_ri(rt_compound, alkane_rts, alkane_carbons):
    """Calculate linear (temperature-programmed) retention index.

    For temperature-programmed GC, retention times increase linearly with
    carbon number, so the simpler linear interpolation formula is used:
    RI = 100 * (n + (tR(x) - tR(n)) / (tR(n+1) - tR(n)))

    Args:
        rt_compound: retention time of the target compound (seconds)
        alkane_rts: list of alkane retention times (seconds), sorted ascending
        alkane_carbons: list of carbon numbers matching alkane_rts
    """
    alkane_rts = np.array(alkane_rts, dtype=np.float64)
    alkane_carbons = np.array(alkane_carbons)

    idx = np.searchsorted(alkane_rts, rt_compound) - 1
    idx = max(0, min(idx, len(alkane_rts) - 2))

    n_z = alkane_carbons[idx]
    rt_z = alkane_rts[idx]
    rt_z1 = alkane_rts[idx + 1]

    ri = 100 * (n_z + (rt_compound - rt_z) / (rt_z1 - rt_z))
    return ri

# Alkane ladder C8-C30
alkane_rts = [180.2, 240.5, 310.1, 385.7, 462.3, 540.8, 621.4, 703.9,
              788.1, 874.2, 962.0, 1051.5, 1142.8, 1235.6, 1330.1,
              1426.0, 1523.4, 1622.1, 1722.3, 1823.8, 1926.5, 2030.7, 2136.2]
alkane_carbons = list(range(8, 31))

ri = linear_ri(rt_compound=450.0, alkane_rts=alkane_rts, alkane_carbons=alkane_carbons)
print(f"Linear RI: {ri:.0f}")
```

### 4. Spectral Deconvolution and Library Matching with matchms

matchms supports GC-MS spectral comparison against reference libraries such as NIST and MassBank.

```python
from matchms.importing import load_from_msp
from matchms.similarity import CosineGreedy
from matchms import calculate_scores
from matchms.filtering import (
    default_filters,
    normalize_intensities,
    select_by_mz,
    select_by_relative_intensity,
)

def process_spectrum(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_mz(spectrum, mz_from=35, mz_to=500)
    spectrum = select_by_relative_intensity(spectrum, intensity_from=0.01)
    return spectrum

# Load reference library (NIST-format MSP export)
references = list(load_from_msp("nist_ei_library.msp"))
references = [process_spectrum(s) for s in references if s is not None]

# Load query spectra (deconvolved peaks exported as MSP)
queries = list(load_from_msp("sample_deconvolved.msp"))
queries = [process_spectrum(s) for s in queries if s is not None]

# Calculate cosine similarity scores
similarity = CosineGreedy(tolerance=0.5)
scores = calculate_scores(references, queries, similarity)

# Extract best matches per query
for i, query in enumerate(queries):
    best = scores.scores_by_query(query, "CosineGreedy_score", sort=True)
    # Each entry is (reference, score_object); access style depends on matchms version
    top_matches = []
    for ref, score_obj in best[:3]:
        sc = float(score_obj["score"] if isinstance(score_obj, dict) else score_obj[0])
        nm = int(score_obj["matches"] if isinstance(score_obj, dict) else score_obj[1])
        if sc > 0.7 and nm >= 6:
            top_matches.append((ref, sc, nm))
    if top_matches:
        ref, score, n_matches = top_matches[0]
        name = ref.get("compound_name", "Unknown")
        print(f"Query RT={query.get('retention_time')}: {name} "
              f"(score={score:.3f}, matched_peaks={n_matches})")
```

### 5. NIST Library Matching Workflow

For production-grade NIST matching, combine RI filtering with spectral similarity:

```python
def nist_match_with_ri(query_spectrum, query_ri, references, ri_tolerance=20):
    """Filter candidates by RI window, then rank by spectral similarity."""
    similarity = CosineGreedy(tolerance=0.5)

    candidates = []
    for ref in references:
        ref_ri = ref.get("retention_index")
        if ref_ri is not None and abs(float(ref_ri) - query_ri) <= ri_tolerance:
            score_pair = similarity.pair(ref, query_spectrum)
            if score_pair is not None:
                # matchms >= 0.18 returns a dict-like object; older versions return a tuple
                score_val = float(score_pair["score"] if isinstance(score_pair, dict) else score_pair[0])
                n_matches = int(score_pair["matches"] if isinstance(score_pair, dict) else score_pair[1])
                if score_val > 0.6 and n_matches >= 5:
                    candidates.append({
                        "name": ref.get("compound_name", "Unknown"),
                        "score": score_val,
                        "matched_peaks": n_matches,
                        "ri_diff": abs(float(ref_ri) - query_ri),
                    })

    return sorted(candidates, key=lambda x: x["score"], reverse=True)
```

### 6. Peak Detection and Integration with PyOpenMS

```python
import pyopenms as ms

exp = ms.MSExperiment()
ms.MzMLFile().load("gcms_sample.mzML", exp)

# Smooth spectra
gaussian = ms.GaussFilter()
params = gaussian.getParameters()
params.setValue("gaussian_width", 0.15)
gaussian.setParameters(params)
gaussian.filterExperiment(exp)

# Peak picking (centroiding)
picker = ms.PeakPickerHiRes()
pick_params = picker.getParameters()
pick_params.setValue("signal_to_noise", 3.0)
picker.setParameters(pick_params)

centroided = ms.MSExperiment()
picker.pickExperiment(exp, centroided)
print(f"Centroided spectra: {centroided.getNrSpectra()}")
```

## Common GC-MS Issues and Solutions

### TIC vs EIC Selection

- **TIC** (Total Ion Chromatogram): use for overall sample overview and retention time alignment. Prone to masking low-abundance compounds.
- **EIC** (Extracted Ion Chromatogram): use for targeted quantification. Select characteristic fragment ions (e.g., m/z 73, 147 for TMS derivatives; m/z 117 for amino acid TMS).

### Derivatization Artifacts

TMS (trimethylsilyl) derivatization is the most common approach for GC-MS metabolomics. Common reagents include MSTFA (N-Methyl-N-(trimethylsilyl)trifluoroacetamide) and BSTFA (N,O-Bis(trimethylsilyl)trifluoroacetamide), often with 1% TMCS as catalyst. TMS derivatization produces multiple derivatives per metabolite:

```python
# Common TMS artifacts to flag and exclude
TMS_ARTIFACTS = {
    "hexamethylcyclotrisiloxane": 207,   # D3 siloxane, column bleed
    "octamethylcyclotetrasiloxane": 281,  # D4 siloxane, column bleed
    "trimethylsilanol": 75,               # TMS reagent byproduct
}

def flag_artifacts(peak_list, artifact_mz_values, tolerance=0.5):
    """Flag peaks matching known derivatization artifacts."""
    flagged = []
    for peak in peak_list:
        base_peak_mz = peak.get("base_peak_mz", 0)
        is_artifact = any(
            abs(base_peak_mz - mz) <= tolerance
            for mz in artifact_mz_values
        )
        flagged.append({**peak, "is_artifact": is_artifact})
    return flagged
```

### Column Bleed Detection

Column bleed produces characteristic siloxane fragments that increase with temperature:

```python
COLUMN_BLEED_MZ = [207.0, 281.0, 355.0, 429.0]

def detect_column_bleed(spectrum_mz, spectrum_intensity, threshold=0.05):
    """Check if column bleed fragments dominate the spectrum."""
    total_intensity = np.sum(spectrum_intensity)
    if total_intensity == 0:
        return False
    bleed_intensity = sum(
        spectrum_intensity[np.argmin(np.abs(spectrum_mz - mz))]
        for mz in COLUMN_BLEED_MZ
        if np.min(np.abs(spectrum_mz - mz)) < 0.5
    )
    return (bleed_intensity / total_intensity) > threshold
```

## Workflow: Complete GC-MS Metabolomics Pipeline

```python
# 1. Load raw data
exp = ms.MSExperiment()
ms.MzMLFile().load("gcms_sample.mzML", exp)

# 2. Smooth and centroid
gaussian = ms.GaussFilter()
g_params = gaussian.getParameters()
g_params.setValue("gaussian_width", 0.15)
gaussian.setParameters(g_params)
gaussian.filterExperiment(exp)

picker = ms.PeakPickerHiRes()
centroided = ms.MSExperiment()
picker.pickExperiment(exp, centroided)

# 3. Export spectra for library matching
ms.MzMLFile().store("gcms_centroided.mzML", centroided)

# 4. Calculate RIs for detected peaks using alkane ladder
# 5. Match against NIST/MassBank using matchms with RI filtering
# 6. Export annotated feature table
```

## Best Practices

- **Always run alkane standards** in the same batch for reliable RI calculation
- **Set mass range to 35-500 m/z** for EI spectra (below 35 is noise)
- **Require minimum 6 matched peaks** and cosine score > 0.7 for confident library matches
- **Check for multiple TMS derivatives** of the same metabolite (mono-, di-, tri-TMS)
- **Use EIC for quantification**, TIC for qualitative screening
- **Blank subtraction** is critical: run solvent blanks and derivatization blanks
- **Monitor column bleed** at high-temperature segments of the temperature program
- **RI tolerance of 10-20 units** for intra-lab matching, 30-50 for cross-lab comparison

## Resources

- matchms documentation: https://matchms.readthedocs.io
- OpenMS documentation: https://www.openms.org
- NIST Chemistry WebBook: https://webbook.nist.gov
- Golm Metabolome Database (GMD): http://gmd.mpimp-golm.mpg.de
- MassBank spectral database: https://massbank.eu
