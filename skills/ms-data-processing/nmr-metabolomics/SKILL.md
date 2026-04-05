---
name: nmr-metabolomics
description: NMR metabolomics processing pipeline for 1D and 2D NMR data. Covers phase correction, baseline correction, chemical shift referencing, spectral binning, peak picking, metabolite identification, and quantification. Use for processing 1H-NMR metabolomics data, identifying metabolites from HMDB/BMRB databases, and performing absolute or relative quantification.
tool_type: python
primary_tool: nmrglue
---

# NMR Metabolomics

## Overview

NMR spectroscopy provides non-destructive, highly reproducible metabolite profiling with minimal sample preparation. 1D 1H-NMR is the workhorse for metabolomics, offering broad metabolite coverage in biofluids (serum, urine, CSF). This skill covers the complete NMR metabolomics pipeline from raw FID processing through metabolite quantification.

## Installation

### Python Tools

```bash
uv pip install nmrglue numpy scipy pandas scikit-learn matplotlib
```

### R Tools

```r
install.packages("speaq")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("PepsNMR")
```

### Verify Installation

```python
import nmrglue as ng
print(f"nmrglue version: {ng.__version__}")
```

## Core Capabilities

### 1. Loading and Processing Raw FID Data

Read Bruker, Varian/Agilent, or JCAMP-DX format NMR data.

```python
import nmrglue as ng
import numpy as np

# Load Bruker 1D 1H-NMR data
dic, fid = ng.bruker.read("bruker_experiment_dir/")

# Apply spectral processing
# Zero-fill to double the points for better digital resolution
fid = ng.proc_base.zf_size(fid, 65536)

# Apply exponential line broadening (0.3 Hz typical for metabolomics)
fid = ng.proc_base.em(fid, lb=0.3)

# Fourier transform
spectrum = ng.proc_base.fft(fid)

# Phase correction (automatic)
spectrum = ng.proc_autophase.autops(spectrum, fn="acme")

print(f"Spectrum points: {len(spectrum)}")
print(f"Data type: {spectrum.dtype}")
```

### 2. Phase Correction

Accurate phasing is critical. Two approaches: automatic and manual.

```python
# Automatic phase correction (ACME algorithm)
spectrum_phased = ng.proc_autophase.autops(spectrum, fn="acme")

# Alternative: peak_minima method (better for crowded spectra)
spectrum_phased = ng.proc_autophase.autops(spectrum, fn="peak_minima")

# Manual phase correction if automatic fails
# p0 = zero-order phase (degrees), p1 = first-order phase (degrees)
spectrum_phased = ng.proc_base.ps(spectrum, p0=30.0, p1=-10.0)
```

### 3. Baseline Correction

Remove broad baseline distortions that affect quantification accuracy.

```python
from scipy import sparse
from scipy.sparse.linalg import spsolve

def baseline_als(spectrum_real, lam=1e7, p=0.001, n_iter=10):
    """Asymmetric least squares baseline correction.

    Args:
        spectrum_real: real part of the spectrum (1D array)
        lam: smoothness parameter (1e5 to 1e9; larger = smoother)
        p: asymmetry parameter (0.001 to 0.01; smaller = more asymmetric)
        n_iter: number of iterations
    """
    length = len(spectrum_real)
    diag = sparse.diags([1, -2, 1], [0, -1, -2], shape=(length, length - 2))
    d_matrix = lam * diag.dot(diag.transpose())
    weights = np.ones(length)
    for _ in range(n_iter):
        w_diag = sparse.spdiags(weights, 0, length, length)
        z_matrix = w_diag + d_matrix
        baseline = spsolve(z_matrix, weights * spectrum_real)
        weights = p * (spectrum_real > baseline) + (1 - p) * (spectrum_real <= baseline)
    return spectrum_real - baseline

spectrum_real = np.real(spectrum_phased)
spectrum_corrected = baseline_als(spectrum_real, lam=1e7, p=0.005)
```

### 4. Chemical Shift Referencing

Reference to internal standard (TSP at 0.00 ppm for aqueous, DSS for D2O).

```python
def reference_to_tsp(ppm_scale, spectrum_real, tsp_region=(-0.05, 0.05)):
    """Reference spectrum to TSP signal at 0.00 ppm.

    Args:
        ppm_scale: chemical shift axis (array)
        spectrum_real: real spectrum intensities (array)
        tsp_region: ppm range to search for TSP peak
    """
    mask = (ppm_scale >= tsp_region[0]) & (ppm_scale <= tsp_region[1])
    tsp_idx = np.argmax(spectrum_real[mask])
    tsp_ppm = ppm_scale[mask][tsp_idx]
    shift = 0.0 - tsp_ppm
    corrected_ppm = ppm_scale + shift
    return corrected_ppm

# Generate ppm scale from Bruker parameters
udic = ng.bruker.guess_udic(dic, fid)
uc = ng.fileiobase.uc_from_udic(udic)
ppm_scale = uc.ppm_scale()

ppm_referenced = reference_to_tsp(ppm_scale, spectrum_corrected)
```

### 5. Spectral Binning

#### Uniform Binning

```python
def uniform_binning(ppm_scale, spectrum_real, bin_width=0.04,
                    ppm_range=(0.5, 10.0), exclude_water=(4.7, 5.0)):
    """Bin spectrum into uniform-width segments.

    Args:
        ppm_scale: chemical shift axis
        spectrum_real: real spectrum intensities
        bin_width: bin width in ppm (0.04 ppm is standard)
        ppm_range: spectral region to bin
        exclude_water: water resonance region to exclude
    """
    bins = np.arange(ppm_range[0], ppm_range[1], bin_width)
    bin_intensities = []
    bin_centers = []

    for bin_start in bins:
        bin_end = bin_start + bin_width
        # Skip water region
        if (bin_start >= exclude_water[0] and bin_start <= exclude_water[1]):
            continue
        mask = (ppm_scale >= bin_start) & (ppm_scale < bin_end)
        if np.any(mask):
            bin_intensities.append(np.sum(spectrum_real[mask]))
            bin_centers.append(bin_start + bin_width / 2)

    return np.array(bin_centers), np.array(bin_intensities)

centers, intensities = uniform_binning(ppm_referenced, spectrum_corrected)
print(f"Number of bins: {len(centers)}")
```

#### Adaptive Intelligent Binning (R with speaq)

```r
library(speaq)

# Load spectra matrix (rows = samples, columns = ppm points)
spectra_matrix <- as.matrix(read.csv("spectra_matrix.csv", row.names = 1))
ppm <- as.numeric(colnames(spectra_matrix))

# Adaptive intelligent binning using speaq
# Groups spectral data points that belong to the same signal
peaks <- getWaveletPeaks(
    Y.spec = spectra_matrix,
    X.ppm = ppm,
    baselineThresh = 1000,
    SNR.Th = 3,
    nCPU = 4
)

# Group peaks across samples
grouped <- PeakGrouper(
    Y.peaks = peaks,
    grouping.window.width = 100,
    verbose = FALSE
)

# Build feature matrix from grouped peaks
feature_matrix <- BuildFeatureMatrix(grouped)
cat("Features:", ncol(feature_matrix), "\n")
cat("Samples:", nrow(feature_matrix), "\n")
```

### 6. Peak Picking and Metabolite Identification

```python
from scipy.signal import find_peaks

def pick_nmr_peaks(ppm_scale, spectrum_real, height_threshold=1e5,
                   min_distance_ppm=0.01):
    """Pick peaks from processed 1H-NMR spectrum.

    Args:
        ppm_scale: chemical shift axis
        spectrum_real: processed real spectrum
        height_threshold: minimum peak height
        min_distance_ppm: minimum distance between peaks in ppm
    """
    # Calculate minimum distance in data points
    ppm_per_point = abs(ppm_scale[1] - ppm_scale[0])
    min_distance_pts = max(1, int(min_distance_ppm / ppm_per_point))

    indices, properties = find_peaks(
        spectrum_real,
        height=height_threshold,
        distance=min_distance_pts,
        prominence=height_threshold * 0.1,
    )

    peaks = [
        {"ppm": float(ppm_scale[idx]), "intensity": float(spectrum_real[idx])}
        for idx in indices
    ]
    return sorted(peaks, key=lambda p: p["ppm"])

peaks = pick_nmr_peaks(ppm_referenced, spectrum_corrected, height_threshold=5e4)
print(f"Peaks found: {len(peaks)}")
```

#### Metabolite Identification Against HMDB/BMRB

```python
# Reference chemical shifts from HMDB for common metabolites (1H, D2O, pH 7.4)
HMDB_REFERENCES = {
    "lactate": {"shifts": [1.33, 4.11], "multiplicity": ["d", "q"]},
    "alanine": {"shifts": [1.48, 3.78], "multiplicity": ["d", "q"]},
    "glucose_alpha": {"shifts": [5.23, 3.53, 3.71, 3.42, 3.84, 3.76],
                      "multiplicity": ["d", "dd", "t", "dd", "ddd", "dd"]},
    "creatinine": {"shifts": [3.04, 4.05], "multiplicity": ["s", "s"]},
    "citrate": {"shifts": [2.52, 2.68], "multiplicity": ["d", "d"]},
    "valine": {"shifts": [0.99, 1.04, 2.27, 3.62], "multiplicity": ["d", "d", "m", "d"]},
}

def match_metabolites(observed_peaks, references, tolerance=0.03):
    """Match observed peaks to reference metabolite chemical shifts.

    Args:
        observed_peaks: list of dicts with 'ppm' keys
        references: dict of metabolite reference shifts
        tolerance: matching tolerance in ppm
    """
    observed_ppm = [p["ppm"] for p in observed_peaks]
    matches = []

    for metabolite, ref_data in references.items():
        ref_shifts = ref_data["shifts"]
        matched_shifts = []
        for ref_ppm in ref_shifts:
            diffs = [abs(obs - ref_ppm) for obs in observed_ppm]
            if min(diffs) <= tolerance:
                matched_shifts.append(ref_ppm)

        fraction = len(matched_shifts) / len(ref_shifts)
        if fraction >= 0.5:
            matches.append({
                "metabolite": metabolite,
                "matched_fraction": fraction,
                "matched_shifts": matched_shifts,
                "total_shifts": len(ref_shifts),
            })

    return sorted(matches, key=lambda m: m["matched_fraction"], reverse=True)

identifications = match_metabolites(peaks, HMDB_REFERENCES)
for hit in identifications:
    print(f"{hit['metabolite']}: {hit['matched_fraction']:.0%} "
          f"({len(hit['matched_shifts'])}/{hit['total_shifts']} peaks)")
```

### 7. 2D NMR for Metabolite Confirmation

2D experiments (HSQC, TOCSY) resolve overlapping 1D signals:

```python
# Load 2D HSQC Bruker data
dic_2d, data_2d = ng.bruker.read("hsqc_experiment_dir/")

# Process indirect dimension (13C)
data_2d = ng.proc_base.zf_size(data_2d, (512, 4096))
data_2d = ng.proc_base.fft2(data_2d)
data_2d = ng.proc_autophase.autops(data_2d, fn="acme", dim=1)

# 2D HSQC cross-peaks confirm metabolite assignments
# e.g., lactate: 1H 1.33 ppm / 13C 20.8 ppm (methyl)
#        lactate: 1H 4.11 ppm / 13C 69.3 ppm (methine)
HSQC_CONFIRMATIONS = {
    "lactate": [(1.33, 20.8), (4.11, 69.3)],
    "alanine": [(1.48, 17.4), (3.78, 51.5)],
    "glucose_alpha": [(5.23, 93.0), (3.53, 72.3)],
}
```

### 8. Quantification

#### Absolute Quantification with Internal Standard

```python
def absolute_quantification(peak_area, is_area, is_concentration_mm,
                            is_protons=9, analyte_protons=1):
    """Calculate absolute metabolite concentration using internal standard.

    Args:
        peak_area: integrated area of the metabolite signal
        is_area: integrated area of the internal standard (e.g., TSP)
        is_concentration_mm: concentration of IS in mM
        is_protons: number of equivalent protons in IS signal (TSP = 9)
        analyte_protons: number of equivalent protons in analyte signal
    """
    concentration = (peak_area / is_area) * (is_protons / analyte_protons) * is_concentration_mm
    return concentration

# Example: quantify lactate doublet at 1.33 ppm using TSP (0.5 mM)
lactate_mm = absolute_quantification(
    peak_area=2.5e6,
    is_area=1.8e6,
    is_concentration_mm=0.5,
    is_protons=9,
    analyte_protons=3,  # 3 equivalent methyl protons
)
print(f"Lactate concentration: {lactate_mm:.2f} mM")
```

#### Probabilistic Quotient Normalization (PQN)

```python
def pqn_normalization(spectra_matrix):
    """Probabilistic Quotient Normalization for relative quantification.

    Args:
        spectra_matrix: numpy array (samples x spectral_points)

    Returns:
        Normalized spectra matrix
    """
    # Step 1: integral normalization to total area
    row_sums = np.sum(spectra_matrix, axis=1, keepdims=True)
    integral_normalized = spectra_matrix / row_sums

    # Step 2: calculate reference spectrum (median across samples)
    reference = np.median(integral_normalized, axis=0)

    # Step 3: calculate quotients for each sample
    # Avoid division by zero
    ref_nonzero = reference > 0
    quotients = np.zeros(spectra_matrix.shape[0])
    for i in range(spectra_matrix.shape[0]):
        sample_quotients = integral_normalized[i, ref_nonzero] / reference[ref_nonzero]
        quotients[i] = np.median(sample_quotients)

    # Step 4: divide each spectrum by its quotient
    normalized = integral_normalized / quotients[:, np.newaxis]
    return normalized

# Apply PQN to a batch of spectra
# spectra_matrix shape: (n_samples, n_points)
spectra_pqn = pqn_normalization(spectra_matrix)
```

## R Pipeline with PepsNMR

```r
library(PepsNMR)

# Load Bruker FIDs
fid_list <- ReadFids(
    path = "bruker_data_dir",
    l = 1,
    subdirs = TRUE
)
fid_data <- fid_list[["Fid_data"]]
fid_info <- fid_list[["Fid_info"]]

# Processing pipeline
fid_ss <- SolventSuppression(fid_data, lambda.ss = 1e6)
fid_apod <- Apodization(fid_ss, DT = fid_info[1, "DT"], lambda.apod = 0.3)
spectrum <- FourierTransform(fid_apod, Fid_info = fid_info)
spectrum <- ZeroOrderPhaseCorrection(spectrum, Fid_info = fid_info)
spectrum <- InternalReferencing(spectrum, Fid_info = fid_info, ppm.ref = 0.0)
spectrum <- BaselineCorrection(spectrum, lambda.bc = 1e7, p.bc = 0.001)
spectrum <- NegativeValuesZeroing(spectrum)

# Bucketing (uniform binning)
bucketed <- Bucketing(spectrum, width = TRUE, mb = 0.04)
cat("Bucketed features:", ncol(bucketed), "\n")
```

## Best Practices

- **Always reference to TSP** (0.00 ppm) for aqueous samples or DSS for organic/D2O samples
- **Exclude the water region** (4.7-5.0 ppm) from all analyses to avoid artifacts
- **Use 0.3 Hz line broadening** as default; increase to 0.5-1.0 Hz for noisy spectra
- **PQN normalization** is preferred over total area normalization for biofluids
- **Confirm assignments with 2D NMR** (HSQC, TOCSY) when possible, especially for overlapping regions
- **Match at least 50% of reference peaks** before calling an identification confident
- **Bin width of 0.04 ppm** is standard; use adaptive binning for better resolution
- **pH affects chemical shifts**: buffer samples or apply pH correction models for urine
- **Relaxation delay of 4-5 seconds** is needed for quantitative 1H-NMR (T1 recovery)
- **Temperature control at 298 K** ensures chemical shift reproducibility

## Resources

- HMDB NMR spectra: https://hmdb.ca
- BMRB metabolomics: https://bmrb.io
- nmrglue documentation: https://nmrglue.readthedocs.io
- speaq R package: https://cran.r-project.org/package=speaq
- PepsNMR: https://bioconductor.org/packages/PepsNMR
