---
name: bio-ms-data-processing-spatial-metabolomics
description: Spatial metabolomics data processing for MALDI-MSI and DESI-MSI imaging experiments. Covers imzML parsing with pyimzML, ion image generation, tissue region segmentation, coregistration with histology, and ROI-based statistical comparison. Use when analyzing mass spectrometry imaging data to map metabolite distributions across tissue sections.
tool_type: python
primary_tool: pyimzml
---

## Version Compatibility

Reference examples tested with: pyimzML 1.5+, scikit-learn 1.3+, scikit-image 0.21+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Metabolomics

**"Map metabolite distributions across tissue sections from MSI data"** -> Parse imzML mass spectrometry imaging files, generate ion images for target m/z values, segment tissue regions by spectral similarity, and compare metabolite profiles between regions of interest.
- Python: `pyimzml` for imzML parsing, `scikit-learn` for spatial segmentation, `matplotlib` for visualization

## Installation

```bash
pip install pyimzml numpy scipy scikit-learn scikit-image matplotlib pandas
```

## Parse imzML Files

**Goal:** Read mass spectrometry imaging data stored in the standard imzML format.

**Approach:** Use pyimzML's ImzMLParser to iterate over spectra at each pixel coordinate, extracting m/z arrays and intensity arrays.

```python
from pyimzml.ImzMLParser import ImzMLParser
import numpy as np

# Parse imzML file (requires both .imzML and .ibd files in same directory)
parser = ImzMLParser("tissue_section.imzML")

# Get image dimensions
coordinates = np.array(parser.coordinates)
x_max = coordinates[:, 0].max()
y_max = coordinates[:, 1].max()
n_pixels = len(parser.coordinates)

print(f"Image dimensions: {x_max} x {y_max} pixels")
print(f"Total spectra: {n_pixels}")

# Read a single spectrum
idx = 0
mzs, intensities = parser.getspectrum(idx)
print(f"Pixel ({coordinates[idx, 0]}, {coordinates[idx, 1]})")
print(f"m/z range: {mzs.min():.2f} - {mzs.max():.2f}")
print(f"Number of peaks: {len(mzs)}")
```

## Ion Image Generation

**Goal:** Generate 2D spatial maps showing the distribution of a target metabolite across the tissue.

**Approach:** For each pixel, extract the intensity at a given m/z value (within a tolerance window), then reshape into a 2D image array.

```python
from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import matplotlib.pyplot as plt

def generate_ion_image(parser, target_mz, tolerance_da=0.05):
    """Generate a 2D ion image for a target m/z value.

    Args:
        parser: ImzMLParser instance
        target_mz: Target m/z value
        tolerance_da: Mass tolerance in Daltons (default 0.05 Da)

    Returns:
        2D numpy array with intensity values at each pixel
    """
    coords = np.array(parser.coordinates)
    x_max, y_max = coords[:, 0].max(), coords[:, 1].max()
    image = np.zeros((y_max, x_max))

    for idx, (x, y, z) in enumerate(parser.coordinates):
        mzs, intensities = parser.getspectrum(idx)
        mask = np.abs(mzs - target_mz) <= tolerance_da
        if mask.any():
            image[y - 1, x - 1] = intensities[mask].max()

    return image


parser = ImzMLParser("tissue_section.imzML")

# Generate ion images for metabolites of interest
metabolites = {
    "Glutamate": 148.0604,
    "Glucose": 203.0526,     # [M+Na]+
    "Phosphocholine": 184.0733,
    "Sphingomyelin": 731.5461,
}

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, (name, mz) in zip(axes.flat, metabolites.items()):
    ion_img = generate_ion_image(parser, mz, tolerance_da=0.03)
    im = ax.imshow(ion_img, cmap="viridis", interpolation="nearest")
    ax.set_title(f"{name} (m/z {mz:.4f})")
    ax.axis("off")
    plt.colorbar(im, ax=ax, fraction=0.046)

plt.tight_layout()
plt.savefig("ion_images.png", dpi=300, bbox_inches="tight")
plt.close()
```

## Build Datacube from imzML

**Goal:** Convert sparse imzML spectra into a regular 3D datacube for downstream analysis.

**Approach:** Bin all spectra onto a common m/z axis, then populate a (rows x cols x bins) array.

```python
from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
from scipy.interpolate import interp1d

def build_datacube(parser, mz_min=100, mz_max=1000, bin_width=0.1):
    """Build a 3D datacube from imzML with uniform m/z binning.

    Args:
        parser: ImzMLParser instance
        mz_min: Lower m/z bound
        mz_max: Upper m/z bound
        bin_width: Width of each m/z bin in Da

    Returns:
        datacube: 3D array (rows, cols, n_bins)
        mz_axis: 1D array of bin centers
    """
    mz_axis = np.arange(mz_min, mz_max, bin_width)
    n_bins = len(mz_axis)

    coords = np.array(parser.coordinates)
    x_max, y_max = coords[:, 0].max(), coords[:, 1].max()
    datacube = np.zeros((y_max, x_max, n_bins), dtype=np.float32)

    for idx, (x, y, z) in enumerate(parser.coordinates):
        mzs, intensities = parser.getspectrum(idx)
        binned = np.histogram(mzs, bins=n_bins, range=(mz_min, mz_max),
                              weights=intensities)[0]
        datacube[y - 1, x - 1, :] = binned

    return datacube, mz_axis


parser = ImzMLParser("tissue_section.imzML")
datacube, mz_axis = build_datacube(parser, mz_min=100, mz_max=1000, bin_width=0.1)
print(f"Datacube shape: {datacube.shape}")
```

## Spatial Segmentation

**Goal:** Segment the tissue into regions with distinct metabolic profiles.

**Approach:** Flatten the datacube to a 2D matrix (pixels x m/z bins), apply dimensionality reduction with PCA, then cluster with k-means or hierarchical clustering.

```python
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def spatial_segmentation(datacube, n_clusters=5, method="kmeans"):
    """Segment tissue by spectral similarity.

    Args:
        datacube: 3D array (rows, cols, n_bins)
        n_clusters: Number of segments
        method: 'kmeans' or 'hierarchical'

    Returns:
        2D label map (rows, cols)
    """
    rows, cols, n_bins = datacube.shape
    flat = datacube.reshape(-1, n_bins)

    # Mask background (zero-intensity pixels)
    mask = flat.sum(axis=1) > 0
    data = flat[mask]

    # Normalize and reduce dimensions
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)

    pca = PCA(n_components=20)
    data_pca = pca.fit_transform(data_scaled)
    print(f"PCA explained variance: {pca.explained_variance_ratio_.sum():.2%}")

    # Cluster
    if method == "kmeans":
        labels = KMeans(n_clusters=n_clusters, n_init=10, random_state=42).fit_predict(data_pca)
    elif method == "hierarchical":
        labels = AgglomerativeClustering(n_clusters=n_clusters).fit_predict(data_pca)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Map back to 2D
    label_map = np.full(rows * cols, -1)
    label_map[mask] = labels
    label_map = label_map.reshape(rows, cols)

    return label_map


label_map = spatial_segmentation(datacube, n_clusters=5, method="kmeans")

# Visualize segmentation
cmap = ListedColormap(["#f0f0f0", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"])
plt.figure(figsize=(8, 6))
plt.imshow(label_map, cmap=cmap, interpolation="nearest")
plt.colorbar(label="Cluster")
plt.title("Spatial Segmentation (k=5)")
plt.axis("off")
plt.savefig("segmentation.png", dpi=300, bbox_inches="tight")
plt.close()
```

## Coregistration with H&E Histology

**Goal:** Overlay MSI ion images onto H&E stained histology for anatomical context.

**Approach:** Load the histology image, apply affine transformation to align MSI coordinates, then blend the ion image as a semi-transparent overlay.

```python
import numpy as np
from skimage import io, transform
import matplotlib.pyplot as plt

def coregister_msi_histology(ion_image, histology_path,
                              scale_factor=1.0, offset_x=0, offset_y=0):
    """Overlay ion image on H&E histology.

    Args:
        ion_image: 2D array from generate_ion_image()
        histology_path: Path to H&E image file
        scale_factor: Scale MSI to histology resolution
        offset_x, offset_y: Translation offsets in pixels

    Returns:
        Blended overlay image
    """
    histology = io.imread(histology_path)

    # Resize ion image to match histology dimensions
    ion_resized = transform.resize(
        ion_image,
        (int(ion_image.shape[0] * scale_factor),
         int(ion_image.shape[1] * scale_factor)),
        anti_aliasing=True
    )

    # Create padded array matching histology size
    overlay = np.zeros(histology.shape[:2])
    h = min(ion_resized.shape[0], histology.shape[0] - offset_y)
    w = min(ion_resized.shape[1], histology.shape[1] - offset_x)
    overlay[offset_y:offset_y + h, offset_x:offset_x + w] = ion_resized[:h, :w]

    # Normalize overlay to [0, 1]
    if overlay.max() > 0:
        overlay = overlay / overlay.max()

    return overlay, histology


ion_img = generate_ion_image(parser, target_mz=184.0733)
overlay, histology = coregister_msi_histology(
    ion_img, "HE_stain.png", scale_factor=8.0, offset_x=50, offset_y=30
)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
axes[0].imshow(histology)
axes[0].set_title("H&E Histology")
axes[1].imshow(overlay, cmap="hot")
axes[1].set_title("Phosphocholine Ion Image")
axes[2].imshow(histology)
axes[2].imshow(overlay, cmap="hot", alpha=0.5)
axes[2].set_title("Overlay")
for ax in axes:
    ax.axis("off")
plt.savefig("coregistration.png", dpi=300, bbox_inches="tight")
plt.close()
```

## ROI Extraction and Statistical Comparison

**Goal:** Compare metabolite profiles between user-defined or segmentation-derived tissue regions.

**Approach:** Extract mean spectra from each ROI, then run statistical tests (t-test, fold change) to identify discriminating m/z features.

```python
import numpy as np
from scipy import stats
import pandas as pd

def compare_rois(datacube, mz_axis, label_map, roi_a, roi_b):
    """Compare metabolite profiles between two ROIs.

    Args:
        datacube: 3D array (rows, cols, n_bins)
        mz_axis: 1D array of m/z bin centers
        label_map: 2D cluster labels from segmentation
        roi_a, roi_b: Cluster labels to compare

    Returns:
        DataFrame with m/z, fold_change, p_value, fdr columns
    """
    mask_a = label_map == roi_a
    mask_b = label_map == roi_b

    spectra_a = datacube[mask_a]
    spectra_b = datacube[mask_b]

    results = []
    for i, mz in enumerate(mz_axis):
        vals_a = spectra_a[:, i]
        vals_b = spectra_b[:, i]

        mean_a = vals_a.mean()
        mean_b = vals_b.mean()
        fc = (mean_a + 1e-10) / (mean_b + 1e-10)

        _, pval = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        results.append({"mz": mz, "mean_roi_a": mean_a, "mean_roi_b": mean_b,
                        "fold_change": fc, "p_value": pval})

    df = pd.DataFrame(results)

    # FDR correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(df["p_value"], method="fdr_bh")
    df["fdr"] = fdr

    # Sort by significance
    df = df.sort_values("p_value")
    return df


comparison = compare_rois(datacube, mz_axis, label_map, roi_a=1, roi_b=3)
significant = comparison[comparison["fdr"] < 0.05]
print(f"Significant features (FDR < 0.05): {len(significant)}")
print(significant.head(20))

significant.to_csv("roi_comparison.csv", index=False)
```

## DESI-MSI Overview

DESI (Desorption Electrospray Ionization) MSI operates under ambient conditions and produces imzML files compatible with the same processing pipeline above. Key differences:
- Spatial resolution typically 50-200 um (vs 10-50 um for MALDI)
- Softer ionization preserves intact lipids and drug metabolites
- No matrix application required
- Data files use the same imzML/ibd format; parse identically with pyimzML

## Best Practices

- Use TIC (total ion current) normalization before comparing spectra across pixels
- Apply mass recalibration if m/z drift exceeds 10 ppm across the acquisition
- For MALDI: check for matrix cluster interference in the low m/z range (< 300 Da)
- Set mass tolerance based on instrument resolution (TOF: 20 ppm, Orbitrap: 5 ppm, FTICR: 2 ppm)
- Remove background pixels (no tissue) before segmentation to avoid spurious clusters
- Export ion images as TIFF for downstream analysis in ImageJ/FIJI
- For large datasets (>10k pixels), use incremental PCA (`sklearn.decomposition.IncrementalPCA`)

## Related Skills

- ms-data-processing/pyopenms - General MS data processing
- ms-data-processing/matchms - Spectral matching and annotation
- metabolomics-analysis/metabolite-annotation - Identify spatial features
