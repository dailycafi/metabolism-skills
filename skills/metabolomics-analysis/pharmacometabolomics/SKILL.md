---
name: bio-metabolomics-analysis-pharmacometabolomics
description: Pharmacometabolomics workflow for drug response profiling, CYP450 substrate prediction, ADME property estimation, metabolic soft spot prediction, and phase I/II metabolite identification from LC-MS data. Integrates DrugBank API, RDKit molecular descriptors, and MetaboAnalystR for differential analysis. Use when studying drug metabolism, predicting drug-metabolite interactions, or profiling pharmacokinetic responses.
tool_type: python
primary_tool: rdkit
---

## Version Compatibility

Reference examples tested with: RDKit 2023.09+, MetaboAnalystR 4.0+, pandas 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion("<pkg>")` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pharmacometabolomics

**"Profile drug response via metabolomics and predict drug metabolism"** -> Compare pre-dose vs post-dose metabolite profiles, predict CYP450 enzyme interactions and ADME properties from molecular structure, identify phase I/II metabolites in LC-MS data, and integrate PK/PD parameters with metabolomics measurements.
- Python: `rdkit` for molecular descriptors, `pandas`/`scipy` for statistics
- R: `MetaboAnalystR` for differential metabolite analysis

## Installation

```bash
pip install rdkit pandas numpy scipy scikit-learn matplotlib requests

# R package (for MetaboAnalystR section) — NOT on CRAN, install from GitHub
Rscript -e 'devtools::install_github("xia-lab/MetaboAnalystR", build_vignettes = FALSE)'
```

## Drug Metabolite Profiling: Pre-dose vs Post-dose

**Goal:** Identify metabolites that change significantly after drug administration.

**Approach:** Perform paired statistical testing between pre-dose and post-dose metabolite levels, apply multiple testing correction, and generate volcano plots.

```python
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

def drug_response_analysis(pre_dose_file, post_dose_file, alpha=0.05):
    """Compare metabolite profiles before and after drug treatment."""
    pre = pd.read_csv(pre_dose_file, index_col=0)
    post = pd.read_csv(post_dose_file, index_col=0)

    common = pre.index.intersection(post.index)
    pre, post = pre.loc[common], post.loc[common]

    results = []
    for metabolite in pre.columns:
        pre_vals, post_vals = pre[metabolite].values, post[metabolite].values
        t_stat, p_val = stats.ttest_rel(post_vals, pre_vals)
        log2fc = np.log2((post_vals.mean() + 1e-10) / (pre_vals.mean() + 1e-10))
        results.append({"metabolite": metabolite, "log2_fold_change": log2fc,
                        "t_statistic": t_stat, "p_value": p_val})

    df = pd.DataFrame(results)
    _, fdr, _, _ = multipletests(df["p_value"], method="fdr_bh")
    df["fdr"] = fdr
    df["significant"] = (df["fdr"] < alpha) & (np.abs(df["log2_fold_change"]) > 0.5)
    return df.sort_values("p_value")


results = drug_response_analysis("pre_dose.csv", "post_dose.csv")
print(f"Significant changes: {results['significant'].sum()}")

# Volcano plot
fig, ax = plt.subplots(figsize=(8, 6))
colors = ["#b2182b" if s else "#cccccc" for s in results["significant"]]
ax.scatter(results["log2_fold_change"], -np.log10(results["p_value"]), c=colors, s=20)
ax.set_xlabel("Log2 Fold Change (Post/Pre)")
ax.set_ylabel("-Log10(p-value)")
plt.savefig("volcano_drug_response.png", dpi=300, bbox_inches="tight")
plt.close()
```

## CYP450 Substrate Prediction

**Goal:** Predict which CYP450 enzymes are likely to metabolize a given drug.

**Approach:** Compute RDKit molecular descriptors relevant to CYP450 binding (lipophilicity, molecular weight, H-bond donors/acceptors, aromatic rings), then apply rule-based filtering for major CYP isoforms.

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def compute_cyp450_descriptors(smiles):
    """Compute molecular descriptors relevant to CYP450 metabolism."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return {
        "MW": Descriptors.MolWt(mol), "LogP": Descriptors.MolLogP(mol),
        "HBA": Descriptors.NumHAcceptors(mol), "TPSA": Descriptors.TPSA(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def predict_cyp_isoforms(descriptors):
    """Rule-based CYP450 isoform prediction based on substrate profiles."""
    predictions = []
    mw, logp, hba = descriptors["MW"], descriptors["LogP"], descriptors["HBA"]
    aromatic = descriptors["AromaticRings"]

    if mw > 300 and logp > 1.0:  # CYP3A4: larger lipophilic molecules
        predictions.append(("CYP3A4", round(min(1.0, (mw - 300) / 500 * 0.5 + logp / 6 * 0.5), 2)))
    if 200 < mw < 500 and 1.0 < logp < 4.0:  # CYP2D6: basic amines
        predictions.append(("CYP2D6", round(0.6 if aromatic >= 1 else 0.3, 2)))
    if 200 < mw < 500 and hba >= 2:  # CYP2C9: acidic compounds
        predictions.append(("CYP2C9", round(0.5 if logp > 2.0 else 0.3, 2)))
    if 200 < mw < 450 and 1.5 < logp < 5.0:  # CYP2C19
        predictions.append(("CYP2C19", 0.4))
    if aromatic >= 2 and mw < 400:  # CYP1A2: planar aromatics
        predictions.append(("CYP1A2", round(0.5 + aromatic * 0.1, 2)))

    return sorted(predictions, key=lambda x: x[1], reverse=True)


# Example: Ibuprofen
desc = compute_cyp450_descriptors("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
for isoform, score in predict_cyp_isoforms(desc):
    print(f"  {isoform}: {score}")
```

## ADME Property Prediction

**Goal:** Estimate absorption, distribution, metabolism, and excretion properties from molecular structure.

**Approach:** Compute Lipinski and Veber rule descriptors, predict Caco-2 permeability class, and estimate aqueous solubility from LogP.

```python
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

def predict_adme(smiles):
    """Predict ADME properties from molecular structure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotbonds = Descriptors.NumRotatableBonds(mol)

    lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_pass = tpsa <= 140 and rotbonds <= 10

    # ESOL solubility model (Delaney 2004)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    arom_prop = aromatic_atoms / max(mol.GetNumHeavyAtoms(), 1)
    log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotbonds - 0.74 * arom_prop

    sol_class = ("Highly soluble" if log_s > -2 else "Soluble" if log_s > -4
                 else "Moderately soluble" if log_s > -6 else "Poorly soluble")

    return {
        "MW": round(mw, 1), "LogP": round(logp, 2),
        "HBD": hbd, "HBA": hba, "TPSA": round(tpsa, 1),
        "Lipinski_pass": lipinski_violations <= 1, "Veber_pass": veber_pass,
        "LogS_ESOL": round(log_s, 2), "Solubility": sol_class,
        "GI_absorption": "High" if tpsa < 140 else "Low",
        "BBB_permeant": "Yes" if tpsa < 80 and logp > 0 else "No",
    }


# Example: Metformin
for prop, val in predict_adme("CN(C)C(=N)NC(=N)N").items():
    print(f"  {prop}: {val}")
```

## Metabolic Soft Spot Prediction

**Goal:** Identify likely sites of metabolic oxidation on a drug molecule.

**Approach:** Compute atomic partial charges and identify electron-rich centers (aromatic carbons, allylic positions, benzylic carbons) that are susceptible to CYP-mediated oxidation.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def predict_soft_spots(smiles):
    """Predict metabolic soft spots (likely sites of CYP oxidation)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    AllChem.ComputeGasteigerCharges(mol)
    soft_spots = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        idx, symbol = atom.GetIdx(), atom.GetSymbol()
        reasons = []

        if atom.GetIsAromatic() and symbol == "C":
            reasons.append("aromatic_hydroxylation")
        if symbol == "C" and not atom.GetIsAromatic():
            if any(n.GetIsAromatic() for n in atom.GetNeighbors()):
                reasons.append("benzylic_oxidation")
        if symbol == "N" and atom.GetDegree() == 3 and not atom.GetIsAromatic():
            reasons.append("N_dealkylation")
        if symbol == "O" and atom.GetDegree() == 2:
            if "C" in [n.GetSymbol() for n in atom.GetNeighbors()]:
                reasons.append("O_dealkylation")
        if symbol == "S" and atom.GetDegree() <= 2:
            reasons.append("S_oxidation")

        if reasons:
            soft_spots.append((idx, symbol, reasons))

    return soft_spots


# Example: Verapamil
for idx, symbol, reasons in predict_soft_spots("COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC"):
    print(f"  Atom {idx} ({symbol}): {', '.join(reasons)}")
```

## Phase I/II Metabolite Identification from LC-MS

**Goal:** Identify expected phase I and phase II metabolites in LC-MS data by computing theoretical m/z shifts.

**Approach:** Define mass shifts for common biotransformations, then search extracted ion chromatograms for matching peaks.

```python
import pandas as pd
import numpy as np

# Common phase I/II biotransformation mass shifts (Da)
BIOTRANSFORMATIONS = {
    "Hydroxylation (+O)": 15.9949, "Dihydroxylation (+2O)": 31.9898,
    "Dehydrogenation (-2H)": -2.0157, "N-dealkylation (-CH2)": -14.0157,
    "Reduction (+2H)": 2.0157,
    "Glucuronidation (+C6H8O6)": 176.0321, "Sulfation (+SO3)": 79.9568,
    "Glutathione (+GSH)": 305.0682, "Acetylation (+C2H2O)": 42.0106,
    "Methylation (+CH2)": 14.0157, "Glycine conjugation": 57.0215,
}


def search_metabolites(parent_mz, ms_data_file, ppm_tolerance=10):
    """Search LC-MS data for predicted drug metabolites."""
    ms_data = pd.read_csv(ms_data_file)
    candidates = []
    for name, shift in BIOTRANSFORMATIONS.items():
        expected = parent_mz + shift
        tol = expected * ppm_tolerance / 1e6
        matches = ms_data[(ms_data["mz"] >= expected - tol) &
                          (ms_data["mz"] <= expected + tol)]
        if len(matches) > 0:
            best = matches.loc[matches["intensity"].idxmax()]
            candidates.append({
                "transformation": name, "expected_mz": round(expected, 4),
                "observed_mz": round(best["mz"], 4),
                "ppm_error": round((best["mz"] - expected) / expected * 1e6, 1),
                "retention_time": round(best["rt"], 2), "intensity": best["intensity"],
            })
    return pd.DataFrame(candidates).sort_values("intensity", ascending=False)


# Example: Search for metabolites of verapamil ([M+H]+ = 455.2904)
metabolites = search_metabolites(455.2904, "lcms_postdose.csv", ppm_tolerance=5)
print(metabolites.to_string(index=False))
```

## DrugBank API Integration

**Goal:** Query DrugBank for known drug-metabolite interactions and enzyme associations.

**Approach:** Use the DrugBank API to retrieve metabolizing enzymes, transporters, and known metabolites for a drug.

```python
import os
import requests

# DrugBank API requires a paid subscription (https://go.drugbank.com/api)
DRUGBANK_API_KEY = os.environ.get("DRUGBANK_API_KEY")
if not DRUGBANK_API_KEY:
    raise EnvironmentError("Set DRUGBANK_API_KEY environment variable (requires paid subscription)")

def query_drugbank(drugbank_id):
    """Query DrugBank API for drug metabolism information."""
    headers = {"Authorization": f"Bearer {DRUGBANK_API_KEY}"}
    base_url = "https://api.drugbank.com/v1"

    resp = requests.get(f"{base_url}/drugs/{drugbank_id}", headers=headers)
    resp.raise_for_status()
    drug = resp.json()

    return {
        "drug_name": drug.get("name"),
        "metabolizing_enzymes": [
            {"name": e.get("name"), "gene": e.get("gene_name")}
            for e in drug.get("enzymes", [])
        ],
        "transporters": [
            {"name": t.get("name"), "gene": t.get("gene_name")}
            for t in drug.get("transporters", [])
        ],
        "half_life": drug.get("half_life"),
    }


info = query_drugbank("DB00661")  # Verapamil
print(f"Drug: {info['drug_name']}, Half-life: {info['half_life']}")
for enz in info["metabolizing_enzymes"]:
    print(f"  Enzyme: {enz['name']} ({enz['gene']})")
```

## MetaboAnalystR: Differential Drug Response Analysis

**Goal:** Perform comprehensive differential analysis of metabolomics data comparing drug responders vs non-responders.

**Approach:** Use MetaboAnalystR to normalize, perform statistical testing, run pathway enrichment, and generate report-quality figures.

```r
library(MetaboAnalystR)

# Initialize and load data (rows = samples, cols = metabolites, first col = group)
mSet <- InitDataObjects("conc", "stat", FALSE)
mSet <- Read.TextData(mSet, "drug_response_metabolomics.csv", "rowu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)

# Normalize: log transform + auto-scaling
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "LogNorm", "AutoNorm", ref = NULL,
                       ratio = FALSE, ratioNum = 20)

# Volcano plot (fold change + t-test combined)
mSet <- FC.Anal(mSet, fc.thresh = 1.5, cmpType = 0)
mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05,
                     paired = FALSE, equal.var = TRUE)
mSet <- Volcano.Anal(mSet, FALSE, 1.5, 0, 0.1, TRUE, 0.05, TRUE, "raw")
mSet <- PlotVolcano(mSet, "volcano", format = "png", dpi = 300, width = NA)

# PLS-DA: classify responders vs non-responders
mSet <- PLSR.Anal(mSet, reg = TRUE)
mSet <- PlotPLS2DScore(mSet, "plsda_score", format = "png", dpi = 300,
                        width = NA, 1, 2, 0.95, 0)
mSet <- PlotPLS.Imp(mSet, "plsda_vip", format = "png", dpi = 300,
                     width = NA, "vip", "Comp. 1", 20)
```

## PK/PD Integration with Metabolomics

**Goal:** Correlate PK parameters (AUC, Cmax, Tmax) with metabolomics profiles to identify metabolic signatures of drug exposure.

```python
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

def pk_metabolomics_correlation(pk_file, metab_file):
    """Correlate PK parameters (AUC, Cmax, Tmax) with metabolite levels."""
    pk = pd.read_csv(pk_file, index_col="subject_id")
    metab = pd.read_csv(metab_file, index_col="subject_id")
    common = pk.index.intersection(metab.index)
    pk, metab = pk.loc[common], metab.loc[common]

    results = []
    for pk_param in pk.columns:
        for met_name in metab.columns:
            rho, pval = spearmanr(pk[pk_param], metab[met_name])
            results.append({"pk_parameter": pk_param, "metabolite": met_name,
                            "spearman_rho": round(rho, 3), "p_value": pval})

    df = pd.DataFrame(results)
    corrected = []
    for param in pk.columns:
        subset = df[df["pk_parameter"] == param].copy()
        _, fdr, _, _ = multipletests(subset["p_value"], method="fdr_bh")
        subset["fdr"] = fdr
        corrected.append(subset)
    return pd.concat(corrected).sort_values("p_value")


correlations = pk_metabolomics_correlation("pk_parameters.csv", "metabolomics.csv")
sig_corr = correlations[correlations["fdr"] < 0.05]
print(f"Significant PK-metabolite correlations: {len(sig_corr)}")
```

## Best Practices

- Always collect pre-dose baseline samples for paired comparisons
- Use at least 3 time points post-dose to capture metabolic dynamics
- Include BMI, age, sex, and CYP genotype as covariates
- Set mass tolerance based on instrument: TOF 10-20 ppm, Orbitrap 3-5 ppm
- Confirm putative metabolites with MS/MS fragmentation matching to parent drug
- Use internal standards for each metabolite class to correct for matrix effects

## Related Skills

- metabolomics-analysis/statistical-analysis - General differential analysis
- metabolomics-analysis/metabolite-annotation - Identify unknown metabolites
- metabolomics-analysis/xcms-preprocessing - LC-MS preprocessing
- ms-data-processing/matchms - Spectral matching for metabolite ID
