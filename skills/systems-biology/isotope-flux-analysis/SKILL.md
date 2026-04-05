---
name: isotope-flux-analysis
description: 13C metabolic flux analysis (13C-MFA) pipeline covering experimental design for isotope labeling, mass isotopomer distribution measurement, natural abundance correction, and flux estimation. Use for quantifying intracellular metabolic fluxes from 13C tracer experiments using LC-MS or GC-MS isotopomer data. Python-based correction with IsoCor; dedicated 13C-MFA solvers (INCA, OpenFLUX) are MATLAB-based.
tool_type: python
primary_tool: isocor
---

# 13C Metabolic Flux Analysis

## Overview

13C metabolic flux analysis (13C-MFA) quantifies intracellular metabolic fluxes by tracing the fate of 13C-labeled substrates through metabolic networks. The workflow involves feeding cells a 13C tracer (e.g., U-13C glucose), measuring mass isotopomer distributions (MIDs) of intracellular metabolites by MS, correcting for natural isotope abundance, and fitting a stoichiometric model to estimate fluxes.

## Installation

### Python Tools

```bash
uv pip install isocor cobra pandas numpy scipy matplotlib
```

### Verify Installation

```python
import isocor
from isocor import MetaboliteCorrectorFactory
import cobra
print(f"IsoCor version: {isocor.__version__}")
print(f"COBRApy version: {cobra.__version__}")
```

## Core Capabilities

### 1. Experimental Design for 13C Labeling

Tracer selection determines which fluxes are identifiable:

```python
# Common 13C tracers and their flux resolution capabilities
TRACER_DESIGNS = {
    "U-13C6-glucose": {
        "labeled_positions": [1, 2, 3, 4, 5, 6],
        "resolves": [
            "glycolysis vs pentose phosphate pathway splitting ratio",
            "TCA cycle flux (via citrate isotopomers)",
            "anaplerotic flux (via OAA/malate labeling)",
        ],
        "cost": "moderate",
        "typical_concentration_mm": 11.1,  # standard DMEM glucose
    },
    "1,2-13C2-glucose": {
        "labeled_positions": [1, 2],
        "resolves": [
            "oxidative vs non-oxidative PPP branch fluxes",
            "glycolysis vs PPP split with high resolution",
            "cycling within PPP (via C3 labeling patterns)",
        ],
        "cost": "high",
        "typical_concentration_mm": 11.1,
    },
    "U-13C5-glutamine": {
        "labeled_positions": [1, 2, 3, 4, 5],
        "resolves": [
            "reductive carboxylation vs oxidative TCA",
            "glutaminolysis flux",
            "malic enzyme flux (via pyruvate labeling from glutamine)",
        ],
        "cost": "moderate",
        "typical_concentration_mm": 2.0,
    },
}

def recommend_tracer(flux_of_interest):
    """Recommend a tracer based on the flux to be resolved."""
    recommendations = []
    for tracer, info in TRACER_DESIGNS.items():
        for capability in info["resolves"]:
            if flux_of_interest.lower() in capability.lower():
                recommendations.append({
                    "tracer": tracer,
                    "capability": capability,
                    "cost": info["cost"],
                })
    return recommendations
```

### 2. Mass Isotopomer Distribution (MID) Measurement

Extract MIDs from LC-MS or GC-MS data:

```python
import numpy as np
import pandas as pd

def extract_mid_from_areas(peak_areas):
    """Convert raw MS peak areas to fractional MID.

    Args:
        peak_areas: list or array of peak areas for M+0, M+1, ..., M+n

    Returns:
        Normalized MID (fractions summing to 1.0)
    """
    areas = np.array(peak_areas, dtype=np.float64)
    total = np.sum(areas)
    if total <= 0:
        raise ValueError("Total peak area must be positive")
    mid = areas / total
    return mid

# Example: citrate (6 carbons) measured by LC-MS in negative mode
# Peak areas for M+0 through M+6
citrate_areas = [15200, 8900, 42300, 18700, 55100, 12400, 89200]
citrate_mid = extract_mid_from_areas(citrate_areas)
print("Citrate MID:")
for i, fraction in enumerate(citrate_mid):
    print(f"  M+{i}: {fraction:.4f}")

# Quality check: verify enrichment above natural abundance
def check_enrichment(mid, n_carbons):
    """Verify that labeling is above natural abundance."""
    natural_m0 = (1 - 0.011) ** n_carbons
    if mid[0] > natural_m0 * 0.95:
        print(f"WARNING: M+0={mid[0]:.4f} near natural abundance ({natural_m0:.4f})")
    return mid[0] < natural_m0 * 0.95

check_enrichment(citrate_mid, n_carbons=6)
```

### 3. Natural Abundance Correction with IsoCor

IsoCor corrects measured MIDs for natural isotope abundance of C, H, N, O, Si, and S.
The `correct()` method returns a 4-tuple: `(corrected_area, isotopologue_fraction, residuum, mean_enrichment)`.

```python
from isocor import MetaboliteCorrectorFactory

# Define the metabolite and correction parameters
# Citrate: C6H8O7, measured as [M-H]- by LC-MS
corrector = MetaboliteCorrectorFactory(
    formula="C6H8O7",
    tracer="13C",
    resolution=60000,           # Orbitrap resolution at m/z 200
    mz_of_resolution=200,
    tracer_purity=[0.0, 0.99],  # 99% 13C purity of tracer
    correct_NA_tracer=True,
    derivative_formula="",       # No derivatization for LC-MS
)

# Measured (uncorrected) peak areas for M+0 through M+6
measured_areas = [15200, 8900, 42300, 18700, 55100, 12400, 89200]

# Correct for natural abundance
# Returns: corrected_area, isotopologue_fraction (MID), residuum, mean_enrichment
corrected_area, isotopologue_fraction, residuum, mean_enrichment = corrector.correct(
    measured_areas
)

print("Corrected MID (isotopologue fractions):")
for i, fraction in enumerate(isotopologue_fraction):
    print(f"  M+{i}: {fraction:.4f}")
print(f"Mean 13C enrichment: {mean_enrichment:.4f}")
print(f"Residuum: {residuum:.6f}")
```

#### Batch Correction for Multiple Metabolites

```python
METABOLITES = {
    "pyruvate": {"formula": "C3H4O3"}, "lactate": {"formula": "C3H6O3"},
    "citrate": {"formula": "C6H8O7"}, "alpha_ketoglutarate": {"formula": "C5H6O5"},
    "succinate": {"formula": "C4H6O4"}, "fumarate": {"formula": "C4H4O4"},
    "malate": {"formula": "C4H6O5"}, "aspartate": {"formula": "C4H7NO4"},
    "glutamate": {"formula": "C5H9NO4"},
}

def batch_correct_mids(measured_data, metabolites, resolution=60000,
                       tracer_purity=0.99):
    """Correct MIDs for a batch of metabolites.

    Args:
        measured_data: dict of {metabolite_name: raw peak areas list}
        metabolites: dict of {metabolite_name: {"formula": str}}
        resolution: MS resolution (use 1 for unit-mass GC-MS)
        tracer_purity: fraction of 13C in tracer (e.g., 0.99)

    Returns:
        dict with corrected MID, mean enrichment, and residuum per metabolite
    """
    corrected = {}
    for name, raw_areas in measured_data.items():
        if name not in metabolites:
            continue
        corrector = MetaboliteCorrectorFactory(
            formula=metabolites[name]["formula"], tracer="13C",
            resolution=resolution, mz_of_resolution=200,
            tracer_purity=[0.0, tracer_purity], correct_NA_tracer=True,
        )
        corrected_area, isotopologue_fraction, residuum, mean_enr = corrector.correct(
            raw_areas
        )
        corrected[name] = {
            "mid": list(isotopologue_fraction),
            "mean_enrichment": mean_enr,
            "residuum": residuum,
        }
    return corrected

measured = {
    "pyruvate": [45000, 12000, 18000, 25000],   # raw peak areas, not fractions
    "lactate": [42000, 14000, 19000, 25000],
    "citrate": [6000, 4000, 17000, 8000, 23000, 5000, 37000],
}
results = batch_correct_mids(measured, METABOLITES)
```

### 4. GC-MS MID Correction (TMS Derivatives)

GC-MS requires accounting for derivatization atoms:

```python
# For GC-MS with TMS derivatization, include derivative formula
# Example: alanine-2TMS has formula C3H7NO2, derivative adds C6H18Si2O
ala_corrector = MetaboliteCorrectorFactory(
    formula="C3H7NO2",
    tracer="13C",
    resolution=1,                # unit mass resolution for GC-MS
    mz_of_resolution=200,
    tracer_purity=[0.0, 0.99],
    correct_NA_tracer=True,
    derivative_formula="C6H18Si2O",  # 2x TMS groups
)

# Fragment-specific correction: use only the fragment formula
# Alanine m/z 116 fragment: [M-COOTMS]+, retains C2 backbone
ala_frag_corrector = MetaboliteCorrectorFactory(
    formula="C2H6NO",            # fragment contains 2 of 3 alanine carbons
    tracer="13C",
    resolution=1,
    mz_of_resolution=200,
    tracer_purity=[0.0, 0.99],
    correct_NA_tracer=True,
    derivative_formula="C3H9Si",  # 1x TMS on fragment
)
```

### 5. Stoichiometric Network Definition with COBRApy

Define the metabolic network stoichiometry that underlies the flux model. Note: COBRApy
performs FBA (linear programming on mass balances), which is fundamentally different from
13C-MFA (nonlinear fitting of isotopomer data). Here we use COBRApy only to define the
network structure; the actual 13C-MFA fitting uses the scipy-based optimizer in section 6.

```python
import cobra

model = cobra.Model("central_carbon")

# Define key metabolites
met_ids = ["g6p", "f6p", "gap", "pyr", "lac", "accoa", "oaa", "cit", "akg", "suc", "mal"]
mets = {m: cobra.Metabolite(m, compartment="c") for m in met_ids}
mets["glc_ext"] = cobra.Metabolite("glc_ext", compartment="e")

# Define reactions (simplified central carbon metabolism)
reaction_defs = [
    ("v_glc_uptake", {mets["glc_ext"]: -1, mets["g6p"]: 1}),
    ("v_glycolysis_upper", {mets["g6p"]: -1, mets["gap"]: 2}),
    ("v_glycolysis_lower", {mets["gap"]: -1, mets["pyr"]: 1}),
    ("v_ldh", {mets["pyr"]: -1, mets["lac"]: 1}),
    ("v_pdh", {mets["pyr"]: -1, mets["accoa"]: 1}),
    ("v_cs", {mets["accoa"]: -1, mets["oaa"]: -1, mets["cit"]: 1}),
    ("v_idh", {mets["cit"]: -1, mets["akg"]: 1}),
    ("v_akgdh_sdh", {mets["akg"]: -1, mets["suc"]: 1}),
    ("v_fum_mdh", {mets["suc"]: -1, mets["mal"]: 1}),
    ("v_mdh", {mets["mal"]: -1, mets["oaa"]: 1}),
]

for rxn_id, stoich in reaction_defs:
    rxn = cobra.Reaction(rxn_id)
    rxn.add_metabolites(stoich)
    rxn.bounds = (0, 1000)
    model.add_reactions([rxn])

print(f"Reactions: {len(model.reactions)}, Metabolites: {len(model.metabolites)}")
```

### 6. Flux Estimation Concepts

13C-MFA estimates fluxes by minimizing the difference between simulated and measured MIDs:

```python
from scipy.optimize import minimize
import numpy as np

def mid_residual(flux_vector, measured_mids, simulate_mids_fn, metabolite_names):
    """Objective function for flux fitting: sum of squared residuals.

    Args:
        flux_vector: free flux values to optimize
        measured_mids: dict of {metabolite: corrected MID array}
        simulate_mids_fn: function(fluxes) -> dict of simulated MIDs
        metabolite_names: list of metabolites to fit
    """
    simulated = simulate_mids_fn(flux_vector)
    ssr = 0.0
    for name in metabolite_names:
        if name in measured_mids and name in simulated:
            diff = np.array(measured_mids[name]) - np.array(simulated[name])
            ssr += np.sum(diff ** 2)
    return ssr

def estimate_fluxes(measured_mids, simulate_fn, n_free_fluxes,
                    n_restarts=50, bounds=None):
    """Multi-start optimization for flux estimation.

    Args:
        measured_mids: dict of corrected MIDs
        simulate_fn: MID simulation function
        n_free_fluxes: number of free flux parameters
        n_restarts: number of random restarts
        bounds: list of (min, max) tuples for each flux
    """
    if bounds is None:
        bounds = [(0.01, 100.0)] * n_free_fluxes

    best_result = None
    metabolite_names = list(measured_mids.keys())

    for _ in range(n_restarts):
        x0 = np.array([
            np.random.uniform(lo, hi) for lo, hi in bounds
        ])
        result = minimize(
            mid_residual,
            x0,
            args=(measured_mids, simulate_fn, metabolite_names),
            method="L-BFGS-B",
            bounds=bounds,
        )
        if best_result is None or result.fun < best_result.fun:
            best_result = result

    return best_result
```

### 7. Confidence Interval Estimation

Confidence intervals are computed by parameter continuation: fix each flux at test values across its range, re-optimize remaining fluxes, and find the boundary where SSR exceeds `best_SSR + chi2_threshold` (3.84 for 95% CI with 1 DOF).

```python
def compute_flux_ci(best_fluxes, measured_mids, simulate_fn, metabolite_names,
                    bounds, chi2_threshold=3.84, n_points=20):
    """Estimate 95% confidence intervals by parameter continuation."""
    best_ssr = mid_residual(best_fluxes, measured_mids, simulate_fn, metabolite_names)
    threshold = best_ssr + chi2_threshold
    cis = []

    for i in range(len(best_fluxes)):
        lo, hi = bounds[i]
        valid = []
        for test_val in np.linspace(lo, hi, n_points):
            reduced_bounds = [b for j, b in enumerate(bounds) if j != i]
            x0 = np.array([best_fluxes[j] for j in range(len(best_fluxes)) if j != i])
            obj = lambda x, tv=test_val: mid_residual(
                np.insert(x, i, tv), measured_mids, simulate_fn, metabolite_names)
            res = minimize(obj, x0, method="L-BFGS-B", bounds=reduced_bounds)
            if res.fun <= threshold:
                valid.append(test_val)
        cis.append((min(valid), max(valid)) if valid else (best_fluxes[i], best_fluxes[i]))
    return cis
```

### 8. INCA and OpenFLUX Workflow Overview

INCA (MATLAB-based, Vanderbilt, closed-source) is the most widely used dedicated 13C-MFA
tool. OpenFLUX is the primary open-source alternative (also MATLAB-based, uses the EMU
decomposition framework). There is currently no widely adopted pure-Python 13C-MFA solver
-- COBRApy does constraint-based FBA but NOT isotope-resolved MFA.

**Key distinction -- FBA vs 13C-MFA:**
- **FBA** uses linear programming on steady-state mass balances (no labeling data needed)
- **13C-MFA** fits mass isotopomer distribution (MID) data using nonlinear least-squares
  optimization, requiring atom (carbon) transition maps for every reaction

**Typical INCA/OpenFLUX workflow:**
1. Define metabolic network with atom (carbon) transition maps for each reaction
2. Import corrected MIDs and measured extracellular fluxes as constraints
3. Run multi-start flux estimation (50-100 restarts)
4. Compute 95% confidence intervals via parameter continuation (chi-square threshold)
5. Export flux maps and statistics

**Python-based partial alternative:** combine IsoCor (natural abundance correction) +
COBRApy (network stoichiometry definition) + scipy (nonlinear optimization) as shown
above. The EMU simulation step (forward labeling simulation from fluxes to predicted MIDs)
is the most complex part and is best handled by INCA or OpenFLUX for real datasets.

## Steady-State Verification

```python
def verify_isotopic_steady_state(mid_timepoints, tolerance=0.02):
    """Check that MIDs have reached isotopic steady state.

    Args:
        mid_timepoints: dict of {time_hours: MID_array}
        tolerance: maximum allowed change between consecutive timepoints

    Returns:
        True if steady state is reached
    """
    times = sorted(mid_timepoints.keys())
    if len(times) < 2:
        raise ValueError("Need at least 2 timepoints to verify steady state")

    for i in range(1, len(times)):
        prev_mid = np.array(mid_timepoints[times[i - 1]])
        curr_mid = np.array(mid_timepoints[times[i]])
        max_change = np.max(np.abs(curr_mid - prev_mid))
        if max_change > tolerance:
            print(f"NOT at steady state: max MID change = {max_change:.4f} "
                  f"between {times[i-1]}h and {times[i]}h")
            return False

    print(f"Steady state reached (max change < {tolerance})")
    return True

# Check with 3 timepoints
mid_over_time = {
    12: [0.44, 0.13, 0.18, 0.25],
    18: [0.43, 0.13, 0.19, 0.25],
    24: [0.42, 0.14, 0.19, 0.25],
}
verify_isotopic_steady_state(mid_over_time)
```

## Best Practices

- **Verify isotopic steady state** by sampling at least 2-3 timepoints before harvesting
- **Use tracer purity correction** (typically 99% 13C, account for the 1% 12C contamination)
- **Multi-start optimization** with at least 50 random restarts to avoid local minima
- **Report 95% confidence intervals** for all estimated fluxes, not just point estimates
- **Measure extracellular fluxes** (glucose uptake, lactate secretion) as independent constraints
- **Fragment-specific correction** for GC-MS: correct only the carbons retained in the measured fragment
- **MS resolution matters**: high-resolution MS (Orbitrap) resolves isobaric interferences that unit-mass GC-MS cannot
- **1,2-13C glucose** provides better PPP flux resolution than U-13C glucose
- **Natural abundance experiments** (unlabeled) should be run as controls to verify correction accuracy
- **Biological replicates** (n>=3) are essential; technical replicates alone are insufficient

## Resources

- IsoCor documentation: https://isocor.readthedocs.io
- COBRApy documentation: https://cobrapy.readthedocs.io
- INCA (Vanderbilt): https://mfa.vueinnovations.com
- OpenFLUX: https://github.com/OpenFLUX
- Antoniewicz 2018 (13C-MFA review): doi:10.1016/j.copbio.2017.11.004
- Zamboni et al. 2009 (MFA methods): doi:10.1016/j.copbio.2009.01.001
