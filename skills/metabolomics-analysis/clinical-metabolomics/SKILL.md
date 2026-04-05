---
name: bio-metabolomics-clinical-reporting
description: Clinical metabolomics for inborn errors of metabolism (IEM) screening, newborn screening panels, and diagnostic reporting. Covers amino acid analysis, acylcarnitine profiling, organic acid analysis, z-score calculation against reference ranges, and clinical decision support with differential diagnosis. Use when interpreting clinical metabolomics results or building IEM screening workflows.
tool_type: python
primary_tool: pandas
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, scipy 1.12+, numpy 1.26+, matplotlib 3.8+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Clinical Metabolomics Reporting

**"Screen newborn dried blood spots for inborn errors of metabolism"** → Process amino acid and acylcarnitine panels from LC-MS/MS, calculate z-scores against age-matched reference ranges, flag abnormal results, and generate differential diagnoses.
- Python: pandas for data processing, scipy for z-scores, matplotlib for metabolite panels
- Analytical: ion-exchange chromatography, LC-MS/MS, GC-MS for organic acids

## Installation

```bash
pip install pandas numpy scipy matplotlib seaborn openpyxl
```

## Clinical Reference Ranges and Z-Score Calculation

**Goal:** Define age-matched reference ranges for newborn screening metabolites and compute z-scores for patient results.

**Approach:** Store reference ranges as structured data, compute z-scores using (value - mean) / SD, and flag results outside clinical thresholds.

```python
import pandas as pd
import numpy as np
from scipy import stats

reference_ranges = {
    'phenylalanine':         {'mean': 55.0,  'sd': 15.0, 'unit': 'umol/L', 'upper_cutoff': 120.0, 'lower_cutoff': 10.0},
    'tyrosine':              {'mean': 100.0, 'sd': 40.0, 'unit': 'umol/L', 'upper_cutoff': 250.0, 'lower_cutoff': 15.0},
    'leucine_isoleucine':    {'mean': 130.0, 'sd': 35.0, 'unit': 'umol/L', 'upper_cutoff': 250.0, 'lower_cutoff': 30.0},
    'methionine':            {'mean': 22.0,  'sd': 8.0,  'unit': 'umol/L', 'upper_cutoff': 50.0,  'lower_cutoff': 5.0},
    'citrulline':            {'mean': 18.0,  'sd': 7.0,  'unit': 'umol/L', 'upper_cutoff': 45.0,  'lower_cutoff': 3.0},
    'valine':                {'mean': 120.0, 'sd': 30.0, 'unit': 'umol/L', 'upper_cutoff': 250.0, 'lower_cutoff': 30.0},
    'C0_carnitine':          {'mean': 30.0,  'sd': 10.0, 'unit': 'umol/L', 'upper_cutoff': 60.0,  'lower_cutoff': 8.0},
    'C3_propionylcarnitine': {'mean': 2.0,   'sd': 0.8,  'unit': 'umol/L', 'upper_cutoff': 5.0,   'lower_cutoff': 0.2},
    'C5_isovalerylcarnitine':{'mean': 0.15,  'sd': 0.06, 'unit': 'umol/L', 'upper_cutoff': 0.45,  'lower_cutoff': 0.01},
    'C8_octanoylcarnitine':  {'mean': 0.10,  'sd': 0.05, 'unit': 'umol/L', 'upper_cutoff': 0.30,  'lower_cutoff': 0.01},
    'C14_1_tetradecenoylcarnitine': {'mean': 0.12, 'sd': 0.06, 'unit': 'umol/L', 'upper_cutoff': 0.60, 'lower_cutoff': 0.01},
    'C16_palmitoylcarnitine':{'mean': 3.5,   'sd': 1.2,  'unit': 'umol/L', 'upper_cutoff': 7.0,   'lower_cutoff': 0.5},
}

ref_df = pd.DataFrame(reference_ranges).T

def compute_z_scores(patient_results, ref_df):
    z_scores = {}
    for analyte, value in patient_results.items():
        if analyte in ref_df.index:
            ref = ref_df.loc[analyte]
            z = (value - ref['mean']) / ref['sd']
            z_scores[analyte] = {
                'value': value,
                'unit': ref['unit'],
                'z_score': round(z, 2),
                'percentile': round(stats.norm.cdf(z) * 100, 1),
                'flag': 'HIGH' if value > ref['upper_cutoff'] else ('LOW' if value < ref['lower_cutoff'] else 'NORMAL'),
            }
    return pd.DataFrame(z_scores).T

patient = {
    'phenylalanine': 280.0,
    'tyrosine': 45.0,
    'leucine_isoleucine': 140.0,
    'methionine': 25.0,
    'citrulline': 15.0,
    'C0_carnitine': 32.0,
    'C3_propionylcarnitine': 1.8,
    'C5_isovalerylcarnitine': 0.12,
}

report = compute_z_scores(patient, ref_df)
```

## Newborn Screening Amino Acid Panel

**Goal:** Process a batch of dried blood spot (DBS) amino acid results from the LC-MS/MS instrument and generate a screening report.

**Approach:** Read instrument export, apply internal standard correction, compute analyte ratios, and flag abnormal patterns.

```python
import pandas as pd
import numpy as np

raw_data = pd.read_csv('nbs_amino_acids_batch.csv')

def normalize_by_istd(raw_data, analyte_col, istd_col):
    return raw_data.assign(
        **{f'{analyte_col}_corrected': raw_data[analyte_col] / raw_data[istd_col]}
    )

istd_pairs = {'phenylalanine': 'phe_d5', 'tyrosine': 'tyr_d4',
              'leucine_isoleucine': 'leu_d3', 'methionine': 'met_d3',
              'citrulline': 'cit_d2', 'valine': 'val_d8'}

corrected = raw_data.copy()
for analyte, istd in istd_pairs.items():
    if analyte in corrected.columns and istd in corrected.columns:
        corrected[f'{analyte}_corrected'] = corrected[analyte] / corrected[istd]

corrected['phe_tyr_ratio'] = corrected['phenylalanine'] / corrected['tyrosine'].replace(0, np.nan)

screening_flags = pd.DataFrame({
    'sample_id': corrected['sample_id'],
    'phe_flag': corrected['phenylalanine'].apply(lambda x: 'ABNORMAL' if x > 120 else 'NORMAL'),
    'phe_tyr_ratio_flag': corrected['phe_tyr_ratio'].apply(lambda x: 'ABNORMAL' if x > 3.0 else 'NORMAL'),
    'cit_flag': corrected['citrulline'].apply(lambda x: 'HIGH' if x > 45 else ('LOW' if x < 3 else 'NORMAL')),
})
```

## Acylcarnitine Profiling

**Goal:** Analyze acylcarnitine profiles from LC-MS/MS for fatty acid oxidation disorders and organic acidemias.

**Approach:** Compute diagnostic ratios (C3/C2, C5/C2, C8/C10), classify acylcarnitine patterns, and match against known disorder profiles.

```python
import pandas as pd
import numpy as np

acylcarnitine_data = pd.read_csv('acylcarnitine_panel.csv', index_col='sample_id')

diagnostic_ratios = pd.DataFrame(index=acylcarnitine_data.index)
diagnostic_ratios['C3_C2'] = acylcarnitine_data['C3'] / acylcarnitine_data['C2'].replace(0, np.nan)
diagnostic_ratios['C5_C2'] = acylcarnitine_data['C5'] / acylcarnitine_data['C2'].replace(0, np.nan)
diagnostic_ratios['C5_C3'] = acylcarnitine_data['C5'] / acylcarnitine_data['C3'].replace(0, np.nan)
diagnostic_ratios['C8_C10'] = acylcarnitine_data['C8'] / acylcarnitine_data['C10'].replace(0, np.nan)
diagnostic_ratios['C14_1_C16'] = acylcarnitine_data['C14_1'] / acylcarnitine_data['C16'].replace(0, np.nan)
diagnostic_ratios['C0_C16_C18'] = acylcarnitine_data['C0'] / (
    acylcarnitine_data['C16'] + acylcarnitine_data['C18']
).replace(0, np.nan)

disorder_profiles = {
    'MCADD': {'C8': ('>',  0.30), 'C8_C10': ('>', 1.5), 'C6': ('>', 0.15)},
    'VLCADD': {'C14_1': ('>', 0.60), 'C14_1_C16': ('>', 0.1), 'C14': ('>', 0.70)},
    'PA': {'C3': ('>', 5.0), 'C3_C2': ('>', 0.25)},
    'MMA': {'C3': ('>', 5.0), 'C3_C2': ('>', 0.25)},
    'IVA': {'C5': ('>', 0.45), 'C5_C2': ('>', 0.03)},
    'GA1': {'C5DC': ('>', 0.15)},
    'CPT2': {'C16': ('>', 7.0), 'C18': ('>', 2.5), 'C0': ('<', 8.0)},
}

def screen_for_disorders(sample_values, ratios, profiles):
    matches = []
    combined = {**sample_values, **ratios}
    for disorder, criteria in profiles.items():
        criteria_met = 0
        total_criteria = len(criteria)
        for analyte, (op, threshold) in criteria.items():
            if analyte in combined:
                val = combined[analyte]
                if op == '>' and val > threshold:
                    criteria_met += 1
                elif op == '<' and val < threshold:
                    criteria_met += 1
        if criteria_met == total_criteria:
            matches.append({'disorder': disorder, 'confidence': 'HIGH'})
        elif criteria_met >= total_criteria * 0.5:
            matches.append({'disorder': disorder, 'confidence': 'POSSIBLE'})
    return matches
```

## Organic Acid Analysis (GC-MS Urine)

**Goal:** Process urine organic acid results from GC-MS and flag diagnostic markers for IEM.

**Approach:** Normalize to creatinine, compute z-scores, and identify patterns characteristic of specific organic acidurias.

```python
import pandas as pd
import numpy as np

organic_acids = pd.read_csv('urine_organic_acids.csv', index_col='sample_id')

creatinine = organic_acids.pop('creatinine_mmol_L')
oa_normalized = organic_acids.div(creatinine, axis=0)

# Reference ranges in mmol/mol creatinine (creatinine-normalized urine values)
oa_reference = {
    'methylmalonic_acid':      {'mean': 2.5, 'sd': 1.5,  'alert_threshold': 15.0},
    'methylcitric_acid':       {'mean': 1.0, 'sd': 0.8,  'alert_threshold': 5.0},
    '3_hydroxypropionic_acid': {'mean': 3.0, 'sd': 2.0,  'alert_threshold': 20.0},
    'glutaric_acid':           {'mean': 1.5, 'sd': 1.0,  'alert_threshold': 10.0},
    '3_hydroxyglutaric_acid':  {'mean': 2.0, 'sd': 1.2,  'alert_threshold': 8.0},
    'isovalerylglycine':       {'mean': 0.3, 'sd': 0.2,  'alert_threshold': 2.0},
    'suberylglycine':          {'mean': 0.1, 'sd': 0.08, 'alert_threshold': 1.0},
    'hexanoylglycine':         {'mean': 0.1, 'sd': 0.07, 'alert_threshold': 0.8},
    'orotic_acid':             {'mean': 1.0, 'sd': 0.6,  'alert_threshold': 5.0},
}

oa_ref_df = pd.DataFrame(oa_reference).T

def flag_organic_acids(sample_normalized, oa_ref_df):
    flags = {}
    for acid in sample_normalized.index:
        if acid in oa_ref_df.index:
            ref = oa_ref_df.loc[acid]
            val = sample_normalized[acid]
            z = (val - ref['mean']) / ref['sd']
            flags[acid] = {
                'value_per_creatinine': round(val, 2),
                'z_score': round(z, 2),
                'flag': 'CRITICAL' if val > ref['alert_threshold'] else ('HIGH' if z > 3 else 'NORMAL'),
            }
    return pd.DataFrame(flags).T
```

## Diagnostic Algorithm: Pattern Matching

**Goal:** Map flagged metabolite patterns to specific IEM diagnoses with differential diagnoses.

**Approach:** Encode diagnostic rules as structured decision trees and return ranked differentials with confirmatory test recommendations.

```python
import pandas as pd

diagnostic_rules = [
    {'condition': 'PKU', 'omim': 261600,
     'primary_markers': {'phenylalanine': ('>', 120)},
     'secondary_markers': {'phe_tyr_ratio': ('>', 3.0), 'tyrosine': ('<', 50)},
     'confirmatory': ['PAH gene sequencing', 'BH4 loading test']},
    {'condition': 'MSUD', 'omim': 248600,
     'primary_markers': {'leucine_isoleucine': ('>', 200)},
     'secondary_markers': {'valine': ('>', 250), 'alloisoleucine': ('>', 5)},
     'confirmatory': ['Alloisoleucine quantification', 'BCKDH enzyme assay']},
    {'condition': 'MMA', 'omim': 251000,
     'primary_markers': {'C3_propionylcarnitine': ('>', 5.0)},
     'secondary_markers': {'C3_C2': ('>', 0.25), 'methylmalonic_acid': ('>', 15)},
     'confirmatory': ['Urine organic acids', 'MUT/MMAA/MMAB gene panel', 'Vitamin B12 trial']},
    {'condition': 'IVA', 'omim': 243500,
     'primary_markers': {'C5_isovalerylcarnitine': ('>', 0.45)},
     'secondary_markers': {'isovalerylglycine': ('>', 2.0)},
     'confirmatory': ['Urine isovalerylglycine', 'IVD gene sequencing']},
    {'condition': 'MCADD', 'omim': 201450,
     'primary_markers': {'C8_octanoylcarnitine': ('>', 0.30)},
     'secondary_markers': {'C8_C10': ('>', 1.5), 'hexanoylglycine': ('>', 0.8)},
     'confirmatory': ['ACADM gene sequencing (c.985A>G common)', 'Urine acylglycines']},
    {'condition': 'Citrullinemia Type I', 'omim': 215700,
     'primary_markers': {'citrulline': ('>', 45)},
     'secondary_markers': {'argininosuccinate': ('<', 2.0)},
     'confirmatory': ['ASS1 gene sequencing', 'Plasma ammonia']},
]

def evaluate_diagnosis(patient_values, rules):
    results = []
    for rule in rules:
        primary_match = 0
        primary_total = len(rule['primary_markers'])
        for marker, (op, threshold) in rule['primary_markers'].items():
            if marker in patient_values:
                val = patient_values[marker]
                if (op == '>' and val > threshold) or (op == '<' and val < threshold):
                    primary_match += 1

        secondary_match = 0
        secondary_total = len(rule['secondary_markers'])
        for marker, (op, threshold) in rule['secondary_markers'].items():
            if marker in patient_values:
                val = patient_values[marker]
                if (op == '>' and val > threshold) or (op == '<' and val < threshold):
                    secondary_match += 1

        if primary_match == primary_total:
            confidence = 'LIKELY' if secondary_match > 0 else 'POSSIBLE'
            results.append({
                'condition': rule['condition'],
                'omim': rule['omim'],
                'confidence': confidence,
                'primary_markers_met': f'{primary_match}/{primary_total}',
                'secondary_markers_met': f'{secondary_match}/{secondary_total}',
                'confirmatory_tests': '; '.join(rule['confirmatory']),
            })
    return pd.DataFrame(results)

patient_values = {
    'phenylalanine': 280.0,
    'tyrosine': 45.0,
    'phe_tyr_ratio': 6.2,
    'leucine_isoleucine': 140.0,
    'C3_propionylcarnitine': 1.8,
    'C5_isovalerylcarnitine': 0.12,
    'C8_octanoylcarnitine': 0.08,
    'citrulline': 15.0,
}

diagnosis_report = evaluate_diagnosis(patient_values, diagnostic_rules)
```

## Quality Control: Internal Standards and External QC Programs

**Goal:** Monitor analytical quality using internal standard recovery, external QC program results (ERNDIM, CDC NSQAP), and Levey-Jennings tracking.

**Approach:** Compute IS recovery rates, compare to external QC target values, and generate control charts.

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

batch_istd = pd.read_csv('batch_istd_areas.csv', index_col='sample_id')
median_istd = batch_istd.median()
recovery = batch_istd.div(median_istd) * 100

istd_pass = (recovery > 70) & (recovery < 130)
failed_samples = recovery.index[~istd_pass.all(axis=1)].tolist()

nsqap_results = pd.DataFrame({
    'analyte': ['phenylalanine', 'C3', 'C8', 'C5', 'citrulline'],
    'lab_value': [62.0, 2.1, 0.12, 0.16, 19.5],
    'peer_mean': [61.5, 2.05, 0.115, 0.155, 19.8],
    'peer_sd': [5.0, 0.3, 0.02, 0.025, 2.5],
})
nsqap_results['z_score'] = (nsqap_results['lab_value'] - nsqap_results['peer_mean']) / nsqap_results['peer_sd']
nsqap_results['acceptable'] = nsqap_results['z_score'].abs() < 2.0

def levey_jennings_chart(qc_values, analyte_name, target, sd):
    runs = np.arange(1, len(qc_values) + 1)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(runs, qc_values, 'ko-', markersize=5)
    ax.axhline(y=target, color='green', linestyle='-', label='Target')
    for mult, color in [(2, 'orange'), (3, 'red')]:
        ax.axhline(y=target + mult * sd, color=color, linestyle='--')
        ax.axhline(y=target - mult * sd, color=color, linestyle='--')
    violations = [i for i, v in enumerate(qc_values) if abs(v - target) > 3 * sd]
    if violations:
        ax.plot([runs[i] for i in violations], [qc_values[i] for i in violations],
                'rv', markersize=10, label='Out of control')
    ax.set_xlabel('Run Number')
    ax.set_ylabel(f'{analyte_name} Concentration')
    ax.set_title(f'Levey-Jennings: {analyte_name}')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f'qc_chart_{analyte_name}.png', dpi=150)
    return fig

qc_phe = [58, 62, 55, 61, 59, 63, 57, 60, 64, 56, 61, 72, 59, 60, 58]
levey_jennings_chart(qc_phe, 'Phenylalanine', target=60.0, sd=5.0)
```

## Metabolite Panel Visualization

**Goal:** Generate a clinical-style metabolite panel plot showing patient values against reference ranges.

**Approach:** Plot each analyte as a point with reference range bars, color-coded by flag status.

```python
import matplotlib.pyplot as plt

def plot_metabolite_panel(report_df, title='Metabolite Screening Panel'):
    analytes = report_df.index.tolist()
    n = len(analytes)
    fig, ax = plt.subplots(figsize=(10, max(6, n * 0.4)))
    colors = {'NORMAL': '#4CAF50', 'HIGH': '#F44336', 'LOW': '#FF9800', 'CRITICAL': '#9C27B0'}
    for i, analyte in enumerate(analytes):
        row = report_df.loc[analyte]
        z, flag = row['z_score'], row['flag']
        ax.barh(i, z, color=colors.get(flag, '#757575'), height=0.6, alpha=0.7)
        ax.text(z + (0.1 if z >= 0 else -0.1), i, f"{row['value']} {row['unit']}",
                va='center', ha='left' if z >= 0 else 'right', fontsize=8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(analytes, fontsize=9)
    ax.axvline(x=0, color='black', linewidth=0.8)
    for threshold, color in [(2, 'orange'), (-2, 'orange'), (3, 'red'), (-3, 'red')]:
        ax.axvline(x=threshold, color=color, linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_xlabel('Z-Score')
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig('metabolite_panel.png', dpi=200)
    return fig

plot_metabolite_panel(report, title='Newborn Screening Panel - Patient 2024-001')
```

## Best Practices

- **Reference ranges:** Always use age-matched and method-specific reference ranges; newborn, pediatric, and adult cutoffs differ substantially for amino acids and acylcarnitines.
- **Second-tier testing:** Never report a positive IEM diagnosis from first-tier screening alone; always recommend confirmatory testing (enzyme assay, molecular genetics, repeat specimen).
- **Internal standards:** Use stable isotope-labeled internal standards for every analyte class; recovery should be 70-130% across the batch.
- **External QC:** Participate in proficiency testing programs (ERNDIM for Europe, CDC NSQAP for US) and track z-scores longitudinally.
- **Creatinine normalization:** Always normalize urine organic acids to creatinine; spot urine concentration varies by hydration status.
- **Ratio markers:** Use diagnostic ratios (Phe/Tyr, C3/C2, C5/C2, C8/C10) in addition to absolute values; ratios improve specificity over single-analyte cutoffs.
- **Westgard rules:** Apply multi-rule QC (1-3s, 2-2s, R-4s) rather than single-rule rejection to minimize false rejection of valid analytical runs.
- **Clinical context:** Metabolite elevations can be transient (prematurity, TPN, maternal medication); always correlate with clinical presentation before reporting.

## Related Skills

- targeted-analysis - Standard curves and absolute quantification
- normalization-qc - Batch QC and normalization methods
- statistical-analysis - Group comparisons and biomarker discovery
- pathway-mapping - Map abnormal metabolites to pathway context
