# Phase 1: Outcome Analysis — Instruction

This notebook builds **Phase 1** of the PBTA_RNA Clinical Deep Dive. It tests associations between clinical variables and patient outcomes (OS and EFS) using Mann-Whitney / Chi-squared (binary STATUS) + KM + log-rank (time-to-event).

**Prerequisite:** Read `/home/alon/menow_home_ass/context/clinical_deep_dive_general.md` for shared methodology (test selection, FDR, effect sizes, formatting, imports, helper functions). This instruction only covers Phase 1 specifics.

---

## Scope

| Section | Variables | Outcome | Approach |
|---------|-----------|---------|----------|
| 1A | AGE | OS & EFS | Mann-Whitney (binary STATUS) + KM + log-rank (time) |
| 1B | TUMOR_FRACTION | OS & EFS | Mann-Whitney (binary STATUS) + KM + log-rank (time) |
| 1C | TUMOR_PLOIDY | OS & EFS | Mann-Whitney (binary STATUS) + KM + log-rank (time) |
| 1D | SEX | OS & EFS | Chi-squared (binary STATUS) + KM + log-rank (time) |
| 1E | CANCER_PREDISPOSITIONS | OS & EFS | KM + log-rank (time) |

**Threshold:** Cancer groups with n >= 20 samples (with the relevant non-missing outcome data) are included in per-group analyses.

---

## Notebook Structure

### 0. Setup Cells

A single cell at the top with:
1. All imports (from general instruction Section 6)
2. All helper functions (from general instruction Section 7)
3. Data loading and cleaning
4. Merge patients -> samples on PATIENT_ID
5. Print merge validation (orphan check)
6. Define and print target cancer groups for per-group analysis (n >= 20 samples with OS data)

```python
# After loading merged
MIN_SAMPLES = 20
# Count samples with non-missing OS data per group
os_available = merged.dropna(subset=['OS_MONTHS', 'os_event'])
per_group_counts = os_available.groupby('CANCER_GROUP').size()
target_groups = per_group_counts[per_group_counts >= MIN_SAMPLES].index.tolist()
print(f"Groups with n >= {MIN_SAMPLES} for OS analysis: {len(target_groups)}")
for g in target_groups:
    print(f"  {g}: {per_group_counts[g]} samples")
excluded = per_group_counts[per_group_counts < MIN_SAMPLES]
if len(excluded) > 0:
    print(f"Excluded groups (< {MIN_SAMPLES}):")
    for g in excluded.index:
        print(f"  {g}: {per_group_counts[g]} samples")
```

7. Initialize the summary results table as an empty list:
```python
results_rows = []
```

### 1. Per-Section Pattern

Every section (1A-1E) follows this structure:

#### Markdown Cell
```
## Section 1X: [Variable Name] x OS/EFS

**Rationale:** [1-2 sentences about why this comparison matters clinically]
```

#### For continuous variables (AGE, TF, TP) — Sections 1A, 1B, 1C:

**1X.1 Global — Binary (OS_STATUS)**
- Split patients into LIVING vs DECEASED
- Test: Mann-Whitney U
- Effect size: Cliff's delta
- Plot: Boxplot (LIVING vs DECEASED) with jittered points
  - Annotate with: `U = X, p = X.XXXX, d = X.XX`
- Record to summary table

**1X.2 Global — Binary (EFS_STATUS)**
- Same structure with EFS (Event vs No Event)
- Record to summary table

**1X.3 Global — Time-to-Event (OS_MONTHS)**
- Dichotomize the variable at the median (high vs low)
- KM curves with log-rank test
- Risk table below the KM curve
- Record to summary table

**1X.4 Global — Time-to-Event (EFS_MONTHS)**
- Same as 1X.3 but for EFS

**1X.5 Per-Group — KM subplots (OS)**
- For each cancer group (n ≥ 20):
  - Dichotomize variable at the group-specific median
  - KM curve with log-rank test
- Arange in a subplot grid (e.g., 3×3 for 9 groups, or 2×3 for top 6)
- Each subplot: KM curve + log-rank p-value annotation
- Apply FDR correction across all per-group log-rank p-values
- Print a summary table with: Cancer Group, N, N events, median survival (high), median survival (low), log-rank p, FDR q

**1X.6 Per-Group — KM subplots (EFS)**
- Same for EFS

#### For categorical variables (SEX) - Section 1D:

**1D.1 Global - Binary (OS_STATUS)**
- Contingency table: SEX (Male/Female) x OS_STATUS (LIVING/DECEASED)
- Test: Chi-squared
- Effect size: Cramer's V
- Plot: Stacked bar chart with chi-squared annotation
- Record to summary table

**1D.2 Global - Binary (EFS_STATUS)**
- Same with EFS

**1D.3 Global - Time-to-Event (OS_MONTHS)**
- KM curves: stratified by SEX
- Log-rank p-value
- Risk table

**1D.4 Global - Time-to-Event (EFS_MONTHS)**
- Same for EFS

**1D.5 Per-Group — KM subplots (OS by SEX)**
- For each cancer group (n ≥ 20):
  - KM curves stratified by SEX
  - Log-rank test
- Arange in a subplot grid (2×3 for top 6 groups)
- Each subplot: KM curves (Male vs Female) + log-rank p
- Apply FDR correction across all per-group log-rank p-values

**1D.6 Per-Group — KM subplots (EFS by SEX)**
- Same for EFS

#### For predisposition - Section 1E:

**1E.1 Binary Predisposition x OS (Time-to-Event)**
- Create binary flag: `has_predisposition`
- KM curves: predisposition vs no predisposition
- Log-rank test
- Record to summary table

**1E.2 Binary Predisposition x EFS (Time-to-Event)**
- Same for EFS

**1E.3 Per-Type Predisposition x OS (Time-to-Event)**
- Top 5 predisposition types + "Other" + "No predisposition"
- Multi-group KM curve
- Global log-rank test
- Pairwise log-rank with FDR correction

**1E.4 Per-Type Predisposition x EFS**
- Same for EFS

---

## 2. Per-Section Detail

### Section 1A: AGE x OS/EFS

**Rationale:** Age at diagnosis is a known prognostic factor in many pediatric brain tumors. Younger patients may have different tumor biology and different outcomes.

**Variables:**
- Predictor: `AGE` (continuous, years)
- Binary outcome: `OS_STATUS` (LIVING/DECEASED), `EFS_STATUS` (binary: Event/No Event)
- Time-to-event: `OS_MONTHS` + `os_event`, `EFS_MONTHS` + `efs_event`

**Plots to produce:**
- P1A.1: Boxplot — AGE × OS_STATUS (global) with Mann-Whitney p + Cliff's d
- P1A.2: Boxplot — AGE × EFS_STATUS (global) with Mann-Whitney p + Cliff's d
- P1A.3: KM — OS by AGE group (young ≤ median, old > median) with log-rank p + risk table
- P1A.4: KM — EFS by AGE group with log-rank p + risk table
- P1A.5: KM subplot grid — OS by AGE group per cancer group + FDR q-values
- P1A.6: KM subplot grid — EFS by AGE group per cancer group + FDR q-values

**Summary table entries:**
- "AGE × OS_STATUS (global)" | Mann-Whitney
- "AGE × OS (KM, global)" | Log-rank
- "AGE × OS (KM, per-group)" | Log-rank (one row per group with FDR)
- ... same for EFS

---

### Section 1B: TUMOR_FRACTION x OS/EFS

**Rationale:** Tumor purity (fraction of tumor cells in the sample) may reflect tumor biology - more aggressive tumors may outgrow their stroma (higher purity), or lower purity may indicate a more desmoplastic/infiltrative phenotype.

**Variables:**
- Predictor: `TUMOR_FRACTION` (continuous, 0-1)
- Same outcome structure as 1A

**Plots to produce:**
- P1B.1: Boxplot — TF × OS_STATUS (global) with MW p + Cliff's d
- P1B.2: Boxplot — TF × EFS_STATUS (global) with MW p + Cliff's d
- P1B.3: KM — OS by TF group (high > median, low ≤ median) + log-rank
- P1B.4: KM — EFS by TF group + log-rank
- P1B.5: KM subplot grid — OS by TF group per cancer group + FDR q
- P1B.6: KM subplot grid — EFS by TF group per cancer group + FDR q

**Validation:**
```python
tf_os = merged.dropna(subset=['TUMOR_FRACTION', 'OS_MONTHS', 'os_event'])
print(f"TF x OS: {len(tf_os)} samples with complete data ({len(tf_os)/len(merged.dropna(subset=['OS_MONTHS','os_event']))*100:.1f}% of OS samples)")
```

**Basic notebook reference:** Step 16 shows TF by CANCER_GROUP and TUMOR_TYPE. This extends to outcome.

---

### Section 1C: TUMOR_PLOIDY x OS/EFS

**Rationale:** Aneuploidy (ploidy != 2) is a hallmark of cancer. Higher ploidy may indicate more genomic instability and potentially worse prognosis.

**Variables:**
- Predictor: `TUMOR_PLOIDY` (continuous)
- Same outcome structure

**Plots to produce:**
- P1C.1: Boxplot — TP × OS_STATUS (global) with MW p + Cliff's d
- P1C.2: Boxplot — TP × EFS_STATUS (global) with MW p + Cliff's d
- P1C.3: KM — OS by TP group (diploid ≈ 2 vs aneuploid ≠ 2) + log-rank
- P1C.4: KM — EFS by TP group + log-rank
- P1C.5: KM subplot grid — OS by TP group per cancer group + FDR q
- P1C.6: KM subplot grid — EFS by TP group per cancer group + FDR q

**Validation:**
```python
tp_os = merged.dropna(subset=['TUMOR_PLOIDY', 'OS_MONTHS', 'os_event'])
print(f"TP x OS: {len(tp_os)} samples")
```

---

### Section 1D: SEX x OS/EFS

**Rationale:** Sex differences in cancer incidence and outcome are well-documented in adults. Less is known in pediatric brain tumors.

**Variables:**
- Predictor: `SEX` (categorical: Male/Female)
- Same outcome structure

**Plots to produce:**
- P1D.1: Stacked bar — SEX × OS_STATUS (global) with χ² p + Cramer's V
- P1D.2: Stacked bar — SEX × EFS_STATUS (global) with χ² p + Cramer's V
- P1D.3: KM — OS by SEX + log-rank
- P1D.4: KM — EFS by SEX + log-rank
- P1D.5: KM subplot grid — OS by SEX per cancer group + FDR q
- P1D.6: KM subplot grid — EFS by SEX per cancer group + FDR q

**Validation:**
```python
sex_os = merged.dropna(subset=['SEX', 'OS_STATUS'])
ct = pd.crosstab(sex_os['SEX'], sex_os['OS_STATUS'].str.replace(r'^\d+:', '', regex=True))
print(f"SEX x OS contingency table:\n{ct}")
# Check for small expected frequencies
chi2, p, dof, expected = chi2_contingency(ct)
print(f"Min expected: {expected.min():.1f} {'-> using Fisher exact' if (expected < 5).any() else '-> chi-squared OK'}")
```

**Basic notebook reference:** Step 15 tests SEX x CANCER_GROUP. This extends to outcome.

---

### Section 1E: PREDISPOSITION x OS/EFS

**Rationale:** Patients with known cancer predisposition syndromes may have different tumor biology and treatment response.

**Variables:**
- Predictor: `CANCER_PREDISPOSITIONS`
  - Binary: has_predisposition (True/False)
  - Multi-type: top individual syndromes + "Other" + "No predisposition" + "Unknown"

**Plots to produce:**
- P1E.1: KM - OS by predisposition status (binary) + log-rank
- P1E.2: KM - EFS by predisposition status (binary) + log-rank
- P1E.3: KM - OS by predisposition type (top 5 + Other + None) + multi-group log-rank
- P1E.4: KM - EFS by predisposition type + multi-group log-rank
- P1E.5: Bar chart - predisposition type frequency with outcome composition

**Validation:**
```python
pred_os = merged.dropna(subset=['OS_MONTHS', 'os_event'])
pred_os['has_pred'] = pred_os['CANCER_PREDISPOSITIONS'].isin(['No predisposition', 'Unknown'])
print(f"With predisposition: {(~pred_os['has_pred']).sum()} / {len(pred_os)}")
```

**Basic notebook reference:** Step 17 shows predisposition x OS (global binary). This extends with per-type breakdown and EFS.

---

## 3. Summary Table

At the end of the notebook, build the results DataFrame and display:

```python
summary = pd.DataFrame(results_rows)
summary = summary.sort_values(['Phase', 'p_value'])
print(f"\n{'='*60}")
print(f"Phase 1 Results Summary: {len(summary)} statistical tests")
print(f"{'='*60}")
display(summary)

# Counts
total = len(summary)
sig = summary['Significant'].str.contains('✅').sum()
sig_fdr = (summary['FDR_BH'] < 0.05).sum()
print(f"\nTotal tests: {total}")
print(f"Significant (raw p < 0.05): {sig} ({sig/total*100:.1f}%)")
print(f"Significant (FDR < 0.05): {sig_fdr} ({sig_fdr/total*100:.1f}%)")

# Save
summary.to_csv('/home/alon/menow_home_ass/notebooks/phase1_deep_dive/phase1_results.csv', index=False)
print("\nSaved: notebooks/phase1_deep_dive/phase1_results.csv")
```

---

## 4. Phase 1 Specific Checklist

- [ ] Both binary (STATUS) AND time-to-event (MONTHS + event) analyses applied for each variable
- [ ] Per-group KM subplot grids include log-rank p-values + FDR q-values
- [ ] KM curves include risk tables + log-rank p-value
- [ ] Effect sizes reported (Cliff's d, Cramer's V) alongside p-values
- [ ] n >= 20 threshold enforced for per-group analyses
- [ ] Summary table populated with all tests
- [ ] CSV summary saved
- [ ] Basic notebook steps referenced where applicable
