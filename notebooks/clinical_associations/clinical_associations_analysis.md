# Phase 2+3: Clinical Associations — Cross-Categorical & Numeric Comparisons

**Prerequisite:** Read `/home/alon/menow_home_ass/context/clinical_deep_dive_general.md` for shared methodology (test selection, FDR, effect sizes, formatting, imports, helper functions).

This notebook builds **Phase 2 (Cross-Categorical)** and **Phase 3 (Numeric Comparisons)** of the PBTA_RNA Clinical Deep Dive. Phase 2 tests per-group enrichment of categorical variables; Phase 3 compares numeric distributions across groups and computes correlations.

---

## Scope

### Phase 2 — Cross-Categorical (per-group, n ≥ 20)

| # | Test | Method | Per-Group | FDR Family |
|---|------|--------|-----------|------------|
| 1 | SEX enrichment | Binomial vs 50:50 | Each group tested | "SEX enrichment" |
| 2 | Predisposition enrichment | Fisher exact 2×2 (has_pred × in_group) | Each pred category × each group | "Predisposition enrichment" |
| 3 | SUBTYPE distribution | Cross-tab heatmap only | Descriptive — no tests | — |
| 4 | SEX × PRED 3-way | Fisher / Chi² 2×2 (Sex × PRED_BINARY) | Each group | "SEX × Predisposition" |
| 5 | Race enrichment | Fisher exact 2×2 (has_race × in_group) | Each race category × each group | "Race enrichment" |

**Threshold:** Only cancer groups with n ≥ 20 samples (with non-missing data for the specific test) are included.

### Phase 3 — Numeric Comparisons

| # | Test | Method | Per-Group | FDR Family |
|---|------|--------|-----------|------------|
| 6 | AGE × CANCER_GROUP | Global KW + per-group MW | Each group vs all others | "AGE × CG" |
| 7 | TF × CANCER_GROUP | Global KW + per-group MW | Each group vs all others | "TF × CG" |
| 8 | TP × CANCER_GROUP | Global KW + per-group MW | Each group vs all others | "TP × CG" |
| 9 | AGE × TF | Spearman ρ + LOESS scatter | Global | "Numeric correlations" |
| 10 | AGE × TP | Spearman ρ + LOESS scatter | Global | "Numeric correlations" |
| 11 | TF × TP | Spearman ρ + LOESS scatter | Global | "Numeric correlations" |

---

## Notebook Structure

### Cell 0 (markdown): Title
### Cell 1 (markdown): Summary — two paragraphs describing Phase 2 and Phase 3
### Cell 2 (code): Imports
### Cell 3 (code): Data loading — read, merge, clean, print sample sizes
### Cell 4 (markdown): `## Phase 2: Cross-Categorical Associations`
### Cell 5 (code): Phase 2 helpers — `sig_star`, `run_enrichment`
### Cell 6 (code): Test 1 — SEX enrichment per group (binomial vs 50:50)
### Cell 7 (code): Test 2 — Predisposition enrichment per group (Fisher exact)
### Cell 8 (code): Test 3 — SUBTYPE descriptive heatmap (no p-values)
### Cell 9 (code): Test 4 — SEX × PRED 3-way per group
### Cell 10 (code): Test 5 — Race enrichment per group (Fisher exact)
### Cell 11 (code): Phase 2 results — FDR per Comparison, display
### Cell 12 (markdown): `## Phase 3: Numeric Comparisons`
### Cell 13 (code): Phase 3 helpers — `kw_per_group`
### Cell 14 (code): Test 6 — AGE × CANCER_GROUP
### Cell 15 (code): Test 7 — TF × CANCER_GROUP
### Cell 16 (code): Test 8 — TP × CANCER_GROUP
### Cell 17 (code): Tests 9-11 — Spearman correlations
### Cell 18 (code): Phase 3 results — FDR per Comparison, display
### Cell 19 (code): Combined results + save CSV
### Cell 20 (markdown): Summary markdown

---

## Test Details

### Test 1 — SEX enrichment per group (binomial vs 50:50)

For each cancer group (n ≥ 20):
- Count Male and Female
- Binomial test against null proportion p = 0.5 (two-sided)
- Report: Group, N, N_Male, N_Female, %Male, %Female, p_value, enriched sex
- FDR family: "SEX enrichment"
- Plot: Horizontal bar chart — %Male and %Female per group, color-coded by significance

### Test 2 — Predisposition enrichment per group (Fisher exact)

For each unique value in `CANCER_PREDISPOSITIONS` (after cleaning):
- For each cancer group (n ≥ 20):
  - Build 2×2: `has_this_pred` (yes/no) × `in_this_group` (yes/no)
  - Fisher exact test
- Report: Predisposition, Group, N_group, N_with_pred_in_group, N_with_pred_total, odds_ratio, p_value
- FDR family: "Predisposition enrichment"
- Plot: Heatmap (predisposition × cancer group) with count annotations

### Test 3 — SUBTYPE descriptive only

- Cross-tabulation: MOLECULAR_SUBTYPE × CANCER_GROUP
- Heatmap with count annotations (text_auto=True)
- Top 30 subtypes for readability
- No p-values, no statistical tests

### Test 4 — SEX × PRED 3-way per group

- Create `PRED_BINARY`: "Any predisposition" / "No predisposition"
- For each cancer group (n ≥ 20):
  - 2×2 contingency: SEX (M/F) × PRED_BINARY
  - Fisher exact if any expected cell < 5, else Chi² with Yates correction
- Report: Group, N, p_value, odds ratio (or Cramer's V for Chi²)
- FDR family: "SEX × Predisposition"
- Plot: Side-by-side %Female bars for predisposed vs non-predisposed per group

### Test 5 — Race enrichment per group (Fisher exact)

For each unique value in `RACE` (after cleaning, excluding Unknown):
- For each cancer group (n ≥ 20):
  - Build 2×2: `has_this_race` (yes/no) × `in_this_group` (yes/no)
  - Fisher exact test
- Report: Race, Group, N, N_with_race, p_value, odds_ratio
- FDR family: "Race enrichment"
- Plot: Heatmap (race × cancer group) with count annotations

### Tests 6-8 — Numeric × CANCER_GROUP (KW + MW)

- Global Kruskal-Wallis test across all groups
- Per-group Mann-Whitney (each group vs all others)
- Effect size: ε² (KW) and rank-biserial r (MW)
- Boxplot ordered by median
- FDR family per variable: "AGE × CG", "TF × CG", "TP × CG"

### Tests 9-11 — Numeric × Numeric (Spearman)

- Spearman rank correlation (ρ)
- Scatter plot with LOESS trendline (`trendline='lowess'`)
- FDR family: "Numeric correlations"

---

## CSV Output Format

```
Phase,Comparison,Test,Group,N,N_events,Statistic,p_value,FDR_WithinFamily,Significant,Effect_Size
```

- `FDR_WithinFamily`: Benjamini-Hochberg FDR applied within each Comparison family
- `Significant`: `*` p<0.05, `**` p<0.01, `***` p<0.001, `ns` not significant
- Saved to: `notebooks/clinical_associations/clinical_associations_results.csv`

---

## FDR Strategy

FDR is applied **within each Comparison family** (i.e., within each row of the Scope table above). For example, all "SEX enrichment" tests share one FDR correction; all "Predisposition enrichment" tests share another. This controls the FDR per conceptual hypothesis family.

---

## Visualization Rules

- All plots use Plotly (plotly.express and plotly.graph_objects)
- All `fig.show()` calls wrapped in `try/except` blocks
- Heatmaps use `text_auto=True` to show actual counts in cells
- No `.background_gradient()` or `.to_html()` on pandas StyleFrames

---

## Output Files

| File | Path |
|------|------|
| Source notebook | `notebooks/clinical_associations/clinical_associations_analysis.ipynb` |
| Executed notebook | `clinical_associations_analysis_executed.ipynb` (root) |
| Results CSV | `notebooks/clinical_associations/clinical_associations_results.csv` |
| Build script | `notebooks/clinical_associations/src/build_nb.py` |
