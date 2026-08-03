# Clinical Deep Dive — Master Plan

This document is the **high-level overview** of the PBTA_RNA clinical deep-dive analysis. It defines the scope, file structure, and cross-references between all components. Detailed methodology and implementation are in the files it references.

---

## File Structure

```
context/
├── clinical_deep_dive_general.md     ← Shared methodology, conventions, imports, helpers
└── clinical_deep_dive_plan.md         ← THIS FILE — high-level reference

notebooks/
├── survival_analysis_executed.ipynb  ← Phase 1: executed notebook (has outputs)
├── age_deciles_analysis_executed.ipynb ← Age decile analysis (has outputs)
├── clinical_associations_analysis_executed.ipynb ← Phase 2+3: executed (has outputs)
├── clinical_analysis_excuted.ipynb   ← Basic exploratory notebook
├── survival_analysis/
│   ├── survival_analysis.ipynb       ← Source (no outputs)
│   ├── survival_analysis.md          ← Instruction file
│   ├── survival_analysis_results.csv
│   └── fdr_significant_summary.csv
├── age_deciles/
│   ├── age_deciles_analysis.ipynb    ← Source (no outputs)
│   ├── age_deciles_results.csv
│   └── src/
│       └── build_nb.py
├── clinical_associations/
│   ├── clinical_associations_analysis.ipynb  ← Source (no outputs)
│   ├── clinical_associations_analysis.md     ← Instruction file
│   ├── clinical_associations_results.csv
│   ├── clinical_associations_summary.csv
│   └── src/
│       └── build_nb.py
└── clinical_analysis/
    └── src/
        └── imports.py                 ← Shared helper functions
```

---

## How to Use These Documents

1. **Start here** for the big picture
2. **Read `clinical_deep_dive_general.md`** for shared statistical methodology, formatting rules, imports, and helper functions — this applies to all phases
3. **Read the relevant phase file** (e.g., `notebooks/clinical_associations_analysis.md`) for the specific implementation of that phase's analyses

---

## Phases Overview

| Phase | Title | Variables Tested | Outcome | Statistical Tests | Threshold | Status |
|-------|-------|-----------------|---------|------------------|-----------|--------|
| **1** | Outcome Analysis | AGE, TF, TP, SEX, PREDISPOSITION | OS & EFS (binary + time) | MW, χ², KM + log-rank, Cliff's d, Cramer's V | n ≥ 20 | ✅ Complete |
| **2 + 3** | Clinical Associations | SEX, PRED, SUBTYPE, RACE, AGE, TF, TP — per-group enrichment & distribution | Association only | Binomial, Chi-squared, KW + MW, Spearman ρ, FDR per-family | n ≥ 20 | 🔲 Plan review |
| **4** | Multivariate Models | All above combined | OS & EFS (time) | Cox PH, HR forest | n ≥ 20 | 🔲 Plan review |
| **5** | Unsupervised | AGE, TF, TP, SEX, CG, SUBTYPE, etc. | Cluster discovery | PCA, t-SNE, K-means, FAMD, silhouette, log-rank validation | n ≥ 50 | 🔲 Plan review |
| **6** | Summary | All above | Consolidated report | — | — | 🔲 Planned |

---

## Key Design Decisions

### Both Binary and Time-to-Event

Every outcome analysis uses **two complementary approaches**:

| Approach | What it tests | Statistical Tool |
|----------|--------------|------------------|
| **Binary (STATUS)** | "Do patients with different outcomes differ in this variable?" | Mann-Whitney U (continuous predictor), Chi-squared (categorical predictor) |
| **Time-to-event (MONTHS + event)** | "Does this variable predict the *timing* of the outcome?" | Cox PH (continuous predictor), KM + log-rank (categorical/dichotomized predictor) |

The binary approach is simpler and more interpretable for clinicians. The time-to-event approach is statistically more powerful and clinically more relevant. Reporting both allows cross-validation — if they agree, the finding is robust. If they disagree, it signals something interesting (e.g., a variable affects early but not late survival).

### Sample Size Thresholds

| Phase | Minimum per group |
|-------|-------------------|
| 1–4 (statistical tests) | **n ≥ 20** |
| 5 (unsupervised) | **n ≥ 50** |

### FDR Correction — Notes for Revisit After All Phases

**Current approach (Phase 1):** FDR (Benjamini-Hochberg) is applied **within each cancer group** — all tests conducted on a single group (AGE-OS, AGE-EFS, TF-OS, TF-EFS, TP-OS, TP-EFS, SEX-OS, SEX-EFS) form one FDR family. Tests on different groups are treated as independent families.

**The open question (deferred):** After all phases are complete, should we define FDR families differently?

Key considerations:
- Phase 1 SEX × OS within a group and Phase 2 SEX × GROUP both test SEX as a variable — they are not fully independent.
- Reasonable alternative: define FDR families **per clinical variable across all phases** (one SEX family covering all SEX tests from Phase 1, 2, etc.; one AGE family covering Phase 1 AGE tests + decile tests + Phase 3 AGE × CG, etc.).
- This would be more conservative but more principled.
- **Decision**: Keep per-phase/per-group FDR for now. Revisit at the end of Phase 6 to apply a unified cross-phase FDR as a sensitivity analysis. Document how results change.

### Effect Sizes

Reported alongside every p-value (Cliff's d, ε², Cramer's V, HR). See general instruction for interpretation thresholds.

### Phase 2+3 Design Details

#### Phase 2 — Cross-Categorical (per-group enrichment tests)

All Phase 2 tests produce **one row per group** in the results CSV, matching the survival_analysis format. FDR is applied **within each test family** (not across families).

| # | Test | Checks | Method | Details |
|---|------|--------|--------|---------|
| 1 | SEX enrichment | "Is there a sex bias in this cancer group compared to 50:50?" | Binomial vs 50:50, per CG | Tests if sex ratio deviates from 50:50. Reports %Male, %Female, which sex is enriched. Plot: horizontal bar chart of %Male/%Female per CG with significance stars. FDR family: `SEX enrichment` |
| 2 | Predisposition profile | "Does the predisposition makeup of this cancer group differ from all other cancers combined?" | Chi-squared per CG (profile vs rest) | Exclude "No predisposition"/"Unknown". For each CG, compare its distribution of predisposition categories to all other CGs combined. Plot: stacked bar of predisposition composition per CG. Also show top predispositions per CG descriptively. FDR family: `Predisposition profile` |
| 3 | SUBTYPE distribution | (descriptive only — no hypothesis test) | Heatmap only | Cross-tab of MOLECULAR_SUBTYPE × CANCER_GROUP with count annotations. No statistical tests. |
| 4 | Race distribution | "Does the racial composition of this cancer group differ from the overall cohort's racial composition?" | Chi-squared per CG (profile vs overall) | For each CG, compare its within-group race distribution to the overall race distribution. Plot: stacked bar of race composition per CG vs overall. FDR family: `Race distribution` |

#### Phase 3 — Numeric Comparisons

| # | Test | Checks | Method | Details |
|---|------|--------|--------|---------|
| 5 | AGE × CANCER_GROUP | "Does AGE differ across cancer groups? Which groups are outliers?" | Global KW + per-group MW (each vs all others) | Kruskal-Wallis across all groups. Then each group vs all others combined via Mann-Whitney to find which groups differ. Boxplot ordered by median. FDR family: `AGE × CG` |
| 6 | TF × CANCER_GROUP | "Does TUMOR_FRACTION differ across cancer groups? Which groups are outliers?" | Global KW + per-group MW (each vs all others) | Same structure as Test 5. FDR family: `TF × CG` |
| 7 | TP × CANCER_GROUP | "Does TUMOR_PLOIDY differ across cancer groups? Which groups are outliers?" | Global KW + per-group MW (each vs all others) | Same structure as Test 5. FDR family: `TP × CG` |
| 8 | AGE × TF correlation | "Is AGE associated with TUMOR_FRACTION?" | Spearman ρ + LOESS scatter | Scatter plot with LOESS smoother. FDR family: `Numeric correlations` |
| 9 | AGE × TP correlation | "Is AGE associated with TUMOR_PLOIDY?" | Spearman ρ + LOESS scatter | Scatter plot with LOESS smoother. FDR family: `Numeric correlations` |
| 10 | TF × TP correlation | "Is TUMOR_FRACTION associated with TUMOR_PLOIDY?" | Spearman ρ + LOESS scatter | Scatter plot with LOESS smoother. FDR family: `Numeric correlations` |

### Workflow

Before building any notebook, the corresponding `.md` instruction file and the relevant section of this plan **must be reviewed and approved**. No notebook is built until the plan reflects the agreed design.

**Notebook convention:** Every test cell must start with a `# ── Checks: "..."` comment stating the hypothesis. See `clinical_deep_dive_general.md §4.5`.

### Phase 4 — Multivariate Models (Stratified Cox PH)

**Goal:** Estimate the *independent* contribution of each clinical variable to OS and EFS, after adjusting for the others.

**Model structure:**
- Stratified Cox PH: `strata(CANCER_GROUP)` — allows different baseline hazards per cancer group without assuming proportional hazards across groups
- One model for OS, one for EFS (identical formula)
- Outcome columns: `OS_MONTHS` / `OS_STATUS` and `EFS_MONTHS` / `EFS_STATUS`

**Covariates:**
| Variable | Type | Notes |
|----------|------|-------|
| AGE | Continuous (years) | Per 1-year increase |
| SEX | Binary (Male/Female) | Female vs Male reference |
| TUMOR_FRACTION | Continuous (0-1) | Per 1-unit increase |
| TUMOR_PLOIDY | Continuous | Per 1-unit increase |
| CANCER_PREDISPOSITIONS | Categorical | Only included if Phase-1 significant (conditional) |

**Outputs:**
| # | Cell Type | Checks / Purpose | Content |
|---|-----------|------------------|---------|
| 1 | 📌 Checks Markdown | "Does each variable independently predict OS after adjusting for the others?" | — |
| 2 | Code | Fit stratified Cox PH (OS) | Print summary: HR, 95% CI, p-value, Wald statistic, concordance index, log-likelihood, AIC |
| 3 | Code | Fit stratified Cox PH (EFS) | Same summary for EFS |
| 4 | Code | Forest plot (OS) | Horizontal: all predictors' HR + 95% CI on log scale. Reference line at HR=1. Color by significance. |
| 5 | Code | Forest plot (EFS) | Same for EFS |
| 6 | Code | Proportional hazards check | Schoenfeld residual test — global + per-variable p-values. Flag violations. |
| 7 | 📌 Checks Markdown | "Is the effect of each predictor consistent across cancer groups?" | — |
| 8 | Code | Subgroup forest plots (OS) | For each predictor: per-group univariate HR + 95% CI vs global stratified HR (reference line). 4 plots. |
| 9 | Code | Subgroup forest plots (EFS) | Same for EFS. 8 plots total across OS and EFS. |
| 10 | Code | Validation | Print sample sizes per group, event counts, verify n≥20 threshold. |

**Statistical details:**
- No FDR correction (single model per outcome — all p-values reported as-is)
- Stratified Cox: `CoxPHFitter(strata=["CANCER_GROUP"])` in lifelines
- PH assumption: `proportional_hazard_test()` (Schoenfeld residuals)
- Effect size: HR with 95% CI (log HR ± 1.96 × SE)
- Concordance: reported as C-index
- PREDISPOSITION included conditionally: check Phase-1 results; if p<0.05 in any univariate test → include

### Phase 5 — Unsupervised Subgroup Discovery

**Goal:** Let the clinical data reveal hidden patient subgroups, then validate them against survival.

**Threshold:** n ≥ 50 per cancer group

**Two parallel feature approaches (compare):**
1. **Numeric-only (PCA + t-SNE):** AGE, TF, TP (standardized)
2. **Mixed via FAMD:** AGE, TF, TP + SEX, CANCER_PREDISPOSITIONS, RACE_GROUP, MOLECULAR_SUBTYPE

**Two stages (both executed):**
- **Stage A — Pooled:** All qualifying patients together → cross-group clusters
- **Stage B — Per-group:** Each cancer group separately → within-type subgroups

#### Stage A — Pooled Discovery

| # | Cell Type | Checks / Purpose | Content |
|---|-----------|------------------|---------|
| 1 | 📌 Checks Markdown | "Can we discover hidden patient subgroups using clinical features?" | — |
| 2 | Code | PCA on numeric (AGE, TF, TP) | Standardize → PCA → scree plot (variance explained per component) + 2D scatter with dropdown to color by CANCER_GROUP / OS_STATUS / MOLECULAR_SUBTYPE |
| 3 | Code | t-SNE on numeric (same raw scaled features) | 2D projection + scatter with same dropdown coloring. Perplexity=30 default. |
| 4 | Code | FAMD on mixed data | AGE, TF, TP (numeric) + SEX, PRED, RACE, SUBTYPE (categorical). 2D scatter with interactive coloring by CG/OS/SUBTYPE. |
| 5 | Code | K-means on PCA embedding | k=2..10: elbow plot (inertia) + silhouette score. Highlight best k. |
| 6 | Code | K-means on FAMD embedding | Same elbow + silhouette. Compare k selections. |
| 7 | Code | Cluster assignment | Assign clusters using best k for PCA. Scatter plots re-colored by cluster label. |
| 8 | Code | Survival validation (PCA clusters) | KM curves for OS + EFS per cluster. Log-rank p-value. |
| 9 | Code | Survival validation (FAMD clusters) | Same for FAMD-based clusters. |
| 10 | Code | Cluster profiles | Table: mean AGE, %Female, top 3 predispositions, race distribution per cluster. |

#### Stage B — Per-Group Discovery

| # | Cell Type | Checks / Purpose | Content |
|---|-----------|------------------|---------|
| 11 | 📌 Checks Markdown | "Within each cancer type, do distinct subpopulations exist?" | — |
| 12 | Code | Loop over CGs ≥ 50 | For each group: standardize AGE, TF, TP → PCA + t-SNE + elbow/silhouette (k=2..10) → pick best k → K-means → scatter with dropdown coloring (by cluster, OS, SUBTYPE) |
| 13 | Code | Survival per group | KM curves (OS + EFS) per cluster within each CG. Log-rank p. |
| 14 | Code | Summary table | One row per CG: N, k chosen, silhouette score, log-rank p-value (OS), log-rank p-value (EFS), top distinguishing features |

#### Plotting details
- All scatter plots (PCA, t-SNE, FAMD) use **Plotly** with **dropdown menus** to switch color-by variable between: CANCER_GROUP / OS_STATUS / MOLECULAR_SUBTYPE / Cluster
- Hover displays: patient ID, all feature values, assigned cluster
- PCA/FAMD: show % variance explained in axis labels (e.g., "PC1 (34.2%)")
- FAMD implemented via: standardize numeric columns, one-hot encode categorical columns, concatenate, then PCA (or use prince.FAMD if available)

#### Statistical details
- No FDR correction (discovery-oriented, not hypothesis-testing)
- Cluster count determined by elbow + silhouette; report both
- Survival validation: log-rank test across clusters within each analysis

---

## References

| File | Location | Purpose |
|------|----------|---------|
| General Methodology | `context/clinical_deep_dive_general.md` | Shared stats, formatting, imports, helpers |
| Phase 1 (Survival Analysis) | `notebooks/survival_analysis/survival_analysis.md` | Outcome Analysis details |
| Phase 2+3 (Associations) | `notebooks/clinical_associations/clinical_associations_analysis.md` | Cross-Categorical + Numeric Comparisons |
| Phase 4+5 (Multivariate + Unsupervised) | `notebooks/clinical_mulltivar_hidden_strcture_analysis/clinical_mulltivar_hidden_strcture_analysis.md` | Multivariate Models + Subgroup Discovery |
| Basic Notebook | `notebooks/clinical_analysis.ipynb` | Exploratory baseline |
| Assignment | `context/assignment.md` | Original task description |
| Reference Paper | `context/referencess.md` | TCGA glioma statistical methods |
