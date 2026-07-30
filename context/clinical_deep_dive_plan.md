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
| **4** | Multivariate Models | All above combined | OS & EFS (time) | Cox PH, HR forest | n ≥ 20 | 🔲 Planned |
| **5** | Unsupervised | AGE, TF, TP, SEX, CG, SUBTYPE, etc. | Cluster discovery | PCA, t-SNE, K-means, FAMD, silhouette, log-rank validation | n ≥ 50 | 🔲 Planned |
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

---

## References

| File | Location | Purpose |
|------|----------|---------|
| General Methodology | `context/clinical_deep_dive_general.md` | Shared stats, formatting, imports, helpers |
| Phase 1 (Survival Analysis) | `notebooks/survival_analysis/survival_analysis.md` | Outcome Analysis details |
| Phase 2+3 (Associations) | `notebooks/clinical_associations/clinical_associations_analysis.md` | Cross-Categorical + Numeric Comparisons |
| Basic Notebook | `notebooks/clinical_analysis.ipynb` | Exploratory baseline |
| Assignment | `context/assignment.md` | Original task description |
| Reference Paper | `context/referencess.md` | TCGA glioma statistical methods |
