# Phase 2 + 3: Clinical Associations — Cross-Categorical & Numeric Comparisons

**Also see:** `context/clinical_deep_dive_general.md` for shared methodology, imports, and conventions.

## Objective

Analyze associations between clinical variables **without outcome** (Phase 1 already handled outcome). This is a two-part analysis:

- **Phase 2 (Cross-Categorical):** Associations between pairs of categorical clinical variables (e.g., SEX × CANCER_GROUP).
- **Phase 3 (Numeric Comparisons):** Distribution of numeric variables (AGE, TF, TP) across categorical groups + correlations between numeric variables.

## Dataset & Preprocessing

Same merged patient+sample dataset as Phase 1. Files at `PBTA_RNA/data_clinical_patient_attributes.txt` and `data_clinical_sample_attributes.txt`.

### Cleaning
- Reuse same cleaning functions from `notebooks/clinical_analysis/src/imports.py`:
  - `read_patients()`, `read_samples()`
  - `clean_os()`, `clean_efs()`, `clean_race_eth()`, `clean_pred()`, `clean_subtype()`, `clean_tf_tp()`
- **No outcome filtering needed** — these analyses use all available samples (not restricted to OS/EFS available cases)
- Drop rows with missing data per specific test (pairwise deletion)

## Phase 2: Cross-Categorical Associations

### Analyses

| # | Test | Predictor | Response | Min n per cell | Statistic | Effect Size |
|---|------|-----------|----------|---------------|-----------|-------------|
| 1 | SEX × CANCER_GROUP | SEX | CANCER_GROUP | 5 | χ² | Cramer's V |
| 2 | CANCER_PREDISPOSITIONS × CANCER_GROUP | PRED | CANCER_GROUP | 5 | χ² | Cramer's V |
| 3 | MOLECULAR_SUBTYPE × CANCER_GROUP | SUBTYPE | CANCER_GROUP | 5 | χ² | Cramer's V |
| 4 | SEX × CANCER_PREDISPOSITIONS | SEX | PRED | 5 | χ² | Cramer's V |
| 5 | RACE × CANCER_GROUP | RACE | CANCER_GROUP | 5 | χ² | Cramer's V |

### For each test:
1. Build cross-tabulation table
2. Run chi-squared test (or Fisher exact if any cell < 5)
3. Compute Cramer's V as effect size
4. Generate heatmap visualization (annotated with counts or percentages)
5. Add results to summary table

### Edge Cases
- If a category appears in <1% of rows, collapse into "Other"
- If more than 50% of expected cells are < 5, fall back to Fisher exact test
- Report the number of excluded "Other" samples

## Phase 3: Numeric Comparisons

### Analyses

#### Part A: Numeric × Categorical (Kruskal-Wallis + Dunn's post-hoc)

| # | Test | Numeric Var | Groups | Min per group | Statistic | Effect Size |
|---|------|-------------|--------|---------------|-----------|-------------|
| 6 | AGE distribution by CANCER_GROUP | AGE | CANCER_GROUP | 20 | KW | ε² |
| 7 | TUMOR_FRACTION by CANCER_GROUP | TF | CANCER_GROUP | 20 | KW | ε² |
| 8 | TUMOR_PLOIDY by CANCER_GROUP | TP | CANCER_GROUP | 20 | KW | ε² |

For each:
1. Boxplot (grouped by CANCER_GROUP, ordered by median)
2. Kruskal-Wallis test
3. If KW significant (p < 0.05), run Dunn's post-hoc with Bonferroni correction
4. Compute ε² effect size (KW equivalent of η²)
5. Add to summary table

#### Part B: Numeric × Numeric (Spearman Correlation)

| # | Test | Variable 1 | Variable 2 | Min n | Statistic | Effect Size |
|---|------|------------|------------|-------|-----------|-------------|
| 9 | AGE × TF correlation | AGE | TUMOR_FRACTION | 20 | Spearman ρ | ρ |
| 10 | AGE × TP correlation | AGE | TUMOR_PLOIDY | 20 | Spearman ρ | ρ |
| 11 | TF × TP correlation | TF | TUMOR_PLOIDY | 20 | Spearman ρ | ρ |

For each:
1. Scatter plot with LOESS smoother
2. Spearman rank correlation
3. Add to summary table

## Results Saved

- `notebooks/clinical_associations/clinical_associations_results.csv` — full results table
- `notebooks/clinical_associations/clinical_associations_summary.csv` — subset of significant findings

## Cell Outline

### Section 0: Setup
1. **Imports** — load standard + custom helpers
2. **Data loading** — read + merge + clean
3. **Section header**

### Section 1: Phase 2 — Cross-Categorical
4. **Phase 2 header** (markdown)
5. **Test 1:** SEX × CANCER_GROUP
6. **Test 2:** PRED × CANCER_GROUP
7. **Test 3:** SUBTYPE × CANCER_GROUP
8. **Test 4:** SEX × PREDISPOSITION
9. **Test 5:** RACE × CANCER_GROUP
10. **Phase 2 results table**

### Section 2: Phase 3 — Numeric Comparisons
11. **Phase 3 header** (markdown)
12. **Test 6:** AGE × CANCER_GROUP boxplot + KW
13. **Test 7:** TF × CANCER_GROUP boxplot + KW
14. **Test 8:** TP × CANCER_GROUP boxplot + KW
15. **Test 9:** AGE × TF scatter + Spearman
16. **Test 10:** AGE × TP scatter + Spearman
17. **Test 11:** TF × TP scatter + Spearman
18. **Phase 3 results table**

### Section 3: Combined
19. **Combined results table** (Phase 2 + Phase 3)
20. **Save results CSV**

## Outputs & File Locations

| File | Path |
|------|------|
| Source notebook (no outputs) | `notebooks/clinical_associations/clinical_associations_analysis.ipynb` |
| Executed notebook (has outputs) | `clinical_associations_analysis_executed.ipynb` (root) |
| Results CSV | `notebooks/clinical_associations/clinical_associations_results.csv` |
| Summary CSV | `notebooks/clinical_associations/clinical_associations_summary.csv` |
| Build script | `notebooks/clinical_associations/src/build_nb.py` |
| Support .py files | `notebooks/clinical_associations/src/` |

## FDR Note

Phase 2+3 tests involve the same variables as Phase 1 (SEX, AGE, TF, TP, PREDISPOSITION). Currently, FDR is applied within each phase independently. After all phases are complete, a unified cross-phase FDR per variable family should be applied as a sensitivity analysis. See plan document for details.
