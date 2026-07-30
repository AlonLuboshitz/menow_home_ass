# Clinical Deep Dive — Master Plan

This document is the **high-level overview** of the PBTA_RNA clinical deep-dive analysis. It defines the scope, file structure, and cross-references between all components. Detailed methodology and implementation are in the files it references.

---

## File Structure

```
context/
├── clinical_deep_dive_general.md     ← Shared methodology, conventions, imports, helper functions, phases overview
└── clinical_deep_dive_plan.md         ← THIS FILE — high-level reference

notebooks/
├── survival_analysis.md              ← Phase 1: Outcome Analysis (OS & EFS comparisons)
├── clinical_associations_analysis.md ← Phase 2+3: Cross-Categorical + Numeric Associations
├── phase4_deep_dive.md               ← (future) Phase 4: Multivariate Models
├── phase5_deep_dive.md               ← (future) Phase 5: Unsupervised Subgroup Discovery
└── phase6_deep_dive.md               ← (future) Phase 6: Summary & Consolidation
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
| **2 + 3** | Cross-Categorical + Numeric | SEX × CG, PRED × CG, SUBTYPE × CG, AGE × CG, TF × CG, TP × CG + correlations | Association only | χ², Cramer's V, heatmaps, KW + Dunn's, Spearman ρ | n ≥ 20 | ✅ Designed |
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

---

## References

| File | Location | Purpose |
|------|----------|---------|
| General Methodology | `context/clinical_deep_dive_general.md` | Shared stats, formatting, imports, helpers |
| Phase 1 (Survival Analysis) | `notebooks/survival_analysis.md` | Outcome Analysis details |
| Phase 2+3 (Associations) | `notebooks/clinical_associations_analysis.md` | Cross-Categorical + Numeric Comparisons |
| Basic Notebook | `notebooks/clinical_analysis.ipynb` | Exploratory baseline |
| Assignment | `context/assignment.md` | Original task description |
| Reference Paper | `context/referencess.md` | TCGA glioma statistical methods |
