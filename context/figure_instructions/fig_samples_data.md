# Figure 2: Sample-level data (missingness, cancer-group sizes, molecular subtype × cancer-group alignment)

## Figure-level info
- Report subsection: "Meta data and clinical analysis"
- LaTeX label / reference: `\ref{fig:samples_data}` (referenced as `Figure~\ref{fig:samples_data}A`, `...B`, `...C`; resolves to "Figure 2" in report order)
- Dummy file to replace: `Figures/dummy_2.png`
- Output file: `Figures/fig_samples_data.png` (300 dpi) [+ `.pdf`]
- Layout: 3 rows × 1 col grid (3 panels). Reading order top→bottom: A, B, C. Each panel carries many categorical labels, so the single-column stacking keeps x-axis labels legible.
- Notebooks relied upon:
  - `notebooks/clinical_analysis_excuted.ipynb` — steps 5–9 (per the `% relay on clinical_analysis_excuted step 5-9` comment in the figure environment)
  - src files: `notebooks/clinical_analysis/src/step_05.py`, `step_05b.py`, `step_06.py`, `step_06c.py`, `step_08.py`, `step_08b.py`, `step_08c.py`, `step_09.py`
- General styling (apply to all panels, per general instruction §4):
  - Panel letter: **bold, uppercase, top-left corner of the axes area** (A, B, C).
  - Informative axis labels with units; legend whenever >1 series; "Unknown"/missing shown as explicit category where the notebook does.
  - Significance notation per `context/clinical_deep_dive_general.md` §3.6; colors per §4.2 (`px.colors.qualitative.Plotly` / `Set1` / `Set2`).
  - No global figure title needed (caption lives in the .tex), but give each panel an informative title.
  - N = 4312 samples total (executed notebook Step 5 prints "Shape: 4312 rows x 24 columns"). 4074 samples carry a CANCER_GROUP (238 null, 5.5%). All counts in this file are sample-level.

## Panel overview
| Panel | Title | Plot type | Notebook step | N | Statistics to annotate |
|-------|-------|-----------|---------------|---|------------------------|
| A | Sample Data — Missingness | Vertical bar (24 columns), value labels on top | clinical_analysis Step 5/5b (`step_05.py`, `step_05b.py`) | 4312 | % null per column; highlight TF 33.3%, TP 31.9%, RNA selection 42.5%, matched-normal ids 22.0% |
| B | Cancer Group — Case Counts (Top 12) | Vertical bar, value labels on top | clinical_analysis Step 6 (`step_06.py`) | 4312 | Top-12 counts (LGG 862 … Oligodendroglioma 70); 55 CGs total |
| C | Molecular Subtype × Cancer Group (13 multi-CG subtypes) | Heatmap, row % (cell text = counts) | clinical_analysis Step 8b (`step_08b.py`) | 4312 | 13 subtypes in ≥2 distinct CGs; per-cell counts |

## Panel A: Sample Data — Missingness
- **Shows:** Sample-level columns carry substantial null values; the columns the text calls out (Tumor fraction, ploidy, RNA selection, matching-id samples) sit in the ~20–40% missing band.
- **Data:** `data_clinical_sample_attributes.txt` (24 columns, `read_samples()`), all 4312 samples; NaN counted per column (`df.isna().mean()*100`). No cleaning.
- **Plot type:** Vertical bar chart, bars sorted by % missing descending, value labels above each bar (mirrors `step_05b`; Reds color scale by value, or a single crimson series).
- **X axis:** Column name (sample attribute); tick angle −45°, all 24 columns.
- **Y axis:** % Missing (0–100).
- **Grouping/colors:** Single series; highlight the 5 text-called-out columns (TUMOR_FRACTION, TUMOR_PLOIDY, RNA_LIBRARY_SELECTION, MATCHED_NORMAL_SAMPLE_ID, MATCHED_NORMAL_SPECIMEN_ID) with a distinct color/outline.
- **Statistics to annotate:** % null per column (no tests). None required.
- **N to show:** N = 4312 samples.
- **Exact values that must appear:** RNA_LIBRARY_SELECTION **42.5%** (1834), TUMOR_FRACTION **33.3%** (1435), TUMOR_PLOIDY **31.9%** (1374), SAMPLE_TYPE 23.7% (1023), MATCHED_NORMAL_SPECIMEN_ID **22.0%** (949), MATCHED_NORMAL_SAMPLE_ID **22.0%** (949), MOLECULAR_SUBTYPE 20.5% (882), ONCOTREE_CODE 19.9% (860), EXTENT_OF_TUMOR_RESECTION 7.2% (310), CANCER_TYPE 6.9% (297), CANCER_GROUP 5.5% (238), CNS_REGION 3.1% (133), CANCER_TYPE_DETAILED 2.9% (125), PATHOLOGY_FREE_TEXT_DIAGNOSIS 2.8% (122), BROAD_HISTOLOGY 1.1% (47), and 0.0% for the remaining 9 columns (PATIENT_ID, SAMPLE_ID, CBTN_TUMOR_TYPE, COLLECTION_EVENT_ID, EXPERIMENT_STRATEGY, SPECIMEN_ID, SUB_COHORT, TUMOR_TISSUE_TYPE, TUMOR_TYPE).
- **Special annotations:** Bold uppercase "A" top-left. Call-out on the 5 columns the text names ("~20-40% Null"): TF 33.3%, TP 31.9%, RNA selection 42.5%, matched-normal ids 22.0%. "Matching id samples" = the MATCHED_NORMAL_SAMPLE_ID / MATCHED_NORMAL_SPECIMEN_ID pair. Note RNA_LIBRARY_SELECTION (42.5%) sits just above the text's "40%" upper bound — show as-is (trust the notebook). All 9 zero-missing columns are sample identity/metadata columns (PATIENT_ID, SAMPLE_ID, SPECIMEN_ID, COLLECTION_EVENT_ID, CBTN_TUMOR_TYPE, SUB_COHORT, TUMOR_TYPE, TUMOR_TISSUE_TYPE, EXPERIMENT_STRATEGY) — missingness is concentrated in the assay/annotation columns, matching the text's emphasis.

## Panel B: Cancer Group — Case Counts (Top 12)
- **Shows:** Low-grade glioma (LGG) is the largest cancer group (862 samples); the cohort spans 55 cancer groups, with a long tail of small groups.
- **Data:** Sample table, `CANCER_GROUP` (non-null, 4074 samples); `value_counts().head(12)` exactly as Step 6 computes it (`step_06.py`). No cleaning.
- **Plot type:** Vertical bar chart, lightgreen bars, value labels outside above each bar; bars in descending count order (`value_counts()` order). Mirrors the Cancer Group subplot of Step 6.
- **X axis:** Cancer group; tick angle −45°.
- **Y axis:** Count (samples).
- **Grouping/colors:** Single series (lightgreen). No legend needed.
- **Statistics to annotate:** Counts per group (value labels). Annotate "55 cancer groups total" in the title or as a note.
- **N to show:** N = 4312 samples; top-12 CGs shown (of 55).
- **Exact values that must appear:** Low-grade glioma **862**, High-grade glioma **512**, Medulloblastoma **440**, Diffuse midline glioma **421**, Ependymoma **341**, Atypical Teratoid Rhabdoid Tumor **152**, Ganglioglioma **151**, Adamantinomatous Craniopharyngioma **115**, Meningioma **100**, Choroid plexus tumor **100**, Dysembryoplastic neuroepithelial tumor **73**, Oligodendroglioma **70**.
- **Special annotations:** Bold uppercase "B" top-left. Per user decision 2026-08-02: show the top-12 vertical bar (NOT a binned histogram). The true bin distribution (for the flag below): 2 groups ≥500, 3 groups 300–500, 7 groups 70–200, 43 groups <70. Note Step 6 actually renders a 2×2 grid (Broad Histology / Cancer Group / CNS Region / Tumor Type); this panel extracts only the Cancer Group subplot. The Step 6c printout "Most common CG: Low-grade glioma (862)" / "Total CGs: 55" supplies the panel annotation; the report's statement "LGG has most cases (862)" is reproduced exactly by this bar.

## Panel C: Molecular Subtype × Cancer Group (13 multi-CG subtypes)
- **Shows:** 13 molecular subtypes belong to 2 or more distinct cancer groups — i.e., most subtypes align to a single CG, but these 13 are shared (suggesting mislabeled or genuinely multi-CG subtypes).
- **Data:** Sample table; `MOLECULAR_SUBTYPE` cleaned via `clean_subtype()` (blank + `To be classified` → `Unclassified`); crosstab subtype × CANCER_GROUP over samples with a non-null CG (4074). Rows = the 13 subtypes present in ≥2 distinct CGs; columns = the cancer groups they touch.
- **Plot type:** Heatmap (rows = subtype, cols = cancer group), cell text = sample counts, color = row proportion (`row %`, `cf.div(cf.sum(axis=1))`) — mirrors `step_08b.py` (Blues color scale, y-axis reversed, x tick angle −45°).
- **X axis:** Cancer group; tick angle −45°.
- **Y axis:** Molecular subtype.
- **Grouping/colors:** One row per subtype; colorscale encodes row proportion (counts as text). No legend needed.
- **Statistics to annotate:** None (descriptive, no tests). Annotate "13 subtypes in ≥2 distinct CGs" (matches the report text).
- **N to show:** N = 4074 samples with a CG (of 4312); 13 multi-CG subtypes.
- **Exact values that must appear (subtype — total samples (samples w/ CG) — distinct CGs):**
  - Unclassified — 882 (644) — 26 distinct CGs
  - LGG, To be classified — 81 (81) — 3 (Low-grade glioma, Pilocytic astrocytoma, Pleomorphic xanthoastrocytoma)
  - HGG, To be classified — 33 (33) — 4 (High-grade glioma, Astrocytoma, Diffuse intrinsic pontine glioma, Astroblastoma)
  - GNT, wildtype — 17 (17) — 2 (Glial-neuronal tumor NOS, Desmoplastic infantile astrocytoma and ganglioglioma)
  - NBL, To be classified — 17 (17) — 2 (Neuroblastoma, Ganglioneuroblastoma)
  - GNT, FGFR — 14 (14) — 2 (Glial-neuronal tumor NOS, Rosette-forming glioneuronal tumor)
  - GNT, BRAF V600E — 13 (13) — 2 (Glial-neuronal tumor NOS, Desmoplastic infantile astrocytoma and ganglioglioma)
  - GNT, RTK — 9 (9) — 2 (Glial-neuronal tumor NOS, Desmoplastic infantile astrocytoma and ganglioglioma)
  - SEGA, wildtype — 9 (9) — 2 (Low-grade glioma, Subependymal Giant Cell Astrocytoma)
  - NBL, MYCN amplified — 7 (7) — 2 (Neuroblastoma, Ganglioneuroblastoma)
  - GNT, KIAA1549-BRAF — 6 (6) — 2 (Glial-neuronal tumor NOS, Diffuse leptomeningeal glioneuronal tumor)
  - GNT, To be classified — 5 (5) — 2 (Glial-neuronal tumor NOS, Diffuse leptomeningeal glioneuronal tumor)
  - NBL, MYCN non-amplified — 5 (5) — 2 (Neuroblastoma, Ganglioneuroblastoma)
  - (Sum 1098 total / 860 with CG.)
- **Special annotations:** Bold uppercase "C" top-left. Keep the `Unclassified` row (merged placeholder for blank/`To be classified` — 882 total, of which 238 have a missing CG). Show all 13 rows even though the notebook's step-8b heatmap restricts to top-15 subtypes (only Unclassified and LGG, To be classified appear there).

## Open questions / flags
- **Panel B bin claims — RESOLVED (2026-08-02):** The report text claims "3 groups 300-500, 7 groups 70-200, 45 groups <70". Verified against the executed data the real distribution is: 2 groups ≥500 (LGG 862, HGG 512), 3 groups 300–500 (Medulloblastoma 440, Diffuse midline glioma 421, Ependymoma 341), 7 groups 70–200 (ATRT 152, Ganglioglioma 151, Adamantinomatous Craniopharyngioma 115, Meningioma 100, Choroid plexus tumor 100, DNET 73, Oligodendroglioma 70), 43 groups <70. User chose the top-12 bar over the binned histogram; the bins are recorded here for reference only.
- **Panel C "13 subtypes" list — RESOLVED (2026-08-02):** The 13-subtype list was computed from `data_clinical_sample_attributes.txt` with `clean_subtype()` (the exact data/cleaning step 8/8b uses); the notebook does not print this list as a table, but the computed count (13 subtypes in ≥2 distinct CGs) matches the report text exactly. All counts in Panel C were verified this way.
- **Panel C step attribution — RESOLVED (2026-08-02):** The task brief said "step 9 / 9a" for Panel C, but the subtype × CG heatmap lives in step 8b (step 9 = sequencing strategy/library; step 9a = multi-CG co-occurrence). Panel C mirrors step 8b's row-% heatmap style. The step-8b selection (top-15 subtypes × top-10 CGs) is extended to the full 13-subtype set per user decision.
- **Panel A "~20-40%" wording — RESOLVED (2026-08-02):** Text says "~20-40% Null" for TF/TP/RNA selection/matching ids. Executed values: TF 33.3%, TP 31.9%, RNA_LIBRARY_SELECTION 42.5% (just above the 40% bound), matched-normal ids 22.0%. Show as-is; highlight these 5 columns.
- **Dummy file name — flag (no decision needed):** The report references `Figures/dummy_2.png`; the file present in `Figures/` is misspelled `dummpy_2.png`. Rendering agent should write `Figures/fig_samples_data.png` regardless.
