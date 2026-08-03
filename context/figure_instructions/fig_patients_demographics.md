# Figure 1: Patient & clinical demographics (missingness, OS/EFS status, AGE distribution)

## Figure-level info
- Report subsection: "Meta data and clinical analysis"
- LaTeX label / reference: `\ref{fig:patients_demographics}` (referenced as `Figure~\ref{fig:patients_demographics}A`, `... B, C`, `...D`, `...E`)
- Dummy file to replace: `Figures/dummy_1.png`
- Output file: `Figures/fig_patients_demographics.png` (300 dpi) [+ `.pdf`]
- Layout: 3 rows × 2 cols grid (5 panels). Reading order left→right, top→bottom: A, B, C, D, E. Panel A (13-category missingness bar) may span both columns of row 1 if legibility requires; otherwise E occupies row 3 col 1 and (3,2) stays empty.
- Notebooks relied upon:
  - `notebooks/clinical_analysis_excuted.ipynb` — steps 1–4 (per the `% relay on clinical_analysis_excuted step 1-4` comment in the figure environment)
  - src files: `notebooks/clinical_analysis/src/step_01.py`, `step_01b.py`, `step_01c.py`, `step_02.py`, `step_02b.py`, `step_02c.py`, `step_03.py`, `step_03b.py`, `step_03c.py`, `step_03d.py`, `step_04*.py`
- General styling (apply to all panels, per general instruction §4):
  - Panel letter: **bold, uppercase, top-left corner of the axes area** (A, B, C, D, E).
  - Informative axis labels with units; legend whenever >1 series; "Unknown"/missing shown as explicit category where the notebook does.
  - Significance notation per `context/clinical_deep_dive_general.md` §3.6 (`❌`/`✅`/`✅✅`/`✅✅✅`); colors per §4.2 (`px.colors.qualitative.Plotly` / `Set1` / `Set2`).
  - No global figure title needed (caption lives in the .tex), but give each panel an informative title.
  - N = 2870 patients total (executed notebook Step 1). NOTE: the general-instruction worked example (§10) says 2871; the executed notebook prints 2870 — trust the notebook (2870), see Open questions.

## Panel overview
| Panel | Title | Plot type | Notebook step | N | Statistics to annotate |
|-------|-------|-----------|---------------|---|------------------------|
| A | Patient Data — Missingness | Vertical bar (13 columns), value labels on top | clinical_analysis Step 1b (`step_01b.py`) | 2870 | % null per column; highlight EFS_MONTHS 28.4% & OS_MONTHS 27.0% |
| B | Overall Survival Status | Donut/pie (hole 0.3), label + percent | clinical_analysis Step 3 (`step_03.py`) | 2870 | LIVING 1875 (65.3%), DECEASED 640 (22.3%), Unknown 355 (12.4%) |
| C | Event-Free Survival Status | Two grouped bars: EFS binary + EFS detailed | clinical_analysis Step 3b (`step_03b.py`) | 2870 | No Event 1286, Event 1222, Unknown 362; detailed event-type counts |
| D | Age Distribution | Histogram, 40 bins (nbinsx=40) | clinical_analysis Step 2 (`step_02.py`) | 2812 (of 2870) | Mean 9.4, median 8.0, min 0.0, max 73.0, missing 58 (2.0%) |
| E | AGE outliers by cancer group (53/62 → oligodendrogliomas) | Horizontal count bar per cancer group (ad-hoc computation) | No notebook step — ad-hoc; Tukey IQR AGE > 27.125 (verified vs data) | 62 outlier patients (of 2812 with AGE) | 53/62 → Oligodendroglioma (85.5%) |

## Panel A: Patient Data — Missingness
- **Shows:** Share of missing (null) values per patient-level column; demonstrates EFS_MONTHS and OS_MONTHS carry ~28% null values.
- **Data:** `data_clinical_patient_attributes.txt` (13 columns, `read_patients()`), all 2870 patients; NaN counted per column (`df.isna().mean()*100`). No cleaning.
- **Plot type:** Vertical bar chart, bars sorted by % missing descending, value labels (`28.4%`, `27.0%`, ...) above each bar.
- **X axis:** Column name (patient attribute); tick angle −45°, all 13 columns.
- **Y axis:** % Missing (0–100).
- **Grouping/colors:** Single series, crimson bars (`marker_color='crimson'`); highlight OS_MONTHS and EFS_MONTHS bars (distinct color/outline) as the two the text calls out.
- **Statistics to annotate:** % null per column (no tests). None required.
- **N to show:** N = 2870 patients.
- **Exact values that must appear:** EFS_MONTHS **28.4%** (816 null), OS_MONTHS **27.0%** (774 null), RACE 25.9% (744), ETHNICITY 14.6% (420), OS_STATUS 12.4% (355), GERMLINE_SEX_ESTIMATE 11.0% (317), AGE_IN_DAYS 2.0% (58), AGE 2.0% (58), SEX 0.5% (14), CANCER_PREDISPOSITIONS 0.0%, PATIENT_ID 0.0%, EXTERNAL_PATIENT_ID 0.0%, EFS_STATUS 0.0%.
- **Special annotations:** Bold uppercase "A" top-left. Call-out on EFS_MONTHS/OS_MONTHS bars (text claims "~28% Null" for both: EFS 28.4%, OS 27.0%). Note EFS_STATUS is 0% null yet 362 EFS_STATUS values are "1:NA" (shown as Unknown in panel C).

## Panel B: Overall Survival Status
- **Shows:** 65.3% of patients have LIVING status via the OS variable; the rest DECEASED or Unknown.
- **Data:** Patient table; `OS_STATUS` cleaned via `clean_os()` → `os_label` (strip numeric prefix `0:`/`1:`); NaN filled as "Unknown".
- **Plot type:** Donut/pie chart (`hole=0.3`), textinfo `label+percent`.
- **X axis:** n/a (pie/donut).
- **Y axis:** n/a (pie/donut).
- **Grouping/colors:** 3 categories — LIVING, DECEASED, Unknown (e.g. green / red / gray; `px.colors.qualitative.Plotly`).
- **Statistics to annotate:** % per category (in-slice labels). No tests.
- **N to show:** N = 2870.
- **Exact values that must appear:** LIVING **1875 (65.3%)**, DECEASED 640 (22.3%), Unknown 355 (12.4%). Total 1875+640+355 = 2870.
- **Special annotations:** Bold uppercase "B" top-left. Percent labels on slices. (Step 3d validation: OS complete = 2096/2870 pairs; OS_STATUS null = 355.)

## Panel C: Event-Free Survival Status
- **Shows:** Roughly half of patients have no event (1286) and the rest split between different event types; also shows Unknown ("NA") category.
- **Data:** Patient table; `EFS_STATUS` cleaned via `clean_efs()` → `efs_event` (0/1) and `efs_detail`. Note: EFS_STATUS is 100% non-null, but 362 rows carry value `1:NA`, which `clean_efs` maps to NaN → "Unknown" (these are not true missing values).
- **Plot type:** Two grouped bar subplots side-by-side: (1) "EFS Binary" and (2) "EFS Detailed".
- **X axis:** (1) No Event / Event / Unknown; (2) detailed status strings (12 categories).
- **Y axis:** Count (patients).
- **Grouping/colors:** (1) lightgreen / lightcoral / lightgray; (2) single lightcoral series. Value labels above bars.
- **Statistics to annotate:** Counts (no tests). Percentages may be added: No Event 44.8%, Event 42.6%, Unknown 12.6%.
- **N to show:** N = 2870.
- **Exact values that must appear:**
  - Binary: No Event **1286**, Event **1222**, Unknown **362**.
  - Detailed: No Event 1286, Progressive **423**, NA **362**, Recurrence **305**, Progressive - Metastatic 177, Deceased-due to disease 158, Recurrence - Metastatic 88, Second Malignancy 45, Deceased-due to other causes 9, Deceased-causes unavailable 6, Deceased-due to unknown causes 6, Second Malignancy - Metastatic 5.
- **Special annotations:** Bold uppercase "C" top-left. Tilt x labels −45° in detailed subplot (long strings). Keep all 12 detailed categories (do not truncate).

## Panel D: Age Distribution
- **Shows:** AGE at diagnosis is right-skewed with a long right tail (max 73.0 yr in a pediatric dataset) — suggests outliers/misleading points.
- **Data:** Patient table, `AGE` column, drop NaN (`2812/2870` with AGE); AGE stats per Step 2b table.
- **Plot type:** Histogram, `nbinsx=40`, steelblue bars, opacity 0.75, bargap 0.05. (Notebook step also provides a SEX/RACE/ETHNICITY dropdown — the static figure needs only the "Overall" trace.)
- **X axis:** Age (years).
- **Y axis:** Count.
- **Grouping/colors:** Single "All" series (steelblue). No legend needed.
- **Statistics to annotate:** Summary stats (annotate in title or corner box): mean 9.4, median 8.0, min 0.0, max 73.0; missing 58 (2.0%).
- **N to show:** n = 2812 (of 2870; AGE missing 58).
- **Exact values that must appear:** Mean 9.4, Median 8.0, Min 0.0, Max 73.0, Missing 58 (2.0%).
- **Special annotations:** Bold uppercase "D" top-left. Highlight ONLY the Tukey IQR outliers (AGE > 27.125 yr, exactly 62 data points) to visually anchor panel E. Per user decision 2026-08-02: the data give 70 points with AGE > 25, but only the 62 Tukey IQR outliers (AGE > 27.125) are kept/highlighted — matching Panel E exactly.

## Panel E: AGE outliers by cancer group (53/62 → oligodendrogliomas)
- **Shows:** Of the AGE outlier patients (Tukey IQR, AGE > 27.125 yr), most (53/62, 85.5%) are oligodendrogliomas — consistent with an adult cancer appearing in a pediatric dataset.
- **Data:** Patient table AGE + sample table merged on PATIENT_ID for CANCER_GROUP. Outlier = AGE > Q3 + 1.5·IQR where Q1=4.0, Q3=13.25, IQR=9.25 (threshold > 27.125 yr). A patient counts as Oligodendroglioma if any of their samples has CANCER_GROUP == 'Oligodendroglioma'.
- **Plot type:** Horizontal bar chart — count of AGE-outlier patients per cancer group, bars sorted by count descending, only groups with ≥1 outlier shown (8 groups).
- **X axis:** Count of AGE-outlier patients.
- **Y axis:** Cancer group.
- **Grouping/colors:** Single series; highlight the Oligodendroglioma bar (distinct color); others muted/gray.
- **Statistics to annotate:** "53/62 (85.5%) of AGE outliers are Oligodendroglioma" as an annotation.
- **N to show:** N = 62 outlier patients (of 2812 with AGE; 2870 total).
- **Exact values that must appear:** Oligodendroglioma 53, High-grade glioma 3, and 1 each: Neurofibroma/Plexiform, Schwannoma, Low-grade glioma, Diffuse midline glioma, Meningioma, Medulloblastoma. AGE range of outliers 28–73.
- **Special annotations:** Bold uppercase "E" top-left. Definition note: "Tukey IQR outlier (AGE > 27.125 yr), verified to reproduce the 53/62 claim."

## Open questions / flags
- **Panel E — RESOLVED:** Tukey IQR (AGE > 27.125 yr) reproduces the 62-outlier / 53-oligodendroglioma claim exactly (verified computationally against the data on 2026-08-02). Panel E now uses this definition and a count-bar-per-group plot.
- **N discrepancy (2871 vs 2870):** The general-instruction worked example (§10) uses N = 2871, but the executed notebook Step 1 prints "Shape: 2870 rows x 13 columns". All N-dependent values in this file (percentages, counts) are computed from the executed notebook's N = 2870 and must take precedence.
- **Text "~28%" vs exact:** Text says EFS and OS months have "~28% Null". Executed notebook: EFS_MONTHS **28.4%** (816/2870), OS_MONTHS **27.0%** (774/2870). Both shown; do not round OS to 28%.
- **Panel B percentages:** Text only quotes the 65.3% living figure; the DECEASED 22.3% and Unknown 12.4% come from the executed pie (1875/640/355) — included so the donut is complete.
- **Panel C "half/half" phrasing:** Text says "half the patients exhibit no event and the other half split between different events". Executed: No Event 1286 (44.8%), Event 1222 (42.6%), Unknown 362 (12.6%). The Unknown ("NA") slice is not mentioned in the text but must appear (per §4.3 "Unknown" category rule).
- **EFS_STATUS vs EFS "NA":** Missingness bar shows EFS_STATUS 0.0% null, yet 362 rows are `1:NA`. In panel C these appear as Unknown/NA. Not a bug; just record that the binary "Unknown" (362) equals the detailed "NA" (362).
- **Step 4 (predispositions) is part of the relied-upon steps 1–4 but does not feed any panel A–E.** It is listed because the figure comment names steps 1–4; no predispositions content is needed for this figure.
- **Figures 2 and 3 instruction files created (2026-08-02):** `context/figure_instructions/fig_samples_data.md` (Figure 2) and `context/figure_instructions/fig_clinical_association.md` (Figure 3). Figure 2 and Figure 3 references in the report text use explicit "Figure 2/3" numbers (not `\ref`), so the figures are numbered in report order: Figure 1 = patients_demographics, Figure 2 = samples_data, Figure 3 = clinical_association.
