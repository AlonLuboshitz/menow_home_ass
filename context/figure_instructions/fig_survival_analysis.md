# Figure 4: Survival analysis (KM per cancer group, variable × outcome, mini-KM, Cox forest)

## Figure-level info
- Report subsection: "Survival analysis" (Results → Survival analysis).
- LaTeX label / reference: literal `Figure 4A`–`Figure 4H`. **FLAG:** the figure environment reuses
  `\label{fig:clinical_association}` (duplicate of Figure 3, report_summary.tex line 143) — fix to
  `\label{fig:survival_analysis}`.
- Dummy file to replace: `Figures/dummy_4.png` (**NOTE:** not present in `Figures/` — only `dummy_1.png`
  and `dummpy_2.png` exist).
- Output file: `Figures/fig_survival_analysis.png` (300 dpi) [+ `.pdf`].
- Layout: 4 rows × 2 cols, panels A–H in reading order (left→right, top→bottom). Row 1: A/B KM by cancer
  group. Row 2: C/D compound panels (2×2 mini-grid). Row 3: E/F compound panels (5 mini-KM). Row 4: G/H
  forest plots.
- Styling: panel letter **bold uppercase** top-left of each axes area (A–H); compound panels use **bold
  lowercase** sub-letters (c1–c4, d1–d4, e1–e5, f1–f5). Significance per `context/clinical_deep_dive_general.md`
  §3.6; FDR (BH) per endpoint for A/B. Colors per §4.2 (CG curves: `px.colors.qualitative.Set1 + Dark2 + Set3`;
  forests: the helper's fixed forestgreen/crimson/lightgray scheme). Informative axis labels with units; legend
  whenever >1 series; "Unknown" as explicit category where the notebook does. No global title (caption lives in
  the .tex); informative title per panel.

## Sources (reuse existing code, combine into one figure)
- **A/B** — `notebooks/clinical_analysis/src/step_12.py` (notebook Step 12) already draws both KM plots:
  legend `{CG} (n={N_OS})` + `*` p<0.05 / `**` FDR<0.05 vs all others, global log-rank annotation, y range
  −0.05–1.05. `step_12b.py` prints the full per-CG p/q table. Split into two subplots, add panel letters A/B.
- **C/D** — `notebooks/survival_analysis/survival_analysis.md` §1A–1D (executed in
  `survival_analysis_executed.ipynb`): box plots (AGE/TF/TP) + stacked bar (SEX) vs OS_STATUS / EFS_STATUS,
  with MW/χ² + effect-size annotations.
- **E/F** — same notebook, §1A–1E time-to-event plots (`.3`/`.4` and 1E.1/1E.2): mini-KM with log-rank p.
- **G/H** — `notebooks/clinical_mulltivar_hidden_strcture_analysis/src/build_nb.py` cells 6/8 (Cox models
  `cph_os`/`cph_efs`) and cells 9/10 (`forest_plot` helper → `fig_os`/`fig_efs`).

## Panel overview
| Panel | Title | Plot type | Source | Key stats |
|-------|-------|-----------|--------|-----------|
| A | OS by Cancer Group | KM step curves, 15 CGs | step_12.py | global log-rank p<0.0001; 15/22 FDR-sig |
| B | EFS by Cancer Group | KM step curves, same 15 CGs | step_12.py | global log-rank p<0.0001; 16/22 FDR-sig |
| C | Variable × OS status | 2×2 mini-grid: box (AGE, TF, TP) + stacked bar (SEX) | survival 1A.1–1D.1 | AGE, TF, SEX sig; TP ns |
| D | Variable × EFS status | 2×2 mini-grid, same layout | survival 1A.2–1D.2 | AGE, TF, TP sig; SEX ns |
| E | Mini KM — OS | 5 mini-KM (AGE, TF, TP, SEX, pred) | survival 1A.3–1E.1 | TF, SEX, pred sig; AGE, TP ns |
| F | Mini KM — EFS | 5 mini-KM | survival 1A.4–1E.2 | AGE, SEX sig; TF, TP, pred ns |
| G | Forest — OS (stratified Cox) | Forest, log HR axis, HR=1 vline | build_nb.py cells 6+9 | only AGE sig; C=0.564, AIC 7777.31 |
| H | Forest — EFS (stratified Cox) | Forest | build_nb.py cells 8+10 | AGE, TP sig; C=0.590, AIC 11145.77 |

## Panel-specific notes
### A — OS by Cancer Group
Reuse step_12.py row 1 as-is. Note "n = complete OS records; * p<0.05, ** FDR<0.05 vs all other groups".
Annotate "global OS log-rank p<0.0001". Groups with n<20 complete OS records are not plotted (22 tested,
15 plotted).

### B — EFS by Cancer Group
Same as A for EFS; same 15 plotted CGs, same labels/colors as A (legend shown once in A). Annotate
"global EFS log-rank p<0.0001".

### C — Variable × OS status (global, binary)
Reuse plots 1A.1/1B.1/1C.1/1D.1; status-colored (green LIVING / red DECEASED). Annotate per mini-panel
`Mann-Whitney U=…, p=…, d=…` or `χ²=…, p=…, V=…`; N = complete cases per test (total cohort 2870 patients).
Verify: AGE U=1,819,966 p=0.0059 d=0.053 (sig); TF U=957,930 p<0.001 d=0.149 (sig); TP U=837,216 p=0.0772
d=−0.033 (ns); SEX χ²=10.35 p=0.0057 V=0.052 (sig).

### D — Variable × EFS status (global, binary)
Same layout vs EFS_STATUS (Event / No Event; "NA"/Unknown excluded or third category per notebook).
Verify: AGE U=2,081,908 p<0.001 d=0.161 (sig); TF U=940,892 p<0.001 d=0.094 (sig); TP U=868,710 p=0.0454
d=−0.036 (sig); SEX χ²=5.81 p=0.0548 V=0.039 (ns).

### E — Mini KM — OS
Reuse 1A.3 (AGE median split: young ≤ median / old > median), 1B.3 (TF), 1C.3 (TP), 1D.3 (SEX), 1E.1 (pred
binary has/no). Log-rank p per subplot; bold the title of significant panels. Verify: AGE p=0.0997 (ns),
TF p=0.0028, TP p=0.0664 (ns), SEX p=0.0032, pred p=0.0019.

### F — Mini KM — EFS
Reuse 1A.4, 1B.4, 1C.4, 1D.4, 1E.2. Verify: AGE p<0.001, TF p=0.5201 (ns), TP p=0.0555 (ns),
SEX p=0.0001, pred p=0.5678 (ns).

### G — Forest — OS
Model (cell 6): `CoxPHFitter(strata=["CANCER_GROUP"])`, formula `AGE + C(SEX) + TUMOR_FRACTION +
TUMOR_PLOIDY`, complete cases (2256 samples, 933 events). Apply to `fig_os` (cell 9):
- HR (95% CI) as **text at the end of each row** (hover-only CI is NOT acceptable in a static PNG — the
  report wants the CI values visible); Wald p per row.
- Add a legend for the three colors: forestgreen = sig & HR<1, crimson = sig & HR>1, lightgray = ns.
- Rename `C(SEX)[T.Male]` → `SEX (Male)` (flag in caption). Title "OS: Stratified Cox PH", subtitle
  "Stratified by CANCER_GROUP"; dashed HR=1 vline labeled "HR=1"; log x-axis (ticks 0.1–10).
Verify: AGE HR=0.966 (0.953–0.979) p<0.001; all others ns. Concordance 0.564, AIC 7777.31.

### H — Forest — EFS
Same as G from `cph_efs` (cell 8) + `fig_efs` (cell 10); 2147 samples, 1384 events. Verify:
AGE HR=0.957 (0.947–0.968) p<0.001; TUMOR_PLOIDY HR=1.096 (1.019–1.180) p=0.014; SEX/TF ns.
Concordance 0.590, AIC 11145.77.

## Open questions / flags
- **Duplicate LaTeX label:** fix to `\label{fig:survival_analysis}` (report_summary.tex line 143); text
  references stay literal "Figure 4A–H".
- **Dummy file missing:** `Figures/dummy_4.png` not present (only `dummy_1.png`, `dummpy_2.png`). Render the
  PNG regardless.
- **N per global test (C–F):** read from executed survival notebook section-1 outputs / Results Summary
  (168 tests, 39 raw sig = 23.2%, 22 FDR sig = 13.1%).
- **Overlapping risk sets (A/B):** 84 patients contribute to ≥2 CG curves (172 entries) — the report should
  acknowledge this (or footnote the figure).
- **TF=1 bias (report_summary.tex lines 147–149):** unresolved note about high-TF samples biasing the TF×OS
  result; keep the raw curves/boxes or add a TF=1-excluded sensitivity before final render.
- **AGE median split (E/F):** young ≤ median vs old > median; the age-deciles notebook is a separate analysis,
  NOT part of Figure 4.
- **Predisposition (E/F):** only the binary has/no comparison is used (OS p=0.0019, EFS p=0.5678); the
  per-type multi-group KM is NOT a panel of Figure 4.
