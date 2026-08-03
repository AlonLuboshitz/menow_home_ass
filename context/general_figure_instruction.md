# General Figure Instruction

## 1. Role & Purpose

This is the master instruction for the **figure-planning agent**. Your job is **NOT** to render figures. Your job is to produce one **figure instruction file** per figure referenced in the analysis report. Each instruction file describes, panel by panel, exactly what the final figure must contain, so that a separate rendering agent can draw it later without re-reading the report.

Pipeline:

```
report_summary.tex  ──┐
notebooks/*_executed.ipynb ──┤→ [you: this task] → context/figure_instructions/fig_*.md → [rendering agent] → Figures/*.png|*.pdf
notebook src files  ──┘
```

## 2. Inputs

| # | Input | Path | What you take from it |
|---|-------|------|----------------------|
| 1 | Analysis report | `/home/alon/menow_home_ass/report_summary.tex` | Which figures exist, which panels (letters) they have, what each panel must show, which notebook it relies on |
| 2 | Dummy figures | `/home/alon/menow_home_ass/Figures/dummy_*.png` | Placeholders to be replaced; the presence of a dummy image marks a figure that must be planned |
| 3 | Executed notebooks | `notebooks/*_executed.ipynb` | The exact plots the analysis produced, with the values/statistics visible in the executed outputs |
| 4 | Notebook source | `notebooks/<name>/src/*.py`, `build_nb.py` | The code behind each step (what each plot computes: variables, tests, annotations) |
| 5 | Methodology conventions | `context/clinical_deep_dive_general.md` | Statistical test selection, significance notation, effect sizes, color conventions |

## 3. Outputs

- **One file per figure**: `context/figure_instructions/fig_<label>.md`
  - If the figure has a LaTeX label (e.g. `fig:patients_demographics`) → `fig_patients_demographics.md`
  - Otherwise → `fig_<N>.md` (N = the figure number as referenced in the report, e.g. `fig_2.md`)
- Every figure referenced in the report gets its own file. Do **not** merge multiple figures into one file.

## 4. General figure conventions (apply to every figure)

These are hard design rules. Every instruction file you write must encode them.

### 4.1 Panels

- Each panel letter (A, B, C, ...) referenced in the report text = one subplot (panel) in the figure.
- Panels are ordered in reading order: left → right, top → bottom.
- Panel letters are **uppercase** (A, B, C, ...), even if the text writes them lowercase ("notebook plot 2b", "notebook plot 3a").
- The panel letter is **bold** and placed in the **top-left corner of the panel**, inside the axes area.
- Panels sharing a row are aligned on the same baseline.

### 4.2 Every panel must have

- **Informative axis labels** — with units when applicable (e.g. "Age at diagnosis (years)", "Months", "Count", "% of patients", "Tumor fraction").
- **Statistics** — any test / p-value / effect size mentioned in the text or computed in the notebook must be annotated on the panel (in the title, an annotation box, or in-plot text). Use the notation from `context/clinical_deep_dive_general.md` §3.6.
- **Sample size N** — annotate N wherever relevant (total N, per-group N, N events).
- **Legend** whenever more than one series/group is shown.

### 4.3 Styling

- Static, publication-ready images for LaTeX `\includegraphics`: **PNG at 300 dpi** (optionally also PDF).
- Consistent font sizes across panels; axis labels legible.
- Colors: `px.colors.qualitative.Plotly` or `Set1` / `Set2`.
- Significance notation: `*` p<0.05, `**` p<0.01, `***` p<0.001 (or ✅/❌ per project convention), plus raw p and FDR q where computed.
- Missing values are shown as an "Unknown" category wherever the notebooks do the same.

### 4.4 Figure level

- The caption lives in the .tex, so panels do not need a global figure title, but each panel should carry an informative title when it aids readability.
- Figure width ~ `0.7\textwidth`; the panel grid (rows × cols) is chosen by you based on the number of panels.

## 5. Workflow

Follow these steps over the whole report:

### Step 1 — Read the report
Read `/home/alon/menow_home_ass/report_summary.tex` end to end. Focus on the `Results` section and its subsections.

### Step 2 — Find figures and panels
For each subsection:

1. Locate every `\begin{figure}...\end{figure}` environment that includes a dummy image (`Figures/dummy_*.png`). Each one is a figure to plan.
2. Locate every text reference to a figure plus panel letters:
   - `Figure ~\ref{fig:patients_demographics}A` → Figure 1, panel A
   - `Figure ~\ref{fig:patients_demographics} B, C` → Figure 1, panels B and C (two separate panels)
   - `Figure 2A` → Figure 2, panel A
   - `(notebook plot 2b)` / `(notebook plot 3a)` → panels B / A of Figures 2 / 3
3. Build the panel list for the figure: `A, B, C, ...` in the order the letters first appear in the text.

### Step 3 — Identify the relied-upon notebook(s)
Priority order:

1. **LaTeX `%` comment inside the figure environment** (e.g. `% relay on clinical_analysis_excuted step 1-4`). This is authoritative.
2. Text references of the form "notebook plot Xn".
3. Content-keyword → notebook mapping (fallback):
   - missingness / demographics / survival status / predispositions → `clinical_analysis_excuted.ipynb`
   - sex / predisposition / race enrichment, AGE/TF/TP per group, correlations → `clinical_associations_analysis_executed.ipynb`
   - AGE / TF / TP / SEX / predisposition × outcome, KM + log-rank → `survival_analysis_executed.ipynb`
   - age deciles → `age_deciles_analysis_executed.ipynb`
   - clustering / hidden structure → `clinical_mulltivar_hidden_strcture_analysis_executed.ipynb`

### Step 4 — Read the notebooks
For each relied-upon notebook:

- Open `notebooks/<name>_executed.ipynb` and read the markdown headings (`## Step N: ...`) to locate the steps that match the panels.
- Read the corresponding step source files under `notebooks/<name>/src/` (or `build_nb.py`) to see exactly what each plot computes.
- Read the executed outputs in the notebook for exact numbers (percentages, counts, p-values, N).

### Step 5 — Map each panel to a notebook step
For every panel letter decide:

- which notebook step / cell produces the plot it needs,
- which data / columns it uses and any cleaning,
- the plot type (bar, histogram, pie/donut, box, KM curve, heatmap, scatter, ...),
- which statistics are computed,
- which N is available.

### Step 6 — Extract exact values
Record the exact values from the executed outputs so the rendering agent does not need to re-run analysis:

- percentages (e.g. ~28% null, 65.3% living),
- counts (e.g. 53/62),
- test statistics, p-values, FDR q-values,
- N per group / total N.

### Step 7 — Write the instruction files
One file per figure, following the template in §6.

### Step 8 — Validate and flag
Run the checklist in §9. Anything unresolved → put it under "Open questions" in the instruction file and report it to the user. Never silently guess.

## 6. Per-figure instruction file template

```markdown
# Figure N: <short title>

## Figure-level info
- Report subsection: <subsection title in report_summary.tex>
- LaTeX label / reference: <\ref{...} or "Figure N">
- Dummy file to replace: Figures/dummy_X.png
- Output file: Figures/<name>.png (300 dpi) [+ .pdf]
- Layout: <rows × cols grid>
- Notebooks relied upon:
  - notebooks/<name>_executed.ipynb — steps <N-M>
  - src files: <paths>

## Panel overview
| Panel | Title | Plot type | Notebook step | N | Statistics to annotate |
|-------|-------|-----------|---------------|---|------------------------|
| A | ... | ... | ... | ... | ... |
| B | ... | ... | ... | ... | ... |

## Panel A: <title>
- **Shows:** <one sentence — what finding this panel conveys>
- **Data:** <table + columns + cleaning used>
- **Plot type:** <bar / histogram / box / KM / heatmap / scatter ...>
- **X axis:** <label (units)>
- **Y axis:** <label (units)>
- **Grouping/colors:** <series, legend, color mapping>
- **Statistics to annotate:** <test, statistic, p/q value, effect size>
- **N to show:** <sample size / per-group N>
- **Exact values that must appear:** <the numbers from the text/notebook>
- **Special annotations:** <threshold lines, significance brackets, "Unknown" categories, sorting instructions>

## Panel B: <title>
...

## Open questions / flags
- <anything unresolved — missing notebook step, text–notebook mismatch, etc.>
```

## 7. Panel → notebook mapping rules (reference examples)

| Panel topic (from text) | Notebook / step | Notes |
|-------------------------|-----------------|-------|
| Patient missingness incl. EFS/OS months ~28% null | clinical_analysis Step 1 (missingness bar) | highlight OS_MONTHS / EFS_MONTHS |
| OS status (65.3% living) | clinical_analysis Step 3 | LIVING / DECEASED / Unknown |
| EFS events (binary + detailed) | clinical_analysis Step 3b | No Event vs event types |
| AGE histogram (right-skewed) | clinical_analysis Step 2 | ~40 bins |
| Age outliers by cancer group (53/62 → oligodendrogliomas) | no dedicated step — ad-hoc | flag it (§8) |
| Sample missingness (TF/TP/RNA selection/matching ids 20–40%) | clinical_analysis Step 5 | sample-level missingness |
| Cancer group counts (LGG 862, ...) | clinical_analysis Step 6 | bar chart |
| Molecular subtype × cancer group alignment | clinical_analysis Step 8b | heatmap, row % |
| SEX enrichment vs 50:50 | clinical_associations Test 1 | binomial; horizontal bar; color by significance |
| Predisposition composition per CG | clinical_associations Test 2 | stacked bar / heatmap |
| Race composition vs overall cohort | clinical_associations Test 5 | Chi² GoF; stacked bar |
| TF vs OS_STATUS, KM (tf=1 bias) | survival_analysis 1B | box + KM |

The mapping table above is a reference, not a limit. The rule is: every panel in every instruction file must name a concrete notebook step.

## 8. Ambiguity and conflict handling

When you hit any of the following, do **NOT** silently guess:

- A panel is referenced in the text but no notebook step produces it → note it under "Open questions" and propose the closest existing plot, or describe what new computation is needed (variables, filter, test).
- The text numbers differ from the executed notebook → trust the executed notebook, record both values, flag the mismatch.
- A `%` comment names steps that do not contain the referenced plot → map by content and flag it.
- Figure referenced with lowercase letters ("notebook plot 3a") → still plan it as an uppercase panel A.
- Panels referenced together ("B, C") → they are separate panels, not one combined panel.

Report every flag to the user in your final summary.

## 9. Quality checklist

- [ ] One instruction file per figure referenced in the report.
- [ ] Every panel letter that appears in the text appears in the file (and vice versa).
- [ ] Panels ordered A, B, C, ... in reading order.
- [ ] Panel letter = bold, uppercase, top-left corner of the panel (encoded in the styling notes of each panel spec).
- [ ] Every panel spec has: plot type, data source, X/Y axis labels with units, statistics (or explicit "none"), N.
- [ ] Exact numbers from the text / executed notebook are embedded in each panel spec.
- [ ] A concrete notebook step is named for every panel.
- [ ] Layout (rows × cols) is specified.
- [ ] Output file name + dpi are specified.
- [ ] Ambiguities are flagged, not guessed.

## 10. Worked example — Figure 1 (fig:patients_demographics)

Source subsection: "Meta data and clinical analysis". The figure environment comment says: `% relay on clinical_analysis_excuted step 1-4`.

Text extraction:
- **A**: "event-free survival (EFS) and overall survival (OS) months variables has ~28% Null values"
- **B**: "65.3% has a living status via the OS variable"
- **C**: "half the patients exhibit no event and the other half splitted between different events" (EFS)
- **D**: "AGE ... right skewed histogram suggests outliers"
- **E**: "most of this data points (53/62) belong to oligodendrogliomas patients"

Mapping:

| Panel | Content | Notebook step | N | Statistics |
|-------|---------|---------------|----|------------|
| A | Missingness bar chart of patient columns, highlighting OS_MONTHS and EFS_MONTHS (~28.4% null per executed Step 1 output) | clinical_analysis Step 1 | 2871 | % null per column |
| B | OS status donut/pie: LIVING 65.3%, DECEASED, Unknown | Step 3 | 2871 | % per category |
| C | EFS status bars: binary (Event / No Event / Unknown) plus detailed event types | Step 3b | 2871 | counts |
| D | AGE histogram, ~40 bins | Step 2 | patients with AGE | — |
| E | AGE outliers (e.g. AGE above a threshold) annotated by CANCER_GROUP → oligodendrogliomas 53/62 | NO matching step in 1–4 → flag (§8); propose ad-hoc: filter AGE, crosstab with CANCER_GROUP | 62 outliers | 53/62 |

The resulting instruction file for Figure 1 would fill in the §6 template with these five panels, list `clinical_analysis_excuted.ipynb` steps 1–4 as the relied-upon notebook, propose a 3×2 grid layout, and put panel E under "Open questions" until the user confirms the intended outlier definition.
