# General Notebook Instruction

## 1. Role & Purpose

This is the master instruction for the **notebook-planning agent**. Your job is **NOT** to render a single multi-panel figure. Your job is to create, for each figure referenced in the analysis report (identified by a dummy figure and/or an existing figure instruction file), a **Jupyter notebook** in which **every panel is produced by its own code cell**, with the cell code hidden so that only the panel output is visible.

Each cell output must show:
1. the panel label (e.g. `1.A`, `1.B`, ...) — figure number + "." + uppercase panel letter, and
2. the rendered panel itself (static matplotlib image).

This notebook path is an **alternative** to the figure-rendering path of `context/general_figure_instruction.md`. The per-figure instruction files in `context/figure_instructions/fig_*.md` remain the **source of truth for panel content** — you translate each panel spec into a notebook cell. Do not modify `general_figure_instruction.md` or any `fig_*.md` file (we may return to the figure path later).

Pipeline:

```
report_summary.tex ──┐
context/figure_instructions/fig_*.md ──┤→ [you: this task] → notebooks/fig_<label>.ipynb (executed, code hidden, panels as cells)
notebooks/*_executed.ipynb ──┘
```

## 2. Inputs

| # | Input | Path | What you take from it |
|---|-------|------|----------------------|
| 1 | Analysis report | `/home/alon/menow_home_ass/report_summary.tex` | Which figures exist, which subsection they belong to, the panel letters referenced |
| 2 | Figure instruction files | `context/figure_instructions/fig_*.md` | **Main source** — per-panel specs: plot type, data, axes, statistics, N, exact values, styling |
| 3 | Executed notebooks | `notebooks/*_executed.ipynb` | Reference computations and exact numbers (verify panel specs) |
| 4 | Notebook source | `notebooks/<name>/src/*.py`, `build_nb.py` | Underlying code per step (reuse logic in cells) |
| 5 | Methodology conventions | `context/clinical_deep_dive_general.md` | Cleaning helpers (§7), significance notation (§3.6), styling (§4.2) |

## 3. Outputs

- **One executed notebook per figure**: `notebooks/fig_<label>.ipynb`, where `<label>` matches the figure instruction file name (e.g. `fig_patients_demographics` → `notebooks/fig_patients_demographics.ipynb`).
- The notebook contains **one code cell per panel**, in reading order A, B, C, ... — **the first cell renders Panel A**.
- The notebook must be **executed** (outputs embedded, including the labels and the rendered panel images).
- Cell code is **hidden** (see §4.4).
- The notebook is a **production deliverable**: it contains **only** the panel cells (labels + panels). No "Open questions", flag, or note cells are ever added. Any ambiguity is raised to the user and resolved *before* delivery (see §8).
- The notebook begins with a **big title markdown cell**: `# Figure N` (N = figure number as used in the report text, e.g. `# Figure 1`), optionally with a short descriptive subtitle on the next line. The first **code** cell renders Panel A; no other cells precede the title.
- Optional: you may keep a build script (e.g. `notebooks/figure_notebooks/src/build_<label>.py`) using the project's build-nb pattern, but the required deliverable is the executed `.ipynb`.

## 4. General notebook conventions (apply to every panel cell)

### 4.1 Cell per panel

- One code cell per panel letter, in reading order A, B, C, ...
- The notebook's **first cell** is the big title markdown cell (`# Figure N`, see §3). The **first code cell** renders Panel A. No other setup cells precede it.
- Every cell is **self-contained**: it repeats imports, data loading, cleaning, computation and plotting (same convention as the project's step-based notebooks), so any cell can be run independently.

### 4.2 Output label + panel

- The cell's output begins with the label line: `<figure_number>.<UPPERCASE_LETTER>` — e.g. `print("1.A")` for Figure 1 panel A, `print("2.B")` for Figure 2 panel B. The figure number is the number used in the report text ("Figure N").
- The rendered panel follows the label in the same cell output (matplotlib inline image).

### 4.3 Panel letter on the plot

- The panel letter is **also drawn on the plot itself**: **bold, uppercase, top-left corner, ABOVE the axes** (outside the plot area), e.g.:
  ```python
  ax.text(0.0, 1.02, "A", transform=ax.transAxes, fontsize=15, fontweight="bold", va="bottom", ha="left")
  ```
- Both the printed label (`1.A`) and the in-plot letter (`A`) appear.

### 4.4 Code hidden

- Every code cell must carry metadata that hides its input when rendered:
  ```json
  "metadata": {"jupyter": {"source_hidden": true}, "tags": ["hide-input"]}
  ```
- Verify hidden inputs by converting the executed notebook to HTML:
  ```bash
  jupyter nbconvert --to html --TagRemovePreprocessor.remove_input_tags=hide-input notebooks/fig_<label>.ipynb
  ```
  The HTML must show no `In [n]:` code prompts — only labels and panels.

### 4.5 Static inline outputs

- Use `%matplotlib inline` so each panel is embedded as a static PNG output (no `savefig` needed for the notebook; the image is stored in the notebook).
- Styling of each panel follows the figure instruction file: informative axis labels with units, statistics annotations, N, legends, "Unknown" categories where applicable, colors from `px.colors.qualitative.Plotly` / `Set1` / `Set2`.

### 4.6 Exact values

- Use the exact values recorded in the `fig_*.md` panel specs (they were verified against the executed notebooks). Compute from data where possible and assert/print the values so mismatches are visible.

## 5. Workflow

### Step 1 — Read the report
Read `/home/alon/menow_home_ass/report_summary.tex`. Identify the figure(s) and the subsection(s) they belong to, and the panel letters referenced.

### Step 2 — Read the figure instruction file(s)
For each figure to build, read `context/figure_instructions/fig_*.md`. This is the authoritative panel-by-panel spec. If no instruction file exists for a figure, **flag it (§8)** — do not invent panel content.

### Step 3 — Read the relied-upon notebooks
Read the executed notebook(s) and src files named in the figure instruction file to confirm the computations and numbers behind each panel.

### Step 4 — Translate each panel into a cell
For every panel letter in the figure instruction file's overview table, write one self-contained code cell implementing that panel's spec (data → computation → label → plot).

### Step 5 — Build the notebook
Create `notebooks/fig_<label>.ipynb` (directly with `nbformat`, or via a small build script). Prepend the big title markdown cell (`# Figure N`, §3), then one code cell per panel; add the hide-input metadata to every code cell.

### Step 6 — Execute the notebook
Run the notebook end-to-end (e.g. `jupyter nbconvert --to notebook --execute --inplace notebooks/fig_<label>.ipynb`) so labels and panel images are embedded in the outputs.

### Step 7 — Verify
- Convert to HTML with `hide-input` tags removed; confirm no code prompts appear.
- Confirm each cell output starts with the correct label (`1.A`, `1.B`, ...).
- Confirm the values in the outputs match the figure instruction file (or print verification asserts).
- Confirm cell order = panel order A, B, C, ...
- Confirm the notebook contains **no** "Open questions" / flag / note cells (it is production-ready).

### Step 8 — Validate, question the user, resolve
The delivered notebook must be production-ready: **nothing unresolved stays in it**. If anything is ambiguous or unresolved → **do not guess and do not add it to the notebook**. Raise it to the user as a question/flag (ask during the build or in your final report), apply the user's decision to the notebook, re-execute, and re-verify. Only deliver once every flag is resolved and no "Open questions" content remains.

## 6. Cell template

```python
# Panel 1.A — <title>   (from context/figure_instructions/fig_patients_demographics.md)
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- data: self-contained (reuse cleaning helpers from clinical_deep_dive_general.md §7) ---
PATIENT_FILE = "/home/alon/menow_home_ass/PBTA_RNA/data_clinical_patient_attributes.txt"
SAMPLE_FILE  = "/home/alon/menow_home_ass/PBTA_RNA/data_clinical_sample_attributes.txt"
pat = pd.read_csv(PATIENT_FILE, sep="\t", header=4)
...

# --- computation: exactly as the fig_*.md panel spec describes ---
...

# --- output: label first, then the panel ---
print("1.A")                       # <figure_number>.<UPPERCASE_LETTER>
fig, ax = plt.subplots(...)        # panel per fig_*.md spec
...
ax.text(0.0, 1.02, "A", transform=ax.transAxes, fontsize=15,
        fontweight="bold", va="bottom", ha="left")   # bold uppercase above axes
plt.show()
```

Every code cell gets:
```json
"metadata": {"jupyter": {"source_hidden": true}, "tags": ["hide-input"]}
```

The notebook's first cell is a title markdown cell, e.g.:
```markdown
# Figure 1
```
(no hide-input metadata needed on the title cell).

## 7. Panel → cell rules

- Map each row of the figure instruction file's "Panel overview" table to exactly one cell.
- Follow the panel spec fields literally: plot type, data source, X/Y axis labels with units, grouping/colors, statistics to annotate, N to show, exact values, special annotations.
- In the notebook, each panel is its **own standalone plot** (ignore the multi-panel "Layout" field of the figure instruction file — that field applies only to the single-figure path).
- Keep the figure instruction file's panel-letter styling (bold uppercase, above the axes).
- If a panel spec contains unresolved "Open questions", **do not** mirror them into the notebook. Surface them to the user, apply the resolution, and re-verify the affected panels.

## 8. Ambiguity and conflict handling

- A figure has no `fig_*.md` instruction file → flag it; do not build a notebook for it yet.
- A panel is referenced in the text but missing from the figure instruction file → flag it.
- Values in the figure instruction file differ from the executed notebook → follow the figure instruction file (it was already reconciled) and flag the mismatch.
- Unresolved "Open questions" inside the figure instruction file → do not guess and do not add them to the notebook; raise them to the user and apply the decision before delivery.
- Panels referenced together ("B, C") remain separate cells.

Report every flag/question to the user in your final summary and confirm each was **resolved and applied** to the notebook. The final notebook must contain **no** unresolved flags or "Open questions" cells.

## 9. Quality checklist

- [ ] One notebook per figure that has a `fig_*.md` instruction file.
- [ ] Notebook begins with a big title markdown cell `# Figure N` (e.g. `# Figure 1`); first code cell = Panel A; one code cell per panel; order A, B, C, ...
- [ ] Every cell output begins with the label `<figure_number>.<UPPERCASE_LETTER>` (e.g. `1.A`).
- [ ] Every cell output includes the rendered panel (static image embedded).
- [ ] Panel letter drawn bold uppercase, top-left, ABOVE the axes.
- [ ] All code cells carry `{"jupyter": {"source_hidden": true}}` + `"tags": ["hide-input"]`.
- [ ] HTML conversion shows no `In [n]:` prompts (code hidden).
- [ ] Values match the figure instruction file (verified after execution).
- [ ] Cells are self-contained (imports + loading + cleaning + computation + plotting).
- [ ] Notebook contains no "Open questions" / flag / note cells (production-ready); every ambiguity was raised to the user, resolved, and applied; nothing guessed.

## 10. Worked example — Figure 1 (fig_patients_demographics)

Source spec: `context/figure_instructions/fig_patients_demographics.md` (5 panels A–E). Notebook: `notebooks/fig_patients_demographics.ipynb`.

| Cell | Label output | Panel (from fig spec) |
|------|--------------|-----------------------|
| 1    | `1.A`        | Patient Data — Missingness bar (13 cols; EFS_MONTHS 28.4%, OS_MONTHS 27.0% highlighted) |
| 2    | `1.B`        | Overall Survival Status donut (LIVING 1875 / 65.3%, DECEASED 640 / 22.3%, Unknown 355 / 12.4%) |
| 3    | `1.C`        | Event-Free Survival Status — EFS Binary + EFS Detailed bars (No Event 1286, Event 1222, Unknown 362; detailed counts) |
| 4    | `1.D`        | Age Distribution histogram, 40 bins (mean 9.4, median 8.0, min 0.0, max 73.0, missing 58) |
| 5    | `1.E`        | AGE outliers per cancer group — horizontal bar; Tukey IQR outlier (AGE > 27.125 yr); Oligodendroglioma 53, High-grade glioma 3, six groups ×1; annotate "53/62 (85.5%)" |

The notebook opens with the title markdown cell `# Figure 1`. Each panel cell is self-contained, prints its label first, draws the panel per the fig spec, and carries hide-input metadata. The delivered notebook contains only the title cell and the five panel cells — no "Open questions" cell (all flags were resolved with the user and applied). Panel E uses the resolved Tukey IQR definition (AGE > 27.125) recorded in the fig instruction file.
