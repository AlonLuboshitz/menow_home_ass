# Instruction for Building `clinical_analysis.ipynb`

## Context

You are building a Jupyter notebook that performs **exploratory clinical analysis** on pediatric brain tumor datasets from the **PBTA_RNA** study. The ultimate goal (per `assignment.md`) is to select a pediatric brain tumor type and investigate a biological/clinical question — but this notebook focuses purely on **clinical data exploration** (no mRNA expression yet).

## What to Read First

Before writing any code, thoroughly read these files:

1. **`/home/alon/menow_home_ass/basic_clinical_analysis_plan.md`** — Contains the full schema of both clinical datasets, value summaries, column descriptions, and analysis suggestions. This is your primary reference.
2. **`/home/alon/menow_home_ass/referencess.md`** — Describes the statistical methodology used in a related TCGA glioma paper (Kruskal-Wallis, Mann-Whitney, Spearman correlations, KM curves, normality testing). Use these as guidance for statistical tests in the notebook.
3. **`/home/alon/menow_home_ass/assignment.md`** — Explains the broader challenge (not directly needed for this notebook, but good for context).

## Data Location

All data is under `/home/alon/menow_home_ass/PBTA_RNA/`:

| File | Description |
|------|-------------|
| `data_clinical_patient_attributes.txt` | Patient-level clinical data (13 columns, ~2,871 rows) |
| `data_clinical_sample_attributes.txt` | Sample-level clinical data (24 columns, ~4,313 rows) |
| `data_mrna_expression_continuous_rna_seq_v2_mrna.txt` | Gene expression matrix (NOT used in this notebook) |

## File Format

Both clinical files share a **4-line metadata header** before the actual column header:

```
Line 1 (#): Short column description
Line 2 (#): Long description / definition
Line 3 (#): Data type (STRING / NUMBER)
Line 4 (#): Number of unique values per column (diagnostic metadata)
Line 5:     Machine-readable column names (e.g., PATIENT_ID, AGE)
Lines 6+:   Data rows (tab-separated)
```

When reading with pandas, skip the first 4 rows (metadata) and use row 5 as the header. The row number for header in pandas would be `header=4` (0-indexed, so row index 4 is the 5th row).

Example:
```python
import pandas as pd
patients = pd.read_csv(
    'PBTA_RNA/data_clinical_patient_attributes.txt',
    sep='\t',
    header=4,
    dtype={'AGE': float, 'AGE_IN_DAYS': float, 'OS_MONTHS': float, 'EFS_MONTHS': float}
)
```

## Missing Value Strategy (User-Approved)

The user has decided the following. Apply these **consistently throughout the notebook**.

| Column(s) | Strategy |
|-----------|----------|
| `AGE` | Leave as NaN; exclude from age-specific plots individually |
| `AGE_IN_DAYS` | Same as AGE |
| `RACE` | Merge blank + `Not Reported` + `Reported Unknown` → `Unknown` category |
| `ETHNICITY` | Merge blank + `Reported Unknown` → `Unknown` category |
| `OS_MONTHS` / `OS_STATUS` | Bar charts include an `Unknown` category for rows missing either time or status; **KM curves** exclude rows without a complete `(time, status)` pair |
| `EFS_MONTHS` / `EFS_STATUS` | Same as OS |
| `CANCER_PREDISPOSITIONS` | Merge blank + `Not Reported` → `Unknown`; `None documented` stays `No predisposition` |
| `TUMOR_FRACTION` | Show missing as separate `Unknown` group in boxplots |
| `TUMOR_PLOIDY` | Show missing as separate `Unknown` group in boxplots |
| `MOLECULAR_SUBTYPE` | Merge blank + `To be classified` → `Unclassified` |
| `SPECIMEN_ID` | Split multi-value entries on `;` only if needed |
| `SEX` | Keep blanks as NaN |
| All other columns | Drop from specific analyses where missing; document in text |

## Data Cleaning Rules

1. **Parse OS_STATUS labels**: Strip `0:` / `1:` prefixes. Create:
   - `os_event`: 1 if `DECEASED`, 0 if `LIVING`, NaN if blank
   - `os_label`: clean label (`LIVING`, `DECEASED`, or NaN)
2. **Parse EFS_STATUS labels**: Strip `0:` / `1:` prefixes. Create TWO representations:
   - **Binary**: `efs_event` = 1 if any event occurred (`Progressive`, `Recurrence`, `Deceased-due to disease`, etc.), 0 if `No Event`, NaN if `NA` or blank
   - **Detailed**: `efs_detail` = keep the specific event type as a clean label (`Progressive`, `Recurrence`, `Deceased-due to disease`, `Progressive - Metastatic`, etc.), `No Event`, or NaN
   - Throughout the notebook, allow the user to toggle between binary event and detailed categories where relevant
3. **Plotting convention**: When showing categorical plots with an `Unknown` category, use a muted/dashed style to distinguish from real data.

## Interactive Plotting

Wherever possible, use **Plotly** (not static matplotlib) to create interactive plots. This allows the user to:
- Hover for exact values
- Click legend items to toggle categories
- Zoom/pan

```python
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
```

For **KM curves**, use `plotly.figure_factory.create_km_curve` or build custom interactive KM plots with `plotly.graph_objects`.

If Plotly is not available, default to seaborn/matplotlib but prefer interactive.

For the **age histogram by different scenarios** (Step 2), use a Plotly histogram with dropdown buttons to switch between stratifications (by SEX, by RACE, by CANCER_GROUP, overall) — all in one interactive plot.

## Notebook Structure & Steps

The notebook is organized into **20 independent steps** (note: a new Step 9a was added for multi-cancer-group analysis). Each step:
1. Has a markdown cell explaining the purpose
2. Has a validation cell documenting missing-value handling for that step specifically
3. Produces at least one plot
4. Is independently runnable (repeats imports/loading as needed)
5. At the end there should be a final cell printing a **table summarizing all steps** with a mapping of which figures they produced

Here are the detailed step-by-step instructions:

---

### Phase 1 — Independent Dataset Summarization

#### Step 1: Load & Profile Patient Data

**What it shows:** A comprehensive overview of the patient dataset — size, column types, missingness levels.

**Commands/Code:**
1. Read `data_clinical_patient_attributes.txt` with `header=4`, `sep='\t'`
2. Print shape (rows × columns)
3. Print column names, dtypes, non-null counts
4. Generate a **missingness bar plot** using Plotly: for each column, show fraction of non-null values vs missing
5. Print a summary table with: column name, dtype, # non-null, # null, % null, unique values (for strings), min/max (for numbers)

**Validation cell below the plot:**
```python
# Validation: Missing values
AGE_nulls = df['AGE'].isna().sum()
print(f"AGE has {AGE_nulls} missing values ({AGE_nulls/len(df):.1%}) — these will be excluded from age-specific plots")
```

---

#### Step 2: Patient Demographics

**What it shows:** Age distribution, sex balance, race/ethnicity makeup of the cohort.

**Commands/Code:**
1. **Interactive age histogram** with dropdown:
   - Default: overall AGE histogram with KDE overlay
   - Dropdown options: "By SEX", "By RACE", "By CANCER_GROUP" (note: CANCER_GROUP is in the sample file, so only available after merge — initially show just SEX)
   - Use Plotly with `px.histogram` and `updatemenus`
2. **SEX bar chart** (counts)
3. **RACE bar chart** (with merged `Unknown` category)
4. **ETHNICITY bar chart** (with merged `Unknown` category)

Arrange as a 2×2 grid using `plotly.subplots.make_subplots`.

**Validation:**
```python
# Validation: RACE merging
print(f"RACE unique before: {sorted(df['RACE'].dropna().unique())}")
# After merging blanks and 'Not Reported'/'Reported Unknown' → 'Unknown'
```

---

#### Step 3: Patient Survival Overview

**What it shows:** Outcome overview and Kaplan-Meier survival curves.

**Commands/Code:**
1. **OS_STATUS pie chart** — shows proportion LIVING vs DECEASED vs Unknown
2. **EFS_STATUS bar chart** — TWO modes (add a toggle or show side-by-side):
   - Mode A (binary): Event vs No Event vs Unknown
   - Mode B (detailed): All event categories (Progressive, Recurrence, Deceased-due to disease, etc.) + No Event + Unknown
3. **KM curve for OS** using `plotly.figure_factory.create_km_curve`:
   - Event indicator: `os_event` (1 if DECEASED, 0 if LIVING)
   - Time: `OS_MONTHS`
   - Exclude rows where either is NaN
   - Include risk table below
4. **KM curve for EFS**: same approach with `efs_event` and `EFS_MONTHS`

**Validation:**
```python
# Validation: OS missing
complete_os = df[['OS_MONTHS', 'os_event']].dropna()
print(f"OS: {len(complete_os)}/{len(df)} patients have complete (time, status) — usable for KM")
print(f"OS status distribution:\n{df['OS_STATUS'].value_counts(dropna=False)}")

# Validation: EFS binary vs detailed
print(f"EFS binary event counts:\n{df['efs_event'].value_counts(dropna=False)}")
print(f"EFS detailed counts:\n{df['efs_detail'].value_counts(dropna=False)}")
```

---

#### Step 4: Cancer Predispositions — Prevalence & Demographics

**What it shows:** How common each predisposition is, and demographic patterns (age, sex) across predispositions. Multi-syndrome patients are exploded so each syndrome counts individually.

**Commands/Code:**
1. **Clean & explode predispositions:**
   - `None documented` → `No predisposition`
   - blank + `Not Reported` → `Unknown`
   - For multi-syndrome entries (comma-separated, e.g., `Li-Fraumeni syndrome (TP53),Neurofibromatosis, Type 1 (NF-1)`):
     - Split on `,` and explode so each syndrome becomes its own row
     - This means a patient with two syndromes appears twice in the exploded dataframe
     - Keep track of original patient IDs for deduplication when needed
2. **Prevalence bar chart:**
   - Show percentage of ALL patients carrying each predisposition
   - Horizontal bar chart of top syndromes sorted by frequency
   - Annotate each bar with count and percentage
3. **Interactive exploration plot** using Plotly dropdowns:
   - Create a scatter/boxplot where the user can select a predisposition from a dropdown
   - When selected, show:
     - % of patients with this predisposition (overall and within each cancer group)
     - AGE boxplot for patients WITH vs WITHOUT this predisposition
     - SEX breakdown bar chart for patients WITH this predisposition
   - Include options: "Show all top predispositions" as default overview
4. **Summary table**: Print a DataFrame with columns: Predisposition, Count, % of Patients, Median Age, % Female, Top Cancer Group

**Important note on percentages:** When computing "% of patients", use the original un-exploded patient count as denominator. The exploded version is only for per-syndrome analysis.

**Validation:**
```python
# Validation: Explosion check
print(f"Patients before explosion: {df['PATIENT_ID'].nunique()}")
print(f"Rows after explosion: {len(df_exploded)}")
# Check multi-syndrome patients
multi = df[df['CANCER_PREDISPOSITIONS'].str.contains(',', na=False)]
print(f"Patients with multiple predispositions (comma-separated): {len(multi)}")
```

---

#### Step 5: Load & Profile Sample Data

**What it shows:** Overview of the sample-level dataset.

Same structure as Step 1 but for `data_clinical_sample_attributes.txt` (24 columns).

**Commands/Code:**
1. Read with `header=4`, `sep='\t'`
2. Shape, columns, dtypes
3. **Missingness heatmap** (Plotly heatmap or imshow)
4. Column summary table

**Validation:** Print missingness per column.

---

#### Step 6: Sample Cancer Type Distributions

**What it shows:** Histological and anatomical breakdown of all samples.

**Commands/Code:**
1. **BROAD_HISTOLOGY bar chart** (top categories, sorted)
2. **CANCER_GROUP bar chart** (top categories)
3. **CNS_REGION bar chart**
4. **TUMOR_TYPE bar chart** (primary, metastatic, progression, recurrence, etc.)
   - Group rare tumor types into `Other` if there are too many

Arrange as a 2×2 grid in Plotly.

**Validation:** Note the most common and rarest categories.

---

#### Step 7: Tumor Purity & Ploidy

**What it shows:** Distribution of tumor purity and ploidy across samples.

**Commands/Code:**
1. **TUMOR_FRACTION histogram** with `Unknown` category as a separate bar with dashed outline
2. **TUMOR_PLOIDY histogram** with `Unknown` category
3. **Scatter plot**: TUMOR_FRACTION (x) vs TUMOR_PLOIDY (y), colored by CANCER_GROUP, with hover showing SAMPLE_ID

**Validation:**
```python
print(f"TUMOR_FRACTION missing: {df['TUMOR_FRACTION'].isna().sum()} / {len(df)}")
print(f"TUMOR_PLOIDY missing: {df['TUMOR_PLOIDY'].isna().sum()} / {len(df)}")
```

---

#### Step 8: Molecular Subtype Landscape

**What it shows:** The diversity of molecular subtypes and their relationship to cancer groups.

**Commands/Code:**
1. Merge blank + `To be classified` → `Unclassified`
2. **Horizontal bar chart** of top 20 MOLECULAR_SUBTYPE values
3. **Heatmap**: cross-tabulation of MOLECULAR_SUBTYPE (rows, top 15) × CANCER_GROUP (columns), with counts annotated
   - Normalize by row (show proportion within each subtype)
   - This heatmap is key — Step 13 will reference it to decide which cancer groups to analyze by subtype

**Validation:**
```python
unclass_count = (df['MOLECULAR_SUBTYPE'] == 'Unclassified').sum()
print(f"Unclassified subtypes: {unclass_count} ({unclass_count/len(df):.1%})")
```

---

#### Step 9: Sequencing Strategy & RNA Library

**What it shows:** What sequencing methods and library prep were used.

**Commands/Code:**
1. **EXPERIMENT_STRATEGY bar chart** — split multi-method entries (e.g., `WGS;RNA-Seq`) and show frequency of each strategy or keep as combined
2. **RNA_LIBRARY_SELECTION bar chart**

**Validation:**
```python
print(f"Experiment strategy unique values: {df['EXPERIMENT_STRATEGY'].nunique()}")
print(f"RNA library selection unique values: {df['RNA_LIBRARY_SELECTION'].nunique()}")
```

---

#### Step 9a: Multi-Cancer-Group Analysis (NEW)

**What it shows:** How many patients have samples assigned to different cancer groups, and what patterns emerge.

**Commands/Code:**
1. Group samples by PATIENT_ID and collect the distinct CANCER_GROUPs per patient
2. Count how many distinct CANCER_GROUPs each patient has
3. **Histogram**: number of distinct cancer groups per patient
4. Print and highlight patients with >1 cancer group — show their PATIENT_ID and which groups
5. For patients with exactly 2 cancer groups, create a **heatmap**: cancer group A (rows) vs cancer group B (columns) showing co-occurrence counts
6. **Sankey diagram** (optional, if plotly is available): show flow from primary cancer group to secondary cancer group for multi-group patients

**Validation:**
```python
multi_group = patients_with_counts[patients_with_counts > 1]
print(f"Patients with multiple CANCER_GROUPs: {len(multi_group)} ({len(multi_group)/total_patients:.1%})")
print(f"Max groups per patient: {patients_with_counts.max()}")
```

---

### Phase 2 — Cross-Dataset Integration

#### Step 10: Merge Patient + Sample Data

**What it shows:** How well the two datasets connect, and whether any patients/samples are orphaned.

**Commands/Code:**
1. Load both datasets (re-read or reuse)
2. Left join samples → patients on `PATIENT_ID`
3. Print merged shape
4. Check for unmapped:
   - Samples whose PATIENT_ID has no match in the patient file → print count and sample IDs
   - Patients in the patient file with no matching samples → print count and patient IDs
5. Report counts

**Validation:**
```python
# Correct validation: check PATIENT_ID matching, NOT AGE
sample_pids = set(samples['PATIENT_ID'].unique())
patient_pids = set(patients['PATIENT_ID'].unique())

orphan_samples = sample_pids - patient_pids
orphan_patients = patient_pids - sample_pids

print(f"Samples with no matching patient record: {len(orphan_samples)}")
if orphan_samples:
    print(f"  Orphan sample patient IDs: {orphan_samples}")
print(f"Patients with no matching sample: {len(orphan_patients)}")
if orphan_patients:
    print(f"  Orphan patient IDs: {orphan_patients}")

# Note: a patient CAN have AGE missing but still have valid sample data
# AGE NaN does NOT indicate a merge problem
```

---

#### Step 11: Samples per Patient

**What it shows:** How many patients have single vs. multiple samples (longitudinal tracking).

**Commands/Code:**
1. Group by PATIENT_ID, count SAMPLE_ID
2. **Histogram** of sample counts per patient
3. Print patients with most samples (top 10)

**Validation:**
```python
multi = (sample_counts > 1).sum()
print(f"Patients with >1 sample: {multi} ({multi/len(sample_counts):.1%})")
```

---

#### Step 12: Survival by Cancer Group

**What it shows:** OS and EFS stratified by major cancer groups. **Both OS and EFS must be shown.**

**Important rule for patients with multiple cancer groups:** If a patient has samples in multiple CANCER_GROUPs, they will appear in ALL relevant survival curves (one contribution per sample). The KM curve is at the **sample level**, not the patient level. Each sample's survival time comes from its patient's OS/EFS data.

**Commands/Code:**
1. Select top 6 CANCER_GROUPs by frequency
2. For each, extract survival data from the merged table:
   - Filter to samples in that cancer group
   - Use the patient-level OS_MONTHS/os_event and EFS_MONTHS/efs_event
3. **Interactive KM plot** with OS on top subplot, EFS on bottom subplot:
   - One curve per CANCER_GROUP
   - Include log-rank test p-value annotation (pairwise or global)
   - Risk tables below
4. Color palette: distinct colors for each group
5. Also print a table showing: Cancer Group, N samples, N patients, Median OS (months), Median EFS (months)

**Validation:**
```python
# Show how many patients contribute to multiple curves
multi_group_patients = merged.groupby('PATIENT_ID')['CANCER_GROUP'].nunique()
multi = multi_group_patients[multi_group_patients > 1]
print(f"Patients appearing in multiple cancer group KM curves: {len(multi)}")
print(f"These patients contribute to {multi.sum()} curve entries total")
```

---

#### Step 13: Survival by Molecular Subtype (Global)

**What it shows:** OS and EFS stratified by molecular subtype across the entire cohort (not filtered by cancer group).

**Commands/Code:**
1. Use the molecular subtypes from Step 8 (with `Unclassified` merged)
2. Select top 6 most common molecular subtypes (by sample count)
3. For each, compute KM curves at the sample level
4. **Interactive KM plot**: OS (top) + EFS (bottom), one curve per subtype + Unclassified
5. Log-rank test p-value annotation
6. Risk tables

**After the global plot, use Step 8's heatmap to decide which Cancer Groups to also analyze by subtype:**
- From the heatmap in Step 8, identify cancer groups that have high subtype diversity (e.g., Low-grade glioma has many subtypes, Medulloblastoma has SHH/Group3/Group4/WNT)
- For each such cancer group, create a **separate KM plot** showing subtypes within that group
- Print a note: "Based on Step 8 heatmap, we additionally analyze subtypes within: [list of cancer groups]"

**Validation:**
```python
subtype_counts = merged['MOLECULAR_SUBTYPE'].value_counts()
print(f"Top 10 molecular subtypes:\n{subtype_counts.head(10)}")
# Also check which cancer groups have rich subtype annotation
group_subtype_div = merged.groupby('CANCER_GROUP')['MOLECULAR_SUBTYPE'].nunique().sort_values(ascending=False)
print(f"Cancer groups by subtype diversity:\n{group_subtype_div}")
```

---

#### Step 14: Age at Diagnosis by Cancer Group

**What it shows:** Whether different cancer groups occur at different ages (consistent with the reference paper's methodology).

**Commands/Code:**
1. **Boxplot**: AGE by CANCER_GROUP (top 8 groups)
2. Overlay individual points (jittered)
3. Kruskal-Wallis test annotation
4. Post-hoc pairwise comparisons (Mann-Whitney) for significant results
5. Use the sample-level merged data (so patients with samples in multiple groups appear in each)

**Reference paper approach:** Use the same statistical method described in `referencess.md` — test normality first (Kolmogorov-Smirnov), then use Kruskal-Wallis if non-normal.

**Validation:**
```python
from scipy.stats import kruskal, mannwhitneyu
groups = [group['AGE'].dropna() for name, group in merged.groupby('CANCER_GROUP') if len(group) > 10]
stat, p = kruskal(*groups)
print(f"Kruskal-Wallis: H={stat:.2f}, p={p:.4f}")
```

---

#### Step 15: Sex Balance by Cancer Group

**What it shows:** Whether certain cancer groups show sex bias.

**Commands/Code:**
1. **Stacked bar chart**: SEX (Male/Female) by CANCER_GROUP (top 8)
2. Chi-squared test of independence
3. Annotate p-value on plot

**Validation:**
```python
from scipy.stats import chi2_contingency
ct = pd.crosstab(merged['CANCER_GROUP'], merged['SEX'])
chi2, p, dof, expected = chi2_contingency(ct)
print(f"Chi-squared: χ²={chi2:.2f}, p={p:.4f}")
```

---

#### Step 16: Purity by Cancer Group & Tumor Type

**What it shows:** Tumor purity differences across cancer groups and clinical states.

**Commands/Code:**
1. **Boxplot 1**: TUMOR_FRACTION by CANCER_GROUP (with `Unknown` group shown)
2. **Boxplot 2**: TUMOR_FRACTION by TUMOR_TYPE (primary, metastatic, progression, recurrence)
3. Kruskal-Wallis p-values for both

**Validation:**
```python
print(f"Samples with Unknown TUMOR_FRACTION: {(merged['TUMOR_FRACTION'].isna()).sum()}")
```

---

#### Step 17: Predisposition vs Outcome

**What it shows:** Whether known cancer predisposition affects survival or age of onset.

**Commands/Code:**
1. Create binary flag: `has_predisposition` = True if CANCER_PREDISPOSITIONS is not `No predisposition` and not `Unknown`
2. **KM curves** for OS and EFS: stratified by has_predisposition
3. **Boxplot**: AGE by has_predisposition
4. Log-rank test and Mann-Whitney p-values

**Validation:**
```python
pred_counts = df_patients['has_predisposition'].value_counts()
print(f"Patients with predisposition: {pred_counts.get(True, 0)} / {len(df_patients)}")
```

---

#### Step 18: CNS Region vs Cancer Group

**What it shows:** Anatomical distribution patterns of different cancer types.

**Commands/Code:**
1. Cross-tabulation: CNS_REGION (rows) × CANCER_GROUP (columns)
2. **Heatmap** with row-normalization (each row sums to 1)
3. Annotated with counts in each cell

**Validation:**
```python
region_counts = merged['CNS_REGION'].value_counts()
print(f"CNS_REGION missing: {merged['CNS_REGION'].isna().sum()}")
```

---

### Phase 3 — Summary

#### Step 19: Generate Summary Report

**What it shows:** A markdown summary of all findings from the notebook.

**Commands/Code:**
1. Collect key numbers from each step into a dictionary
2. Write a markdown string with sections:
   - Dataset overview (shapes, missingness rates)
   - Key demographic findings (age range, sex balance)
   - Survival highlights (median OS by top groups)
   - Purity findings
   - Molecular subtype landscape
   - Predisposition summary
3. Create a **table mapping step number → figure title → key insight**
4. Save as `basic_clinical_summary.md`

**Validation:** Print confirmation that the summary was saved.

---

#### Step 20: Summary Table of All Figures (NEW)

**What it shows:** A consolidated overview of every figure produced in the notebook, with a brief insight per figure.

**Commands/Code:**
1. Build a DataFrame with columns: Step, Figure Title, Figure Type, Key Insight, Interactive?
2. Fill it manually based on what was generated
3. Display as a styled table
4. Count total figures

**Validation:** Print total figure count.

---

## Library Requirements

Include a cell at the top of the notebook installing/importing:

```python
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.figure_factory import create_km_curve
from scipy.stats import kruskal, mannwhitneyu, chi2_contingency
import warnings
warnings.filterwarnings('ignore')
```

If any library is missing, the agent should install it via `!pip install library_name`.

## Final Quality Checklist

Before finishing, verify:
- [ ] All 20 steps exist as separate cells/sections
- [ ] Each step has a markdown intro cell explaining the step
- [ ] Each step has a validation cell documenting missing-value handling
- [ ] Each step produces at least one plot
- [ ] All plots have: title, axis labels, legend (if applicable)
- [ ] KM plots always show both OS and EFS (side-by-side or top-bottom)
- [ ] EFS_STATUS has both binary (event/no event) and detailed categories
- [ ] Step 2 age histogram has interactive dropdown
- [ ] Step 4 has interactive predisposition explorer (dropdown to select predisposition)
- [ ] Step 9a (multi-cancer-group analysis) exists between Step 9 and Step 10
- [ ] Step 10 validation uses PATIENT_ID matching, NOT AGE
- [ ] Step 12 includes multi-cancer-group patients in all relevant curves
- [ ] Step 13 uses global molecular subtypes + references Step 8 heatmap for group filtering
- [ ] Step 20 summary table exists
- [ ] Missing value strategy matches the table in this instruction
- [ ] Saved as `clinical_analysis.ipynb`
