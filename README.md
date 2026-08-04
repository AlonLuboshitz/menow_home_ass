# PBTA RNA Clinical Data Analysis

Exploratory analysis of patient- and sample-level clinical data from a Pediatric Brain Tumor Atlas (`PBTA_RNA`) dataset.

The project uses a Jupyter notebook to profile the cohort, clean clinical variables, visualize tumor and demographic characteristics, and compare overall survival (OS) and event-free survival (EFS) across cancer groups and molecular subtypes. The current analysis is focused on clinical metadata; the RNA-expression matrix is documented but is not yet analyzed.

## Project goals

This repository addresses the clinical-data portion of the PBTA exploratory challenge:

1. Explore the available pediatric brain tumor datasets.
2. Characterize patient demographics, outcomes, cancer predispositions, and tumor samples.
3. Integrate patient- and sample-level records using `PATIENT_ID`.
4. Identify cancer groups or molecular subtypes suitable for a focused biological or clinical question.
5. Produce reproducible, interactive figures and a concise summary of the cohort.

## Analysis overview

The main notebook contains 20 numbered steps plus Step 9a. The workflow covers:

### Patient-level exploration

- Dataset dimensions, column types, and missingness
- Age, sex, race, and ethnicity distributions
- OS and EFS status summaries
- Kaplan-Meier survival curves
- Cancer-predisposition prevalence and demographic patterns

### Sample-level exploration

- Cancer group, cancer type, broad histology, and CNS-region distributions
- Tumor type and longitudinal sampling patterns
- Tumor fraction and tumor ploidy
- Molecular-subtype landscape
- Sequencing strategies and RNA-library preparation
- Patients represented by multiple cancer groups or multiple samples

### Integrated analyses

- Patient-to-sample merge validation
- OS and EFS by major cancer group
- OS and EFS by molecular subtype
- Age and sex distributions across cancer groups
- Tumor purity by cancer group and tumor type
- Predisposition status versus age and outcome
- CNS-region versus cancer-group associations
- Consolidated figure inventory and generated Markdown summary

The notebook produces 30+ interactive Plotly figures and uses statistical tests including Kruskal-Wallis, Mann-Whitney U, chi-square, and log-rank comparisons.

## Repository structure

```text
menow_home_ass/
├── context/
│   ├── assignment.md
│   ├── clinical_analysis_instruction.md
│   ├── clinical_analysis_plan.md
│   ├── code_guidelines.md
│   └── referencess.md
├── notebooks/
│   ├── clinical_analysis.ipynb
│   └── clinical_analysis/
│       ├── basic_clinical_summary.md
│       └── nb_cells/
│           ├── imports.py
│           ├── build_notebook.py
│           ├── step_01.md / step_01.py / ...
│           └── step_20.md / step_20.py
├── .gitignore
├── opencode.json
└── requirements.txt
```

### Important files

| File | Description |
|---|---|
| `notebooks/clinical_analysis.ipynb` | Main executable analysis notebook with saved outputs. |
| `notebooks/clinical_analysis/nb_cells/imports.py` | Shared imports, data loaders, cleaning functions, Kaplan-Meier helpers, and log-rank functions. |
| `notebooks/clinical_analysis/nb_cells/step_*.py` | Source code for individual analysis steps. |
| `notebooks/clinical_analysis/nb_cells/step_*.md` | Markdown headings and descriptions for notebook steps. |
| `notebooks/clinical_analysis/basic_clinical_summary.md` | Generated high-level cohort summary. |
| `context/clinical_analysis_instruction.md` | Detailed notebook specification and missing-value rules. |
| `context/clinical_analysis_plan.md` | Data schema, field descriptions, and suggested analyses. |
| `context/code_guidelines.md` | Environment, style, validation, and version-control conventions. |

## Data

The raw data are intentionally excluded from version control. Create a `PBTA_RNA/` directory at the repository root and add the following files:

```text
PBTA_RNA/
├── data_clinical_patient_attributes.txt
├── data_clinical_sample_attributes.txt
└── data_mrna_expression_continuous_rna_seq_v2_mrna.txt
```

Only the two clinical files are required for the current notebook. The expression matrix is reserved for future transcriptomic analysis.

### Expected file format

The two clinical files use the cBioPortal-style tab-separated format:

- Lines 1-4: metadata describing the columns
- Line 5: machine-readable column names
- Lines 6 onward: data records

They are therefore loaded with:

```python
pd.read_csv(path, sep="\t", header=4)
```

### Main identifiers

- `PATIENT_ID` uniquely identifies patients and links the two clinical tables.
- `SAMPLE_ID` uniquely identifies tumor samples.
- A patient may have more than one sample or collection event.

## Installation

### Prerequisites

- Python 3.11 recommended; Python newer than 3.10 is expected
- [`uv`](https://docs.astral.sh/uv/)
- JupyterLab or Jupyter Notebook

### Set up the environment

```bash
git clone https://github.com/AlonLuboshitz/menow_home_ass.git
cd menow_home_ass

uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

Register the environment as a Jupyter kernel:

```bash
python -m ipykernel install --user \
  --name pbta_env \
  --display-name "Python (PBTA)"
```

## Configure local paths

Several source cells currently use the original development path:

```text
/home/alon/menow_home_ass
```

When the repository is stored elsewhere, update the following constants before rebuilding or running source cells:

- `DATA_DIR` in `notebooks/clinical_analysis/nb_cells/imports.py`
- `CELLS_DIR` and the output path in `notebooks/clinical_analysis/nb_cells/build_notebook.py`
- The summary-output path in `notebooks/clinical_analysis/nb_cells/step_19.py`

For example:

```python
DATA_DIR = "/absolute/path/to/menow_home_ass/PBTA_RNA"
```

The committed notebook may also contain these absolute paths because the source code is embedded in notebook cells. Search for `/home/alon/menow_home_ass` in `notebooks/clinical_analysis.ipynb` and replace it with your local repository path when needed.

## Run the analysis

Start Jupyter from the repository root:

```bash
jupyter lab
```

Open:

```text
notebooks/clinical_analysis.ipynb
```

Select the **Python (PBTA)** kernel and run the notebook from top to bottom.

### Non-interactive execution

To execute the notebook and save a separate completed copy:

```bash
jupyter nbconvert \
  --to notebook \
  --execute notebooks/clinical_analysis.ipynb \
  --ExecutePreprocessor.timeout=600 \
  --output clinical_analysis_executed.ipynb
```

Interactive Plotly charts are best viewed in JupyterLab or by opening the executed notebook in a compatible notebook viewer.

## Rebuild the notebook from source cells

The modular notebook sources are stored under:

```text
notebooks/clinical_analysis/nb_cells/
```

After correcting the paths in `build_notebook.py`, rebuild with:

```bash
python notebooks/clinical_analysis/nb_cells/build_notebook.py
```

Review the generated notebook path printed by the script before replacing the committed notebook.

## Data-cleaning conventions

The notebook applies documented cleaning rules consistently:

- `RACE` and `ETHNICITY`: missing or unknown-style values are grouped as `Unknown`.
- `OS_STATUS`: prefixes such as `0:` and `1:` are removed; death is encoded as an event.
- `EFS_STATUS`: represented both as a binary event indicator and as detailed event categories.
- `CANCER_PREDISPOSITIONS`: `None documented` becomes `No predisposition`; missing or unreported values become `Unknown`.
- `MOLECULAR_SUBTYPE`: missing values and `To be classified` become `Unclassified`.
- Missing tumor fraction and ploidy values are retained as explicit unknown groups where appropriate.
- Records lacking complete time/status pairs are excluded from the corresponding Kaplan-Meier analysis.

Each notebook section includes validation output for missingness or data integrity.

## Current results snapshot

The committed generated summary reports:

- 2,870 patients
- 4,312 samples
- 55 cancer groups
- Low-grade glioma as the largest cancer group, with 862 samples
- Median age at diagnosis of 8 years
- 2,096 patients with complete OS information
- 313 patients classified as having a known cancer predisposition
- 882 samples with an unclassified molecular subtype
- Approximately one third of samples missing tumor-fraction or tumor-ploidy values

These values reflect the dataset version used to execute the committed notebook and may change with a different or updated export.

## Methods and dependencies

Core analysis packages include:

- `pandas` and `numpy` for data manipulation
- `plotly` for interactive visualizations
- `scipy` for non-parametric and contingency-table tests
- JupyterLab and `nbconvert` for notebook development and execution

Kaplan-Meier estimates and log-rank tests are implemented in project helper functions rather than through a dedicated survival-analysis package.

## Interpretation notes

This repository is intended for exploratory research and assignment work, not clinical decision-making.

Important considerations include:

- Clinical variables contain substantial and non-random missingness.
- Patients may contribute multiple samples; patient-level analyses should avoid counting the same patient more than once unless repeated sampling is part of the question.
- Cancer-group and molecular-subtype classifications may differ across collection events.
- Survival comparisons are unadjusted and do not account for possible confounders such as age, treatment, disease stage, or cohort source.
- Statistical significance should be interpreted together with group size, missingness, and biological plausibility.

## Development conventions

The repository follows these conventions:

- Use `uv pip install`, rather than installing packages inside notebooks.
- Keep raw PBTA data out of Git.
- Prefer Plotly figures with titles, axis labels, and legends.
- Include validation checks in each analysis step.
- Use Conventional Commit prefixes such as `feat:`, `fix:`, `docs:`, and `chore:`.
- Run the notebook end to end before committing analytical changes.

## Possible next steps

- Select a focused `CANCER_GROUP` or `CANCER_TYPE` for the final exploratory question.
- Integrate clinical metadata with the RNA-expression matrix.
- Test pathway or gene-expression differences between relevant subtypes.
- Build and validate a predictive model using patient-level train/test splits.
- Replace hard-coded absolute paths with a portable configuration or project-root discovery mechanism.
- Add automated notebook execution in continuous integration.

## Acknowledgment

The analysis is based on pediatric brain tumor data organized under the PBTA/Open Pediatric Brain Tumor Atlas ecosystem and exported in a cBioPortal-compatible format.
