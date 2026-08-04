# PBTA Clinical and RNA Analysis

Exploratory analysis of the Pediatric Brain Tumor Atlas (PBTA) clinical dataset, with supporting RNA-expression data, to identify clinically meaningful tumor groups and variables for deeper investigation.

The repository contains the source and executed notebooks, exported result tables, figures, and a short LaTeX report. The current analysis covers cohort characterization, clinical associations, survival outcomes, age-related effects, and multivariate hidden-structure analysis.

## Project objective

The project was developed for an exploratory PBTA data-analysis assignment. Its main goals are to:

1. characterize the available pediatric brain-tumor cohort;
2. compare clinical variables across cancer groups;
3. evaluate overall-survival and event-free-survival patterns;
4. examine age-related differences;
5. explore multivariate structure within sufficiently large diagnostic groups; and
6. select a tumor group and research direction for further study.

The committed summary report prioritizes high-grade glioma (HGG) for follow-up analysis, while also highlighting diffuse midline glioma (DMG), low-grade glioma (LGG), and age as clinically important dimensions of the dataset.

## Repository structure

```text
.
├── PBTA_RNA.zip
├── Cites/
│   └── *.bib
├── Figures/
│   ├── fig_patients_demographics.pdf
│   ├── fig_patients_demographics.png
│   └── additional figure assets
├── build/
│   ├── report_summary.pdf
│   ├── survival_tables.tex
│   └── generated LaTeX outputs
├── context/
│   └── assignment, methodology, and analysis notes
├── notebooks/
│   ├── age_deciles/
│   ├── clinical_analysis/
│   ├── clinical_associations/
│   ├── clinical_mulltivar_hidden_strcture_analysis/
│   ├── figure_notebooks/
│   ├── survival_analysis/
│   └── *_executed.ipynb
├── scripts/
│   └── build_survival_tables.py
├── figure_summary.html
├── report_summary.tex
├── requirements.txt
├── .gitattributes
└── .gitignore
```

The directory name `clinical_mulltivar_hidden_strcture_analysis` is reproduced exactly as it currently appears in the repository.

## Analysis modules

| Module | Purpose | Main outputs |
|---|---|---|
| Basic clinical analysis | Cohort size, sample availability, diagnosis distribution, age, sex, and missingness | Summary Markdown and executed notebook |
| Clinical associations | Tests associations between tumor groups and clinical variables, with multiple-testing correction and effect sizes | Result and summary CSV files |
| Survival analysis | Overall-survival and event-free-survival comparisons, post-hoc testing, and Cox proportional-hazards analysis | Survival tables, summaries, and notebook |
| Age-decile analysis | Examines age distributions and age-associated differences across cancer groups | Notebook and result CSV |
| Hidden-structure analysis | Uses dimensionality reduction and clustering within sufficiently large diagnostic groups | Per-group clustering results and notebook |
| Figure notebooks | Produces report-ready demographic, clinical-association, sample, and survival visualizations | HTML and notebook exports |

## Data and Git LFS

`PBTA_RNA.zip` is approximately 102 MB and is tracked with **Git Large File Storage (Git LFS)**. The `.gitattributes` file associates this specific archive with the LFS filter.

Install Git LFS before cloning or pulling the repository. Without it, Git may retrieve only the small LFS pointer instead of the actual archive.

### Clone and retrieve the data

```bash
git lfs install
git clone https://github.com/AlonLuboshitz/menow_home_ass.git
cd menow_home_ass
git lfs pull
```

Extract the archive in the repository root:

```bash
unzip PBTA_RNA.zip
```

On PowerShell:

```powershell
Expand-Archive PBTA_RNA.zip -DestinationPath .
```

After extraction, the notebooks expect a directory named `PBTA_RNA/`. The unpacked raw text files are intentionally ignored by Git, while the ZIP archive itself remains versioned through Git LFS.

Expected clinical-data paths include:

```text
PBTA_RNA/data_clinical_patient_attributes.txt
PBTA_RNA/data_clinical_sample_attributes.txt
```

Some analyses also use the RNA-expression file contained in the extracted data directory.

## Environment setup

The project uses Python and Jupyter notebooks. Python 3.11 is recommended because that is the version referenced by the existing project environment.

### Using `uv`

```bash
uv venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
uv venv --python 3.11
.venv\Scripts\Activate.ps1
uv pip install -r requirements.txt
```

Register the environment as a Jupyter kernel:

```bash
python -m ipykernel install --user --name menow-home-ass --display-name "Python (menow_home_ass)"
```

### Additional notebook dependencies

The current notebooks and methodology files reference several packages that are not all listed in the committed `requirements.txt`. Until the requirements file is consolidated, install the additional analysis packages explicitly:

```bash
uv pip install statsmodels scikit-posthocs lifelines scikit-learn prince
```

## Running the analyses

Run commands and launch Jupyter from the repository root so that the project-relative data and output paths resolve consistently:

```bash
jupyter lab
```

Open the source notebook for the desired module and select the `Python (menow_home_ass)` kernel. The top level of `notebooks/` also contains executed notebook copies that can be reviewed without rerunning every analysis.

A practical execution order is:

1. `notebooks/clinical_analysis/`
2. `notebooks/clinical_associations/`
3. `notebooks/survival_analysis/`
4. `notebooks/age_deciles/`
5. `notebooks/clinical_mulltivar_hidden_strcture_analysis/`
6. `notebooks/figure_notebooks/`

Some module directories contain a `src/` folder with scripts or Markdown cells used to construct the corresponding notebook.

## Generated outputs

Important committed outputs include:

- `build/report_summary.pdf` — compiled short report;
- `report_summary.tex` — report source;
- `figure_summary.html` — browser-viewable figure summary;
- `Figures/` — exported report figures;
- `build/survival_tables.tex` — generated LaTeX survival tables; and
- module-specific CSV result files under `notebooks/`.

To rebuild the survival tables:

```bash
python scripts/build_survival_tables.py
```

Compile `report_summary.tex` with a local LaTeX installation or an online LaTeX editor after regenerating any required tables and figures.

## Current result snapshot

The committed exploratory analysis reports:

- 2,870 patients;
- 4,312 samples;
- 55 cancer groups;
- low-grade glioma as the largest diagnostic group;
- age as the strongest recurring clinical variable across several analyses; and
- HGG as the primary group selected for deeper follow-up in the summary report.

These findings are exploratory and depend on the available PBTA annotations, missing-data patterns, group sizes, and analysis choices documented in the notebooks.

## Reproducibility notes

- The raw extracted data are not committed directly; retrieve the LFS archive and extract it locally.
- Analysis code uses repository-relative paths; run notebooks and scripts from the repository root.
- The current `requirements.txt` does not yet capture every package imported by all deep-dive notebooks.
- Generated outputs are committed for inspection, but the source notebooks and data should be treated as the reproducible source of truth.
- The project is an exploratory research analysis and is not intended for clinical decision-making.

## License

No `LICENSE` file is currently included in the repository. Add an explicit license before reusing or distributing the project beyond its intended context.
