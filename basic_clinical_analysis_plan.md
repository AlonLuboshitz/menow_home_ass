# Basic Clinical Analysis Plan — PBTA_RNA

## Overview

The PBTA_RNA dataset contains three files:

| File | Type | Rows | Key Identifier |
|------|------|------|----------------|
| `data_clinical_patient_attributes.txt` | **Patient-level clinical data** | ~2,871 patients | `PATIENT_ID` |
| `data_clinical_sample_attributes.txt` | **Sample-level clinical data** | ~4,313 samples | `SAMPLE_ID` |
| `data_mrna_expression_continuous_rna_seq_v2_mrna.txt` | **Gene expression matrix** | ~37,679 genes × ~2,516 samples | `Hugo_Symbol` |

The two clinical files are linked by `PATIENT_ID`. Each patient can contribute multiple samples (e.g., tumor progression over time).

---

## 1. data_clinical_patient_attributes.txt — Patient Demographics & Outcomes

### Format

The file uses a **4-line metadata header** (lines 1–4) describing each column, then a column header row (line 5), followed by data rows.

- **Line 1 (#):** Column description (human-readable)
- **Line 2 (#):** Long description / definition
- **Line 3 (#):** Data type (`STRING`, `NUMBER`)
- **Line 4 (#):** Number of unique values per column (diagnostic)
- **Line 5:** Machine-readable column names (e.g., `PATIENT_ID`, `AGE`)
- **Lines 6+:** Data

### Columns

| Column | Type | Description |
|--------|------|-------------|
| `PATIENT_ID` | STRING | Unique patient identifier |
| `AGE` | NUMBER | Age at diagnosis (years) |
| `AGE_IN_DAYS` | NUMBER | Age at diagnosis (days) |
| `CANCER_PREDISPOSITIONS` | STRING | Known germline cancer predispositions |
| `EFS_MONTHS` | NUMBER | Event-free survival (months since initial treatment) |
| `EFS_STATUS` | STRING | Event-free survival status (event type) |
| `ETHNICITY` | STRING | Ethnicity |
| `EXTERNAL_PATIENT_ID` | STRING | Original/external patient ID |
| `GERMLINE_SEX_ESTIMATE` | STRING | Estimated germline sex |
| `OS_MONTHS` | NUMBER | Overall survival (months since initial diagnosis) |
| `OS_STATUS` | STRING | Overall survival status |
| `RACE` | STRING | Race |
| `SEX` | STRING | Patient sex |

### Key Value Summaries

#### Categorical Columns

**CANCER_PREDISPOSITIONS** (top values):
| Value | Count | Notes |
|-------|-------|-------|
| `None documented` | ~2,557 | Majority — no known predisposition |
| `Neurofibromatosis, Type 1 (NF-1)` | ~94 | Most common predisposition |
| `Not Reported` | ~68 | |
| `Other inherited conditions NOS` | ~53 | |
| `Li-Fraumeni syndrome (TP53)` | ~32 | |
| `Neurofibromatosis, Type 2 (NF-2)` | ~20 | |
| `Tuberous Sclerosis (TSC1, TSC2)` | ~14 | |
| Others | ~33 | Rare syndromes (VHL, Lynch, CHEK2, etc.) |

Some patients have **multiple predispositions** (comma-separated).

**OS_STATUS** (overall survival):
| Value | Count |
|-------|-------|
| `0:LIVING` | 1,875 (~74%) |
| `1:DECEASED` | 640 (~25%) |
| (blank) | 355 |

**EFS_STATUS** (event-free survival):
| Value | Count |
|-------|-------|
| `0:No Event` | ~1,286 |
| `1:Progressive` | ~423 |
| `1:NA` | ~362 |
| `1:Recurrence` | ~305 |
| `1:Progressive - Metastatic` | ~177 |
| `1:Deceased-due to disease` | ~158 |
| `1:Recurrence - Metastatic` | ~88 |
| `1:Second Malignancy` | ~45 |
| Other death causes | ~26 |

**RACE**:
| Value | Count |
|-------|-------|
| `White` | ~1,724 |
| (blank) | ~744 |
| `Black or African American` | ~208 |
| `Asian` | ~102 |
| `Reported Unknown` | ~31 |
| `More Than One Race` | ~22 |
| `Other` | ~13 |
| `American Indian or Alaska Native` | ~11 |
| `Native Hawaiian or Other Pacific Islander` | ~9 |
| `Not Reported` | ~6 |

**SEX**: `Male` / `Female` (plus blanks).

**ETHNICITY**: `Not Hispanic or Latino` (majority), `Hispanic or Latino`, blanks, `Reported Unknown`.

**GERMLINE_SEX_ESTIMATE**: `Male`, `Female`, `Unknown`, blanks.

#### Numeric Columns

| Column | Range | Notes |
|--------|-------|-------|
| `AGE` | 0–65 years | Some blanks; includes pediatric through adult |
| `AGE_IN_DAYS` | 1–24,106 days | Matches AGE; some blanks |
| `OS_MONTHS` | 0–291 months | Many blanks (censored); tied to OS_STATUS |
| `EFS_MONTHS` | 0–201 months | Many blanks; tied to EFS_STATUS |

### Suggested Summarization Steps

1. **Counts & missingness**: Count patients per column; flag columns with >10% blanks.
2. **Age distribution**: Histogram of AGE; summary stats (mean, median, range, IQR).
3. **Sex balance**: Bar chart of SEX × GERMLINE_SEX_ESTIMATE; check concordance.
4. **Race/Ethnicity**: Bar charts; note high proportion of blanks.
5. **Survival outcomes**:
   - Counts by OS_STATUS and EFS_STATUS.
   - KM curves for OS and EFS (create \"event\" flags from STATUS columns).
6. **Cancer predispositions**: Frequency table; group rare syndromes; flag patients with multiple predispositions.

---

## 2. data_clinical_sample_attributes.txt — Tumor Sample Annotations

### Format

Same 4-line metadata header + header row + data.

### Columns

| Column | Type | Description |
|--------|------|-------------|
| `PATIENT_ID` | STRING | Link to patient file |
| `SAMPLE_ID` | STRING | Unique sample identifier |
| `BROAD_HISTOLOGY` | STRING | High-level histological category |
| `CANCER_GROUP` | STRING | Cancer group classification |
| `CANCER_TYPE` | STRING | Study-defined cancer type |
| `CANCER_TYPE_DETAILED` | STRING | Detailed cancer type |
| `CBTN_TUMOR_TYPE` | STRING | CBTN-assigned tumor type (initial, progressive, recurrence, etc.) |
| `CNS_REGION` | STRING | Anatomical CNS region |
| `COLLECTION_EVENT_ID` | STRING | Links sequential samples from same biological event |
| `EXPERIMENT_STRATEGY` | STRING | Sequencing strategies (e.g., WGS;RNA-Seq) |
| `EXTENT_OF_TUMOR_RESECTION` | STRING | Surgical resection extent |
| `MATCHED_NORMAL_SAMPLE_ID` | STRING | Normal sample used for comparison |
| `MATCHED_NORMAL_SPECIMEN_ID` | STRING | Normal specimen identifier |
| `MOLECULAR_SUBTYPE` | STRING | Molecular subtype per WHO 2016 classification |
| `ONCOTREE_CODE` | STRING | OncoTree code for cancer type |
| `PATHOLOGY_FREE_TEXT_DIAGNOSIS` | STRING | Free-text pathological diagnosis |
| `RNA_LIBRARY_SELECTION` | STRING | RNA library preparation method |
| `SAMPLE_TYPE` | STRING | Solid Tissue, Derived Cell Line, etc. |
| `SPECIMEN_ID` | STRING | KFDRC tumor biospecimen IDs (may be multiple, `;` separated) |
| `SUB_COHORT` | STRING | Study sub-cohort |
| `TUMOR_FRACTION` | NUMBER | Estimated tumor fraction |
| `TUMOR_PLOIDY` | NUMBER | Estimated tumor ploidy |
| `TUMOR_TISSUE_SITE` | STRING | Anatomical tumor location |
| `TUMOR_TYPE` | STRING | Primary, metastatic, recurrence, progression, etc. |

### Key Value Summaries

#### Categorical Columns

**BROAD_HISTOLOGY** (top):
| Value | Approx. Count |
|-------|---------------|
| `Low-grade astrocytic tumor` | ~1,000 |
| `Diffuse astrocytic and oligodendroglial tumor` | ~700 |
| `Embryonal tumor` | ~600 |
| `Ependymal tumor` | ~400 |
| `Tumor of cranial and paraspinal nerves` | ~200 |
| Others | Meningioma, Neuronal & mixed, etc. |

**CANCER_GROUP** (top):
| Value | Count |
|-------|-------|
| `Low-grade glioma` | ~862 |
| `High-grade glioma` | ~512 |
| `Medulloblastoma` | ~440 |
| `Diffuse midline glioma` | ~421 |
| `Ependymoma` | ~341 |
| Others | See full data |

**TUMOR_TYPE**:
| Value | Count |
|-------|-------|
| `primary` | ~2,456 |
| `metastatic` | ~532 |
| `progression` | ~458 |
| `recurrence` | ~328 |
| `Deceased` | ~190 |
| `Primary Tumor` | ~176 |
| `Diagnosis` | ~83 |
| Other variants | ~90 |

**CNS_REGION**:
| Value | Count |
|-------|-------|
| `Hemispheric` | ~1,017 |
| `Posterior fossa` | ~866 |
| `Mixed` | ~763 |
| `Midline` | ~508 |
| `Other` | ~379 |
| `Ventricles` | ~220 |
| `Spine` | ~219 |
| `Suprasellar` | ~174 |

**MOLECULAR_SUBTYPE**: Contains detailed WHO 2016 subtypes, e.g.:
- `HGG, H3 wildtype`
- `LGG, KIAA1549-BRAF`
- `MB, Group3` / `MB, Group4` / `MB, SHH`
- `DMG, H3 K28`
- `EPN, PF A` / `EPN, ST ZFTA`
- `CNS Embryonal, NOS`
- Many blanks / `To be classified`

**EXPERIMENT_STRATEGY**: Combinations like `WGS;RNA-Seq`, `WXS;RNA-Seq`, `Targeted Sequencing`, `RNA-Seq`, `Fusion_Panel`, etc.

**SAMPLE_TYPE**: Primarily `Solid Tissue`, some `Derived Cell Line`.

**SUB_COHORT**: `CBTN` (majority), `MiOncoSeq`, `DGD`, `PNOC`, `Oligo Nation`.

#### Numeric Columns

| Column | Range | Notes |
|--------|-------|-------|
| `TUMOR_FRACTION` | ~0.02–1.0 | Estimated purity; many values near 1.0; some blanks |
| `TUMOR_PLOIDY` | 2–4 | Integer-like but continuous; some blanks |

### Suggested Summarization Steps

1. **Sample counts per patient**: Histogram of samples/patient to see how many patients have multiple timepoints.
2. **Cancer type distribution**: Bar chart of CANCER_GROUP and BROAD_HISTOLOGY.
3. **Tumor type breakdown**: Primary vs. metastatic vs. recurrence vs. progression counts.
4. **Anatomical distribution**: CNS_REGION and TUMOR_TISSUE_SITE frequencies.
5. **Molecular subtypes**: Frequency table; note proportion of unclassified / blanks.
6. **Tumor purity**: Histogram of TUMOR_FRACTION; check relationship with TUMOR_TYPE.
7. **Ploidy**: Distribution; check relationship with cancer group.
8. **Sequencing strategy**: EXPERIMENT_STRATEGY breakdown.
9. **Link back to patients**: Merge with patient data on PATIENT_ID for survival-by-subtype analyses.

---

## 3. Data Relationships & Integration Plan

```
PATIENT_LEVEL                    SAMPLE_LEVEL
┌─────────────────────┐          ┌──────────────────────────┐
│ PATIENT_ID          │◄────────┤ PATIENT_ID               │
│ AGE                 │   1:N   │ SAMPLE_ID                │
│ SEX                 │          │ CANCER_GROUP             │
│ RACE                │          │ MOLECULAR_SUBTYPE        │
│ OS_MONTHS / STATUS  │          │ TUMOR_FRACTION / PLOIDY  │
│ EFS_MONTHS / STATUS │          │ TUMOR_TYPE               │
│ CANCER_PREDISPOSITION│          │ CNS_REGION              │
└─────────────────────┘          └──────────┬───────────────┘
                                             │
                                             │ SAMPLE_ID matches
                                             │ column names in
                                             ▼
                              ┌──────────────────────────────┐
                              │ mRNA Expression Matrix       │
                              │ Genes (rows) × Samples (cols)│
                              └──────────────────────────────┘
```

### Key Analysis Questions (Clinical)

1. **Survival by cancer group**: OS and EFS stratified by CANCER_GROUP / MOLECULAR_SUBTYPE.
2. **Predisposition prevalence**: How many patients have known cancer predisposition syndromes?
3. **Demographics**: Age at diagnosis, sex, race/ethnicity distributions across cancer types.
4. **Tumor purity**: Does TUMOR_FRACTION vary by CANCER_GROUP or TUMOR_TYPE?
5. **Event patterns**: What fraction of samples are from progressive vs. recurrent vs. primary disease?

---

## 4. Known Quality Notes

- **Blank cells exist** across both clinical files — handle missing data consistently.
- **Some columns have semantic prefixes** (e.g., `0:LIVING`, `1:DECEASED`) — strip prefix or parse as label.
- **MOLECULAR_SUBTYPE** has many entries like `To be classified` — treat as missing/uncategorized.
- **SPECIMEN_ID** sometimes contains multiple IDs separated by `;` — may need splitting.
- **CANCER_PREDISPOSITIONS** can list multiple syndromes comma-separated — consider splitting for one-hot encoding.
- **TUMOR_FRACTION** and **TUMOR_PLOIDY** have blanks — investigate missingness pattern.
