# mRNA Expression Data Reference

This document describes the structure of the mRNA expression matrix and how it connects to the clinical data in the PBTA_RNA project.

---

## 1. File Location

```
/home/alon/menow_home_ass/PBTA_RNA/data_mrna_expression_continuous_rna_seq_v2_mrna.txt
```

## 2. Matrix Dimensions

| Dimension     | Value  |
| ------------- | ------ |
| Genes (rows)  | 37,680 |
| Samples (cols)| 2,518  |
| - _Hugo_Symbol_ | 1 column (gene name) |
| - _Entrez_Gene_Id_ | 1 column (Entrez ID) |
| - Sample IDs | 2,516 columns |

## 3. Structure

**Rows:** One per gene, identified by `Hugo_Symbol` (e.g., `TP53`, `EGFR`, `BRAF`). All 37,680 symbols are unique — no duplicate genes.

**Columns:** Tab-separated values:
```
Hugo_Symbol  Entrez_Gene_Id  16510-1  16510-10  16510-11  ...  7316-272  7316-2577  ...
```

**Values:** Continuous expression values (likely RSEM-normalized or similar). Example for TP53:
```
TP53  7157  1.43  5.83  14.52  5.9  13.33  ...
```

## 4. How Expression Connects to Clinical Data

The expression matrix column headers (excluding `Hugo_Symbol` and `Entrez_Gene_Id`) are **`SAMPLE_ID`s** — they match directly to the `SAMPLE_ID` column in `data_clinical_sample_attributes.txt`.

**Join path:**
```
Expression columns (sample IDs)
       ↓ (exact match on SAMPLE_ID)
data_clinical_sample_attributes.txt  →  has PATIENT_ID
       ↓ (exact match on PATIENT_ID)
data_clinical_patient_attributes.txt  →  patient-level data
```

### 4.1 Sample Count Discrepancy

| Dataset | Sample Count |
|---|---|
| Clinical sample file | **4,312** |
| Expression matrix | **2,516** |
| Samples with clinical data but NO expression | **1,796 (42%)** |

### 4.2 Why the Gap?

The 1,796 clinical samples without expression data were sequenced using **DNA-based methods**, not RNA-Seq:

| Sequencing Strategy | Missing Count |
|---|---|
| WGS (Whole Genome Sequencing) | 759 |
| Targeted Sequencing / Fusion Panel | 1,021 |
| WXS (Whole Exome Sequencing) | 16 |

The 2,516 samples in the expression matrix all had **RNA-Seq** (with or without complementary WGS).

## 5. Implications for Downstream Analysis

### 5.1 Sample Reduction
Any analysis that uses expression data is limited to **2,516 of 4,312** clinical samples (58%). This has several consequences:

### 5.2 Selection Bias
The 2,516 RNA-Seq samples may not be a random subset. Check whether certain cancer groups, tumor types, or patient demographics are over/under-represented in the expression data compared to the full clinical cohort.

### 5.3 Survival Analysis
- Survival curves can use all 4,312 clinical samples (more power)
- But expression-based survival analysis (e.g., high vs. low TP53) is limited to 2,516 samples
- When merging, document the dropout explicitly

### 5.4 Rare Cancer Groups
Small cancer groups may lose too many samples when restricting to RNA-Seq data. Check minimum sample sizes before doing group-wise expression comparisons.

### 5.5 Paired Sample Analysis
A patient may have multiple samples (e.g., primary + recurrence). When merging:
- A patient is retained if **any** of their samples have expression data
- For patient-level analysis, decide how to handle multiple samples (average? pick one?)
- For sample-level analysis, each sample is independent

### 5.6 No Matched Normals
The expression matrix contains only **tumor samples** (primary, recurrence, metastatic, etc.). There are no normal tissue samples in this file. The `MATCHED_NORMAL_SAMPLE_ID` column in the clinical sample file references a normal sample ID, but its expression data is not in this matrix.

### 5.7 Practical Merging Code

```python
import pandas as pd

# Load data
samples = pd.read_csv("/home/alon/menow_home_ass/PBTA_RNA/data_clinical_sample_attributes.txt", sep="\t", header=4)
patients = pd.read_csv("/home/alon/menow_home_ass/PBTA_RNA/data_clinical_patient_attributes.txt", sep="\t", header=4)
expr = pd.read_csv("/home/alon/menow_home_ass/PBTA_RNA/data_mrna_expression_continuous_rna_seq_v2_mrna.txt", sep="\t")

# Transpose expression: samples as rows, genes as columns
expr_t = expr.set_index('Hugo_Symbol').T.iloc[2:]  # drop Entrez_Gene_Id row
expr_t.index.name = 'SAMPLE_ID'
expr_t = expr_t.reset_index()

# Ensure SAMPLE_ID types match
expr_t['SAMPLE_ID'] = expr_t['SAMPLE_ID'].astype(str)
samples['SAMPLE_ID'] = samples['SAMPLE_ID'].astype(str)

# Merge: expression + clinical sample data → patient data
merged = samples.merge(expr_t, on='SAMPLE_ID', how='inner')  # or 'left'
merged = merged.merge(patients, on='PATIENT_ID', how='left')

print(f"Merged shape: {merged.shape}")
print(f"Clinical samples with expression: {merged['SAMPLE_ID'].nunique()}")
print(f"Unique patients: {merged['PATIENT_ID'].nunique()}")
```

## 6. Quick Facts Summary

- **37,680 genes**, all with unique Hugo_Symbol identifiers
- **2,516 tumor samples** with expression data (no normals)
- **All 2,516 map 1:1** to SAMPLE_ID in the clinical file
- **42% of clinical samples** lack expression data (WGS/targeted only)
- Each patient can have **1 or more** samples in the expression matrix

## 7. Quantification & Normalization

The exact quantification method for this expression matrix is **not documented** in the file headers or available metadata. The following is inferred purely from empirical analysis of the values.

### 7.1 Observed Properties

| Property | Value | Interpretation |
|---|---|---|
| Value range | 0 – ~695,000 | Linear scale (not log-transformed) |
| Mean expression | 25.75 | Dominated by lowly-expressed genes |
| Median expression | 0.06 | >50% of gene-sample pairs are near-zero |
| Zero fraction | 43.54% | Many genes not expressed in a given tissue |
| CV of sample sums | **2.98%** | Very low — strong between-sample normalization applied |
| Values are integers? | No | Not raw FASTQ counts |

### 7.2 What the Data Likely Represents

Given the continuous values (not integers), the range, and the very low CV of sample sums, these are **normalized expression estimates** — not raw counts.

Features consistent with **length-normalized quantification** (e.g., TPM, FPKM, or RSEM expected counts):
- Continuous values with many zeros
- Sample sums are nearly equal (CV ~3%), suggesting effective per-sample scaling
- No negative values

### 7.3 What We Don't Know

- ❌ Which quantifier was used (RSEM? Kallisto? Salmon? STAR/FeatureCounts?)
- ❌ What exact normalization pipeline was applied
- ❌ Whether these are TPM, FPKM, or another unit
- ❌ What reference genome / annotation was used

### 7.4 Practical Handling

```python
import numpy as np
import pandas as pd

expr = pd.read_csv("/home/alon/menow_home_ass/PBTA_RNA/data_mrna_expression_continuous_rna_seq_v2_mrna.txt", sep="\t")

# For clustering / correlation / PCA — log-transform to stabilize variance
log_expr = np.log2(expr.iloc[:, 2:] + 1)

# For differential expression — since we don't know the exact tool, prefer
# non-parametric methods (Mann-Whitney/Wilcoxon) over parametric ones (t-test)
# that assume a specific distribution

# For survival analysis (high/low groups) — median-split on log values
median_val = log_expr.median(axis=1)

# For comparing to external datasets — note that the unknown normalization
# may limit cross-study comparability
```

### 7.5 Recommended Approach

1. **Log-transform** (`log2(x + 1)`) for all exploratory work
2. **Non-parametric tests** for group comparisons (Mann-Whitney, Kruskal-Wallis)
3. **No cross-study comparisons** without confirming normalization compatibility
4. If the exact pipeline is needed, check the PBTA study publication or contact the data providers
