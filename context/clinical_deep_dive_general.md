# Clinical Deep Dive — General Methodology & Conventions

This document defines the **shared methodology, formatting conventions, imports, and phase structure** for all deep-dive clinical analyses of the PBTA_RNA dataset. It is the master reference for any agent building analysis notebooks in this series.

---

## 1. Relation to the Basic Notebook

The basic notebook (`clinical_analysis.ipynb`) provides exploratory overviews (distributions, missingness, basic KM curves). The deep-dive notebooks **extend** that work with:

- Deeper statistical comparisons (per-group, post-hoc, effect sizes, FDR)
- Both **binary outcome** (OS_STATUS/EFS_STATUS) and **time-to-event** (OS_MONTHS/EFS_MONTHS + event indicator) analyses
- Multivariate models (Cox PH)
- Unsupervised subgroup discovery (PCA, t-SNE, K-means, FAMD)

Each deep-dive section that overlaps with a basic notebook step will reference it: *"See basic notebook Step N for the exploratory version."*

---

## 2. Data Sources

All data from `/home/alon/menow_home_ass/PBTA_RNA/`:

| File | Description | Rows |
|------|-------------|------|
| `data_clinical_patient_attributes.txt` | Patient-level clinical data (13 columns) | ~2,871 |
| `data_clinical_sample_attributes.txt` | Sample-level clinical data (24 columns) | ~4,313 |

### File Format

Both files share a 4-line metadata header. Read with `header=4` (0-indexed, row 4 is the 5th row) and `sep='	'`.

### Cleaning Functions

Reuse the **exact same functions** from the basic notebook's `imports.py`:

```python
# Located at: notebooks/clinical_analysis/nb_cells/imports.py
from imports import (
    read_patients, read_samples,
    clean_os, clean_efs, clean_race_eth,
    clean_pred, clean_subtype, clean_tf_tp,
    kaplan_meier, add_km, logrank2, logrank_multi
)
```

### Missing Value Handling

| Column(s) | Strategy |
|-----------|----------|
| `AGE` | Leave as NaN; exclude from age-specific plots individually |
| `RACE` | Merge blank + `Not Reported` + `Reported Unknown` → `Unknown` |
| `ETHNICITY` | Merge blank + `Reported Unknown` → `Unknown` |
| `OS_MONTHS` / `OS_STATUS` | KM curves/Cox: drop rows without complete (time, status) pair; binary tests: use available labels |
| `EFS_MONTHS` / `EFS_STATUS` | Same as OS |
| `CANCER_PREDISPOSITIONS` | Merge blank + `Not Reported` → `Unknown`; `None documented` → `No predisposition` |
| `TUMOR_FRACTION` / `TUMOR_PLOIDY` | For boxplots: show `Unknown` as separate category; for correlations/PCA/t-SNE: drop missing and report dropout |
| `MOLECULAR_SUBTYPE` | Merge blank + `To be classified` → `Unclassified` |
| `SEX` | Keep blanks as NaN |
| `SPECIMEN_ID` | Split multi-value entries on `;` only if needed |

### Sample Size Thresholds

- **Statistical tests (Phases 1–4):** Only include groups with **n ≥ 20** samples. Print which groups pass/fail.
- **Unsupervised learning (Phase 5):** Only include groups with **n ≥ 50** samples.

---

## 3. Statistical Methodology

### 3.1 Normality Testing

For every continuous-variable–group combination:
- If n < 5000: **Shapiro-Wilk** test
- If n ≥ 5000: **Kolmogorov-Smirnov** test
- Non-normal if p < 0.05
- Expect most biomedical data to be non-normal → prefer non-parametric tests

### 3.2 Test Selection

| Scenario | Binary Outcome (STATUS) | Time-to-Event (MONTHS + event) |
|----------|------------------------|-------------------------------|
| 1 continuous × 2 groups | Mann-Whitney U | KM (dichotomized at median) + log-rank |
| 1 continuous × 3+ groups | Kruskal-Wallis + Dunn's post-hoc | KM (dichotomized) + pairwise log-rank |
| 2 categorical variables | Chi-squared (or Fisher's exact if expected < 5) | KM + log-rank |
| Correlation (2 continuous) | Spearman's ρ | — |
| Multivariate (Phase 4) | Logistic regression (if needed) | **Cox Proportional Hazards** |

Both approaches should be applied wherever possible — the binary test gives a simple answer ("do survivors differ from non-survivors?"), while the time-to-event test (KM + log-rank) answers the clinically relevant question: "does this factor predict survival timing?" The more powerful multivariate Cox PH analysis is reserved for Phase 4.

### 3.3 Effect Sizes

Report alongside every p-value:

| Test | Effect Size | Interpretations |
|------|-------------|----------------|
| Mann-Whitney U | **Cliff's delta** (d) | small ≥ 0.147, medium ≥ 0.33, large ≥ 0.474 |
| Kruskal-Wallis | **Epsilon-squared** (ε²) | small ≥ 0.01, medium ≥ 0.06, large ≥ 0.14 |
| Chi-squared | **Cramer's V** (φc) | small ≥ 0.1, medium ≥ 0.3, large ≥ 0.5 |
| Cox PH | **Hazard Ratio (HR)** with 95% CI | HR = 1: no effect, HR > 1: increased hazard, HR < 1: protective |
| Spearman | **ρ** | weak ≥ 0.1, moderate ≥ 0.3, strong ≥ 0.5 |

#### Cliff's Delta Formula
```python
def cliffs_delta(x, y):
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0: return 0
    more = sum(1 for xi in x for yi in y if xi > yi) / (n_x * n_y)
    less = sum(1 for xi in x for yi in y if xi < yi) / (n_x * n_y)
    return more - less
```

#### Epsilon-squared Formula
```python
def epsilon_sq(H, k, n):
    return (H - k + 1) / (n - k) if n > k else 0
```

#### Cramer's V Formula
```python
def cramers_v(contingency_table):
    chi2 = chi2_contingency(contingency_table)[0]
    n = contingency_table.sum().sum()
    phi2 = chi2 / n
    r, k = contingency_table.shape
    return np.sqrt(phi2 / min(k-1, r-1))
```

### 3.4 Post-hoc Pairwise Comparisons

After a significant Kruskal-Wallis (p < 0.05), use **Dunn's test** with Benjamini-Hochberg correction — NOT Mann-Whitney on all pairs (inflates Type I error).

```python
import scikit_posthocs as sp
dunn_results = sp.posthoc_dunn(
    data, val_col='AGE', group_col='CANCER_GROUP',
    p_adjust='fdr_bh'
)
```

After a significant global log-rank (3+ groups), use **pairwise log-rank tests** with BH correction.

### 3.5 FDR Correction

- **Within each family of related tests** (e.g., all per-cancer-group AGE × OS comparisons), apply **Benjamini-Hochberg** FDR correction.
- Display both raw p-value and FDR-adjusted q-value in the summary table.
- For post-hoc Dunn's: use `p_adjust='fdr_bh'`.
- For pairwise log-rank: collect p-values into array, apply `multipletests(pvals, method='fdr_bh')`.

```python
from statsmodels.stats.multitest import multipletests
rejected, q_values, _, _ = multipletests(p_values, method='fdr_bh')
```

### 3.6 Significance Notation

| p / q range | Symbol | Label |
|-------------|--------|-------|
| ≥ 0.05 | ❌ | Not significant |
| 0.01 – 0.05 | ✅ | p < 0.05 |
| 0.001 – 0.01 | ✅✅ | p < 0.01 |
| < 0.001 | ✅✅✅ | p < 0.001 |

---

## 4. Formatting & Coding Conventions

### 4.1 General

- Clean, readable Python with comments for non-obvious logic
- Pandas method chaining where appropriate
- All steps independently runnable (repeat imports if needed)
- File paths: absolute, rooted at `/home/alon/menow_home_ass/`

### 4.2 Plotting

- **Plotly** preferred for all plots (interactive: hover, zoom, toggle legend)
- All plots must have: **title, axis labels, legend** (if applicable)
- Statistical annotations: place p-value / effect size in the plot title or as an annotation box
- KM curves: show both OS and EFS (side-by-side or top-bottom), include risk tables, annotate with log-rank p-value
- Forest plots: show HR / effect size as dot, 95% CI as line, reference line at 1 (for HR) or 0 (for difference), color by significance
- Boxplots: include individual jittered points, show statistical test p-value in title
- Use consistent color palettes: `px.colors.qualitative.Plotly` or `Set1`/`Set2`

### 4.3 Validation Cells

Every analysis section must include a validation cell that prints:
- Sample sizes used
- Missing value counts
- Which groups passed/failed the n ≥ 20 threshold

```python
# Example validation
print(f"AGE × OS: {sum(df['OS_STATUS'].notna())} patients with OS data")
print(f"  LIVING: {(df['OS_STATUS'].str.contains('LIVING', na=False)).sum()}")
print(f"  DECEASED: {(df['OS_STATUS'].str.contains('DECEASED', na=False)).sum()}")
print(f"  Missing OS_STATUS: {df['OS_STATUS'].isna().sum()}")
```

### 4.4 Summary Table Design

A running results table is maintained throughout every phase notebook. Initialize as an empty list, append after each test, convert to DataFrame at the end.

**Columns:**

| Column | Type | Example |
|--------|------|---------|
| `Phase` | str | "1A" |
| `Comparison` | str | "AGE × OS (global)" |
| `Test` | str | "Mann-Whitney" |
| `Group` | str | "global" or "Low-grade glioma" |
| `N` | int | 2096 |
| `N_event` | int | 640 |
| `Statistic` | float | U=584032 |
| `p_value` | float | 0.002 |
| `FDR_BH` | float | 0.008 |
| `Significant` | str | "✅ p<0.01" |
| `Effect_Size` | str | "d=0.18 (small)" |
| `Basic_Ref` | str | "Step 14" or "—" |

At the end of the notebook:
- Sort by Phase then p_value
- Print summary with significance counts
- Save as CSV + Markdown

### 4.5 Step & Hypothesis Annotations

Every analytical step **must be preceded by a markdown cell** that clearly, in big words, states what the step is testing and what hypothesis is being checked. This is not optional — it ensures every test in the notebook is self-documenting and readable at a glance.

**Format (as a standalone markdown cell before the code cell):**

```markdown
## 📌 What We're Checking
**Does [variable] differ/associate/correlate with [outcome] in [population]?**
```

**Examples (as markdown cells):**

```markdown
## 📌 What We're Checking
**Is there a sex bias in this cancer group compared to 50:50?**
```

```markdown
## 📌 What We're Checking
**Does the predisposition makeup of this cancer group differ from all other cancers combined?**
```

```markdown
## 📌 What We're Checking
**Does AGE differ across cancer groups? Which groups are outliers?**
```

```markdown
## 📌 What We're Checking
**Is AGE associated with TUMOR_FRACTION?**
```

For non-statistical cells (descriptive plots, setup, summary), use a markdown cell with:

```markdown
## 📌 Purpose
**Describe the goal of this cell (no tests performed).**
```

The markdown cell must appear **immediately before** the corresponding code cell. The code cell itself should still include a `# ── Checks: ...` comment at the top for machine-readability and quick scanning.

---

## 5. Phases Overview

The full deep-dive analysis is organized into 6 phases. Each phase is designed to be **buildable as a standalone notebook** or concatenated into one large notebook.

| Phase | Title | Purpose | Threshold | Basic Ref |
|-------|-------|---------|-----------|-----------|
| **1** | Outcome Analysis (OS & EFS) | Test AGE, TF, TP, SEX, PREDISPOSITION against OS and EFS (binary + time-to-event), globally and per cancer group | n ≥ 20 | Steps 12–17 |
| **2** | Cross-Categorical Associations | Test associations between categorical clinical variables (SEX × CANCER_GROUP, CNS × TUMOR_TYPE, SUBTYPE × RACE, etc.) | n ≥ 20 | Steps 8, 15, 18 |
| **3** | Numeric Comparisons & Correlations | Test AGE × CANCER_GROUP (pairwise), AGE × SUBTYPE, TF × SUBTYPE, TP × SUBTYPE; Spearman correlations between all numeric pairs | n ≥ 20 | Step 14 |
| **4** | Multivariate Models | Cox PH (global + per-group) for OS and EFS with multiple covariates | n ≥ 20 per group for OS/EFS | Steps 12, 17 |
| **5** | Unsupervised Subgroup Discovery | PCA + t-SNE + K-means on numeric features, FAMD on mixed features, per cancer group; cluster survival validation | n ≥ 50 | — |
| **6** | Summary | Aggregate results table across all phases, significance report | — | Step 19–20 |

### How Phases Connect

- **Phase 1** identifies which clinical variables are associated with outcome (univariate). These findings feed into...
- **Phase 2 & 3** explore associations between variables themselves (not outcome). They inform which confounders to adjust for in...
- **Phase 4** multivariate models, which ask: what is the *independent* contribution of each variable?
- **Phase 5** takes a complementary approach: instead of testing predefined hypotheses, let the data reveal hidden structure, then validate against survival.
- **Phase 6** consolidates everything.

---

## 6. Core Imports

All deep-dive notebooks share these imports. Include them in a single cell at the top of the notebook:

```python
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import kruskal, mannwhitneyu, chi2_contingency, spearmanr
from scipy.stats import fisher_exact, ks_2samp, shapiro
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
import scikit_posthocs as sp
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.compose import ColumnTransformer
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = "/home/alon/menow_home_ass/PBTA_RNA"
PATIENT_FILE = f"{DATA_DIR}/data_clinical_patient_attributes.txt"
SAMPLE_FILE = f"{DATA_DIR}/data_clinical_sample_attributes.txt"
```

For Phase 5 only, additionally import:
```python
from prince import FAMD
```

---

## 7. Common Helper Functions

These should be defined once at the top of each notebook and reused throughout:

```python
def read_patients():
    return pd.read_csv(PATIENT_FILE, sep="	", header=4,
                       dtype={"AGE": float, "AGE_IN_DAYS": float,
                              "OS_MONTHS": float, "EFS_MONTHS": float})

def read_samples():
    return pd.read_csv(SAMPLE_FILE, sep="	", header=4)

def clean_os(df):
    df = df.copy()
    df["OS_STATUS"] = df["OS_STATUS"].str.strip()
    df["os_label"] = df["OS_STATUS"].str.replace(r"^\d+:", "", regex=True)
    df["os_event"] = df["OS_STATUS"].apply(
        lambda x: 1 if pd.notna(x) and x.startswith("1:") else (0 if pd.notna(x) and x.startswith("0:") else np.nan))
    return df

def clean_efs(df):
    df = df.copy()
    df["EFS_STATUS"] = df["EFS_STATUS"].str.strip()
    df["efs_detail"] = df["EFS_STATUS"].str.replace(r"^\d+:", "", regex=True)
    df["efs_event"] = df["EFS_STATUS"].apply(
        lambda x: 0 if pd.notna(x) and x == "0:No Event"
        else (1 if pd.notna(x) and x != "1:NA" else np.nan))
    return df

def clean_race_eth(df):
    df = df.copy()
    df["RACE"] = df["RACE"].fillna("Unknown").replace({"Not Reported":"Unknown","Reported Unknown":"Unknown"})
    df["ETHNICITY"] = df["ETHNICITY"].fillna("Unknown").replace({"Not Reported":"Unknown","Reported Unknown":"Unknown"})
    return df

def clean_pred(df):
    df = df.copy()
    df["CANCER_PREDISPOSITIONS"] = df["CANCER_PREDISPOSITIONS"].fillna("Unknown")
    df["CANCER_PREDISPOSITIONS"] = df["CANCER_PREDISPOSITIONS"].replace("Not Reported","Unknown")
    df["CANCER_PREDISPOSITIONS"] = df["CANCER_PREDISPOSITIONS"].replace("None documented","No predisposition")
    return df

def clean_subtype(df):
    df = df.copy()
    df["MOLECULAR_SUBTYPE"] = df["MOLECULAR_SUBTYPE"].fillna("Unclassified")
    df["MOLECULAR_SUBTYPE"] = df["MOLECULAR_SUBTYPE"].replace("To be classified","Unclassified")
    return df

def clean_tf_tp(df):
    df = df.copy()
    df["TF_group"] = np.where(df["TUMOR_FRACTION"].isna(), "Unknown", "Measured")
    df["TP_group"] = np.where(df["TUMOR_PLOIDY"].isna(), "Unknown", "Measured")
    return df

# --- Statistical Helpers ---

def cliffs_delta(x, y):
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0: return 0
    more = sum(1 for xi in x for yi in y if xi > yi) / (n_x * n_y)
    less = sum(1 for xi in x for yi in y if xi < yi) / (n_x * n_y)
    return more - less

def epsilon_sq(H, k, n):
    return (H - k + 1) / (n - k) if n > k else 0

def cramers_v(ct):
    chi2 = chi2_contingency(ct)[0]
    n = ct.sum().sum()
    phi2 = chi2 / n
    r, k = ct.shape
    return np.sqrt(phi2 / min(k - 1, r - 1))

# --- Kaplan-Meier (from basic notebook) ---

def kaplan_meier(times, events):
    d = pd.DataFrame({"t": times, "e": events}).dropna().sort_values("t")
    surv = 1.0
    n = len(d)
    res = []
    for t, grp in d.groupby("t", sort=False):
        ne = int(grp["e"].sum())
        if ne > 0:
            surv *= (1 - ne / n)
        res.append({"t": t, "s": surv, "n": n, "ne": ne})
        n -= len(grp)
    return pd.DataFrame(res)

def add_km(fig, km, label, color):
    fig.add_trace(go.Scatter(
        x=km["t"], y=km["s"], mode="lines", name=label,
        line=dict(color=color, width=2, shape="hv"),
        legendgroup=label,
        hovertemplate=f"Time: %{{x}}<br>Survival: %{{y:.3f}}<extra>{label}</extra>"))
    return fig

# --- Log-rank tests (from basic notebook) ---

def logrank2(t1, e1, t2, e2):
    from scipy.stats import chi2
    all_t = sorted(set(pd.concat([pd.Series(t1.dropna()), pd.Series(t2.dropna())]).dropna()))
    if len(all_t) < 2: return 1.0
    o1e = 0; v = 0
    d1 = pd.DataFrame({"t": t1, "e": e1}).dropna()
    d2 = pd.DataFrame({"t": t2, "e": e2}).dropna()
    for t in all_t:
        r1 = (d1["t"] >= t).sum(); r2 = (d2["t"] >= t).sum(); nr = r1 + r2
        if nr == 0: continue
        o1 = int(((d1["t"] == t) & (d1["e"] == 1)).sum())
        o2 = int(((d2["t"] == t) & (d2["e"] == 1)).sum())
        ot = o1 + o2
        if ot == 0: continue
        e1 = ot * r1 / nr; o1e += (o1 - e1)
        if nr > 1: v += ot * (r1 / nr) * (r2 / nr) * (nr - ot) / (nr - 1)
    if v <= 0: return 1.0
    return 1 - chi2.cdf(o1e ** 2 / v, 1)

def logrank_multi(groups):
    from scipy.stats import chi2
    import numpy as np
    ng = len(groups)
    if ng < 2: return 1.0
    all_t = sorted(set(pd.concat([pd.Series(g[0].dropna()) for g in groups]).dropna()))
    if len(all_t) < 2: return 1.0
    O = np.zeros(ng); E = np.zeros(ng); V = np.zeros((ng, ng))
    for t in all_t:
        ar = np.array([(g[0] >= t).sum() for g in groups]); nr = ar.sum()
        if nr == 0: continue
        ev = np.array([((g[0] == t) & (g[1] == 1)).sum() for g in groups]); ot = ev.sum()
        if ot == 0: continue
        O += ev; E += ot * ar / nr
        if nr > 1:
            for i in range(ng):
                for j in range(ng):
                    if i == j: V[i, j] += ot * ar[i] / nr * (1 - ar[i] / nr) * (nr - ot) / (nr - 1)
                    else: V[i, j] -= ot * ar[i] / nr * ar[j] / nr * (nr - ot) / (nr - 1)
    try: return 1 - chi2.cdf((O - E) @ np.linalg.pinv(V) @ (O - E), ng - 1)
    except: return 1.0
```

---

## 8. Shared Quality Checklist

Every deep-dive notebook must pass these checks:

### Structure
- [ ] Summary table initialized at top (empty), updated throughout
- [ ] Each analysis independently runnable
- [ ] Validation cells present for every analysis

### Statistical
- [ ] Both binary (STATUS) AND time-to-event (MONTHS + event) analyses applied where applicable
- [ ] Effect size reported alongside p-value
- [ ] FDR correction applied within each family of tests
- [ ] Both raw p and FDR-adjusted q shown in summary table
- [ ] Per-group analyses respect thresholds (≥20 for tests, ≥50 for unsupervised)
- [ ] Dunn's post-hoc (not Mann-Whitney) after significant Kruskal-Wallis
- [ ] Normality test (Shapiro-Wilk / KS) run before parametric tests
- [ ] Cox PH proportional hazards assumption checked

### Visualization
- [ ] All plots: title, axis labels, legend, statistical annotation
- [ ] Plotly used for interactive plots
- [ ] KM curves show risk tables + log-rank p-value
- [ ] Forest plots show effect size + 95% CI

### Data Integrity
- [ ] Cleaning functions from basic notebook reused (not re-implemented)
- [ ] Missing value counts printed per analysis
- [ ] Merge validation (orphan check) at load time
- [ ] Per-group analyses print sample sizes and threshold status

### Summary
- [ ] Summary table saved as CSV and Markdown
- [ ] High-level counts printed (total tests, significant fraction, significant-after-FDR fraction)
- [ ] Basic notebook steps referenced where applicable
