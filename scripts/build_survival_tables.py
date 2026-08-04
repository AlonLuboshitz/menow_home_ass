"""Build the per-cancer-group survival comparison table (OS and EFS merged).

Data source: the per-group KM log-rank results produced by
  notebooks/survival_analysis_executed.ipynb  (Phase 1, sections 1A-1D),
saved as notebooks/survival_analysis/survival_analysis_results.csv.

Table layout
  Rows   = the 22 cancer groups tested (n >= 20) in any per-group KM comparison.
  Cols   = each variable (AGE, TF, TP, SEX) for both the OS and EFS outcomes
           (columns AGE-OS | AGE-EFS | TF-OS | TF-EFS | TP-OS | TP-EFS | SEX-OS | SEX-EFS).
  Cell   = FDR_WithinGroup (BH-FDR computed within each cancer group),
           bolded when q < 0.05 (significant); "--" when the group was not
           tested for that variable (< 20 samples with complete data).
  Last row = N significant per variable (count of q < 0.05 cells).

PREDISPOSITION is intentionally omitted: the executed notebook has no
per-group PREDISPOSITION analysis (only global binary/per-type KM).

Caveat: FDR_WithinGroup is not reproducible from the executed notebook
(cell 46 recomputes a per-Phase FDR_BH); this column came from the later
"within-cancer-group FDR" rerun and is read verbatim.

Outputs
  build/survival_tables.tex                        -- LaTeX (single merged table)
  notebooks/survival_analysis/survival_tables_os.csv
  notebooks/survival_analysis/survival_tables_efs.csv
"""
import os

import pandas as pd

RESULTS = "/home/alon/menow_home_ass/notebooks/survival_analysis/survival_analysis_results.csv"
OUT_TEX = "/home/alon/menow_home_ass/build/survival_tables.tex"
OUT_OS = "/home/alon/menow_home_ass/notebooks/survival_analysis/survival_tables_os.csv"
OUT_EFS = "/home/alon/menow_home_ass/notebooks/survival_analysis/survival_tables_efs.csv"

VAR_COL = {"AGE": "AGE", "TUMOR_FRACTION": "TF",
           "TUMOR_PLOIDY": "TP", "SEX": "SEX"}

# Interleaved variable x outcome columns: variable with an OS/EFS suffix.
COLUMNS = [
    ("OS", "AGE"), ("EFS", "AGE"),
    ("OS", "TUMOR_FRACTION"), ("EFS", "TUMOR_FRACTION"),
    ("OS", "TUMOR_PLOIDY"), ("EFS", "TUMOR_PLOIDY"),
    ("OS", "SEX"), ("EFS", "SEX"),
]
HEADERS = [f"{VAR_COL[var]}-{outcome}" for outcome, var in COLUMNS]

COMPARISONS = {
    outcome: {var: f"{VAR_COL[var]} \u00d7 {outcome} (KM, per-group)"
              for var in ("AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY", "SEX")}
    for outcome in ("OS", "EFS")
}


def fmt_q(q):
    """Format an FDR q value (or '--' for untested)."""
    if pd.isna(q):
        return "--"
    if q < 0.0001:
        return "<0.0001"
    return f"{q:.4f}"


def build_matrix(outcome):
    per = df_all[df_all["Group"] != "global"]
    groups = sorted(per["Group"].unique())
    mat = pd.DataFrame(index=groups)
    for var in ("AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY", "SEX"):
        sub = per[per["Comparison"] == COMPARISONS[outcome][var]][["Group", "FDR_WithinGroup"]]
        mat[VAR_COL[var]] = sub.set_index("Group")["FDR_WithinGroup"]
    sig = mat < 0.05
    mat.loc["N significant", :] = sig.sum()
    return mat


def build_merged(mat_os, mat_efs):
    """Merge the OS and EFS matrices into one interleaved (variable x outcome) matrix."""
    merged = pd.DataFrame(index=mat_os.index)
    for outcome, var in COLUMNS:
        src = mat_os if outcome == "OS" else mat_efs
        merged[f"{VAR_COL[var]}-{outcome}"] = src[VAR_COL[var]]
    return merged


def to_latex(mat, headers, label, caption):
    nvar = len(headers)
    rows = []
    rows.append(r"\begin{table}[H]")
    rows.append(r"\centering")
    rows.append(f"\\caption{{{caption}}}\\label{{{label}}}")
    rows.append(r"\resizebox{\textwidth}{!}{%")
    rows.append(r"\begin{tabular}{l" + "c" * nvar + "}")
    rows.append(r"\toprule")
    rows.append("Cancer group & " + " & ".join(headers) + r" \\")
    rows.append(r"\midrule")
    body = mat.drop(index="N significant")
    for group in body.index:
        cells = []
        for col in headers:
            q = mat.loc[group, col]
            val = fmt_q(q)
            if q < 0.05:
                val = f"\\textbf{{{val}}}"
            cells.append(val)
        rows.append(" & ".join([group] + cells) + r" \\")
    rows.append(r"\midrule")
    n_row = mat.loc["N significant"]
    rows.append("N significant (q$<$0.05) & "
                + " & ".join(f"{int(v)}" for v in n_row.values) + r" \\")
    rows.append(r"\bottomrule")
    rows.append(r"\end{tabular}}")
    rows.append(r"\end{table}")
    return "\n".join(rows)


df_all = pd.read_csv(RESULTS)

mat_os = build_matrix("OS")
mat_efs = build_matrix("EFS")
merged = build_merged(mat_os, mat_efs)

for outcome, mat in (("OS", mat_os), ("EFS", mat_efs)):
    out = OUT_OS if outcome == "OS" else OUT_EFS
    os.makedirs(os.path.dirname(out), exist_ok=True)
    mat.round(4).to_csv(out)

os.makedirs(os.path.dirname(OUT_TEX), exist_ok=True)
caption = ("Per-cancer-group survival differences for the OS and EFS outcomes. "
           "Each variable is shown with an outcome suffix (AGE-OS, AGE-EFS, TF-OS, TF-EFS, "
           "TP-OS, TP-EFS, SEX-OS, SEX-EFS); cells show the BH-FDR q-value computed within "
           "each cancer group; bold = significant (q$<$0.05); `--' = not tested "
           "(group has $<$20 samples with complete data for that variable).")
with open(OUT_TEX, "w") as f:
    f.write(to_latex(merged, HEADERS, "Table x", caption))
    f.write("\n")

print("saved:", OUT_OS)
print("saved:", OUT_EFS)
print("saved:", OUT_TEX)

print("\n=== VERIFICATION ===")
print("\nMerged table:")
n_row = merged.loc["N significant"]
print(f"  rows (groups) = {len(merged) - 1}")
for outcome, var in COLUMNS:
    col = f"{VAR_COL[var]}-{outcome}"
    qs = merged.loc[merged.index[:-1], col].dropna()
    print(f"  {col:>6s}: tested={len(qs):2d}  significant={int(n_row[col]):2d}")
print("  last row (N significant):", " | ".join(f"{c}={int(n_row[c])}" for c in HEADERS))
