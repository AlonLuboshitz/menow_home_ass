#!/usr/bin/env python3
"""Build the Phase 2+3 clinical associations notebook programmatically."""

import sys, os, nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT_DIR = '/home/alon/menow_home_ass/notebooks/clinical_associations'
NB_PATH = os.path.join(OUT_DIR, 'clinical_associations_analysis.ipynb')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'src'), exist_ok=True)

nb = new_notebook()
nb.metadata = {
    'kernelspec': {
        'display_name': 'Python 3',
        'language': 'python',
        'name': 'python3'
    },
    'language_info': {'name': 'python', 'version': '3.11.0'}
}

cells = []

# ── Cell 0: Title ──
cells.append(new_markdown_cell(
    '# Phase 2+3: Clinical Associations — Cross-Categorical & Numeric Comparisons'
))

# ── Cell 1: Summary ──
cells.append(new_markdown_cell(
    '## Summary\n\n'
    'This notebook implements **Phase 2 (Cross-Categorical Associations)** and '
    '**Phase 3 (Numeric Comparisons)** of the clinical deep-dive analysis.\n\n'
    '**Phase 2** tests for per-group enrichment of categorical variables across '
    'cancer groups using binomial tests, Fisher exact tests, and chi-squared tests. '
    'Each cancer group gets its own p-value for every comparison. Results are '
    'visualized with horizontal bar charts and annotated heatmaps.\n\n'
    '**Phase 3** compares numeric variables (AGE, TUMOR_FRACTION, TUMOR_PLOIDY) '
    'across cancer groups using Kruskal-Wallis and per-group Mann-Whitney tests, '
    'and computes Spearman correlations between numeric pairs with LOESS-smoothed '
    'scatter plots.\n\n'
    'All results are saved in a unified CSV format matching `survival_analysis_results.csv`.'
))

# ── Cell 2: Imports ──
cells.append(new_code_cell("""import sys, os, warnings, math
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, binom_test, kruskal, spearmanr, mannwhitneyu
import statsmodels.stats.multitest as smm
from itertools import combinations
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly_white"
from IPython.display import display

sys.path.insert(0, "/home/alon/menow_home_ass/notebooks/clinical_analysis/src")
from imports import (
    read_patients, read_samples,
    clean_os, clean_efs, clean_race_eth,
    clean_pred, clean_subtype, clean_tf_tp
)

print("✅ All imports loaded")"""))

# ── Cell 3: Data loading ──
cells.append(new_code_cell("""pat = read_patients()
smp = read_samples()
df = pat.merge(smp, on="PATIENT_ID", how="inner")
print(f"Merged dataset: {df.shape[0]} rows, {df.shape[1]} columns")

df = clean_race_eth(df)
df = clean_pred(df)
df = clean_subtype(df)
df = clean_tf_tp(df)
df = clean_os(df)
df = clean_efs(df)

for col in ["SEX", "CANCER_GROUP", "CANCER_PREDISPOSITIONS", "RACE", "MOLECULAR_SUBTYPE"]:
    if col in df.columns:
        df[col] = df[col].astype(str)

print("\\nSample sizes per cancer group:")
print(df["CANCER_GROUP"].value_counts().to_string())
print(f"\\nTotal unique patients: {df['PATIENT_ID'].nunique()}")
print(f"Total samples: {len(df)}")"""))

# ── Cell 4: Phase 2 header ──
cells.append(new_markdown_cell(
    '## Phase 2: Cross-Categorical Associations\n\n'
    'Each test produces **per-group p-values** — one row per cancer group per comparison. '
    'Tests include binomial tests (Test 1), Fisher exact enrichment (Tests 2, 4, 5), '
    'and descriptive heatmaps (Test 3). FDR is applied within each Comparison family.'
))

# ── Cell 5: Phase 2 helpers ──
cells.append(new_code_cell("""def sig_star(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

def run_enrichment(df, value_col, group_col="CANCER_GROUP", min_n=20):
    results = []
    values = [v for v in df[value_col].unique() if pd.notna(v) and str(v).strip() not in ("", "nan", "Unknown")]
    total_n = len(df)
    for val in sorted(values):
        val_mask = df[value_col] == val
        n_val_total = val_mask.sum()
        for group in sorted(df[group_col].unique()):
            g = df[df[group_col] == group]
            n_group = len(g)
            if n_group < min_n:
                continue
            in_group_and_val = val_mask[df[group_col] == group].sum()
            in_group_not_val = n_group - in_group_and_val
            not_in_group_and_val = n_val_total - in_group_and_val
            not_in_group_not_val = total_n - n_group - n_val_total + in_group_and_val
            table = [[in_group_and_val, in_group_not_val],
                     [not_in_group_and_val, not_in_group_not_val]]
            try:
                odds, p = fisher_exact(table)
            except Exception:
                odds, p = 1.0, 1.0
            results.append({
                "value": val, "group": group,
                "n_group": n_group,
                "n_with_val_in_group": in_group_and_val,
                "n_with_val_total": n_val_total,
                "odds_ratio": odds, "p_value": p
            })
    return results

results_ph2 = []
print("✅ Phase 2 helpers defined")"""))

# ── Cell 6: Test 1 — SEX enrichment per group vs 50:50 ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 1: SEX enrichment per group (binomial vs 50:50)")
print("=" * 70)

sex_sub = df[df["SEX"].isin(["Female", "Male"])].copy()
overall_counts = sex_sub["SEX"].value_counts()
print(f"Overall sex distribution: {overall_counts.to_dict()}")
print(f"Overall %Female: {overall_counts.get('Female', 0) / len(sex_sub) * 100:.1f}%")

for group in sorted(sex_sub["CANCER_GROUP"].unique()):
    g = sex_sub[sex_sub["CANCER_GROUP"] == group]
    n = len(g)
    if n < 20:
        continue
    n_f = (g["SEX"] == "Female").sum()
    n_m = (g["SEX"] == "Male").sum()
    pct_f = n_f / n * 100
    pct_m = n_m / n * 100
    p = binom_test(n_f, n, 0.5, alternative="two-sided")
    enriched = "Female" if n_f > n_m else "Male"
    effect = (n_f - n_m) / n
    results_ph2.append({
        "FDR_Family": "SEX enrichment",
        "Phase": "Phase 2",
        "Comparison": "SEX enrichment",
        "Test": "Binomial",
        "Group": group,
        "N": n,
        "N_events": n_f,
        "Statistic": f"{pct_f:.1f}%F",
        "p_value": round(p, 6),
        "Effect_Size": round(effect, 4),
        "Significant": sig_star(p)
    })
    print(f"  {group:40s} N={n:4d}  M={n_m:4d}  F={n_f:4d}  %M={pct_m:5.1f}  %F={pct_f:5.1f}  p={p:.4f}  {sig_star(p)}  enriched={enriched}")

df_test1 = pd.DataFrame([r for r in results_ph2 if r["Comparison"] == "SEX enrichment"])
if len(df_test1) > 0:
    df_test1 = df_test1.sort_values("Effect_Size")
    fig = go.Figure()
    pct_m_vals = [100 - float(r["Statistic"].replace("%F", "")) for _, r in df_test1.iterrows()]
    pct_f_vals = [float(r["Statistic"].replace("%F", "")) for _, r in df_test1.iterrows()]
    fig.add_trace(go.Bar(
        y=df_test1["Group"], x=pct_m_vals,
        name="%Male", orientation="h",
        marker_color="steelblue",
        text=[f"{v:.1f}%" for v in pct_m_vals],
        textposition="inside",
    ))
    fig.add_trace(go.Bar(
        y=df_test1["Group"], x=pct_f_vals,
        name="%Female", orientation="h",
        marker_color="lightcoral",
        text=[f"{v:.1f}%" for v in pct_f_vals],
        textposition="inside",
    ))
    sig_map = {"***": "darkred", "**": "red", "*": "lightcoral", "ns": "grey"}
    star_colors = [sig_map.get(s, "grey") for s in df_test1["Significant"]]
    fig.add_trace(go.Scatter(
        y=df_test1["Group"],
        x=[105] * len(df_test1),
        mode="text",
        text=df_test1["Significant"],
        textfont=dict(color=star_colors, size=14),
        name="significance",
        showlegend=False
    ))
    fig.update_layout(
        title="SEX composition per cancer group (binomial vs 50:50)",
        xaxis_title="Percent",
        barmode="overlay",
        height=max(400, 25 * len(df_test1)),
        width=900,
        bargap=0.15,
        xaxis=dict(range=[0, 115])
    )
    fig.add_vline(x=50, line_dash="dash", line_color="black", opacity=0.5)
    try:
        fig.show()
    except Exception:
        pass
    print(f"\\nTest 1 complete: {len(df_test1)} groups tested")"""))

# ── Cell 7: Test 2 — Predisposition enrichment per group (Fisher exact) ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 2: Predisposition enrichment per group (Fisher exact)")
print("=" * 70)

pred_counts = df["CANCER_PREDISPOSITIONS"].value_counts()
print(f"\\nPredisposition categories found ({len(pred_counts)}):")
for v, c in pred_counts.items():
    print(f"  {str(v)[:50]:50s} N={c}")

pred_results = run_enrichment(df, "CANCER_PREDISPOSITIONS")
print(f"\\nTotal enrichment tests: {len(pred_results)}")

for r in pred_results:
    results_ph2.append({
        "FDR_Family": "Predisposition enrichment",
        "Phase": "Phase 2",
        "Comparison": "Predisposition enrichment",
        "Test": "Fisher exact",
        "Group": str(r["value"]) + " -> " + str(r["group"]),
        "N": r["n_group"],
        "N_events": r["n_with_val_in_group"],
        "Statistic": round(r["odds_ratio"], 4),
        "p_value": round(r["p_value"], 6),
        "Effect_Size": round(r["odds_ratio"], 4),
        "Significant": sig_star(r["p_value"])
    })
    if r["p_value"] < 0.05:
        print(f"  {str(r['value'])[:40]:40s} x {str(r['group'])[:35]:35s} "
              f"N={r['n_group']:4d}  n_pred={r['n_with_val_in_group']:3d}  "
              f"OR={r['odds_ratio']:.2f}  p={r['p_value']:.4f}  {sig_star(r['p_value'])}")

# Heatmap: predisposition x cancer group with counts
pred_cross = pd.crosstab(df["CANCER_PREDISPOSITIONS"], df["CANCER_GROUP"])
n20 = df["CANCER_GROUP"].value_counts()
n20 = n20[n20 >= 20].index
pred_cross = pred_cross[[c for c in pred_cross.columns if c in n20]]

fig = px.imshow(
    pred_cross.values,
    x=pred_cross.columns.tolist(),
    y=pred_cross.index.tolist(),
    text_auto=True,
    aspect="auto",
    color_continuous_scale="YlOrRd",
    title="Predisposition x Cancer Group: Counts",
    height=max(400, 20 * len(pred_cross)),
    width=max(600, 100 * len(pred_cross.columns))
)
fig.update_layout(xaxis_tickangle=45)
try:
    fig.show()
except Exception:
    pass"""))

# ── Cell 8: Test 3 — SUBTYPE descriptive heatmap ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 3: SUBTYPE descriptive heatmap (no p-values)")
print("=" * 70)

subtype_cross = pd.crosstab(df["MOLECULAR_SUBTYPE"], df["CANCER_GROUP"])
n20 = df["CANCER_GROUP"].value_counts()
n20 = n20[n20 >= 20].index
subtype_cross = subtype_cross[[c for c in subtype_cross.columns if c in n20]]

subtype_totals = subtype_cross.sum(axis=1).sort_values(ascending=False)
if len(subtype_totals) > 30:
    top30 = subtype_totals.head(30).index
    subtype_cross = subtype_cross.loc[top30]
    print(f"Showing top 30 subtypes (out of {len(subtype_totals)} total)")

fig = px.imshow(
    subtype_cross.values,
    x=subtype_cross.columns.tolist(),
    y=subtype_cross.index.tolist(),
    text_auto=True,
    aspect="auto",
    color_continuous_scale="YlGnBu",
    title="Molecular Subtype x Cancer Group: Counts (descriptive, no tests)",
    height=max(500, 25 * len(subtype_cross)),
    width=max(800, 100 * len(subtype_cross.columns))
)
fig.update_layout(xaxis_tickangle=45)
try:
    fig.show()
except Exception:
    pass
print(f"Heatmap: {len(subtype_cross)} subtypes x {len(subtype_cross.columns)} groups")"""))

# ── Cell 9: Test 4 — SEX x PRED 3-way per group ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 4: SEX x Predisposition per group")
print("=" * 70)

def has_pred(val):
    if pd.isna(val) or str(val).strip() in ("", "nan", "Unknown", "No predisposition"):
        return "No predisposition"
    return "Any predisposition"

df["PRED_BINARY"] = df["CANCER_PREDISPOSITIONS"].apply(has_pred)
sub = df[df["SEX"].isin(["Female", "Male"])].copy()
sub = sub[sub["PRED_BINARY"].isin(["Any predisposition", "No predisposition"])].copy()

sex_pred_data = []
for cg in sorted(sub["CANCER_GROUP"].unique()):
    g = sub[sub["CANCER_GROUP"] == cg]
    n = len(g)
    if n < 20:
        continue
    ctab = pd.crosstab(g["SEX"], g["PRED_BINARY"])
    if ctab.shape != (2, 2):
        continue
    try:
        if ctab.min().min() < 5:
            odds_ratio, p = fisher_exact(ctab)
            test_name = "Fisher"
            eff = odds_ratio
        else:
            chi2, p, dof, expected = chi2_contingency(ctab, correction=True)
            test_name = "Chi2"
            n_total = ctab.values.sum()
            cramer_v = math.sqrt(chi2 / (n_total * (min(ctab.shape) - 1))) if n_total > 0 else 0
            eff = cramer_v
    except Exception:
        p, eff = 1.0, 0.0
        test_name = "Fisher"

    pct_f_pred = (ctab.loc["Female", "Any predisposition"] / ctab["Any predisposition"].sum() * 100
                  if "Any predisposition" in ctab.columns and ctab["Any predisposition"].sum() > 0 else 0)
    pct_f_nopred = (ctab.loc["Female", "No predisposition"] / ctab["No predisposition"].sum() * 100
                    if "No predisposition" in ctab.columns and ctab["No predisposition"].sum() > 0 else 0)

    results_ph2.append({
        "FDR_Family": "SEX x Predisposition",
        "Phase": "Phase 2",
        "Comparison": "SEX x Predisposition",
        "Test": test_name,
        "Group": cg,
        "N": n,
        "N_events": n,
        "Statistic": round(eff, 4),
        "p_value": round(p, 6),
        "Effect_Size": round(eff, 4),
        "Significant": sig_star(p)
    })
    sex_pred_data.append({
        "Group": cg, "pct_f_pred": pct_f_pred,
        "pct_f_nopred": pct_f_nopred, "p": p
    })
    print(f"  {cg:40s} N={n:4d}  %F(Pred)={pct_f_pred:5.1f}  %F(NoPred)={pct_f_nopred:5.1f}  p={p:.4f}  {sig_star(p)}")

# Side-by-side %Female bars
if sex_pred_data:
    sp_df = pd.DataFrame(sex_pred_data).sort_values("p")
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=sp_df["Group"], y=sp_df["pct_f_pred"],
        name="%Female (Any predisposition)",
        marker_color="indianred",
        text=sp_df["pct_f_pred"].apply(lambda v: f"{v:.1f}%"),
        textposition="outside"
    ))
    fig.add_trace(go.Bar(
        x=sp_df["Group"], y=sp_df["pct_f_nopred"],
        name="%Female (No predisposition)",
        marker_color="steelblue",
        text=sp_df["pct_f_nopred"].apply(lambda v: f"{v:.1f}%"),
        textposition="outside"
    ))
    fig.update_layout(
        title="%Female by Predisposition status per cancer group",
        xaxis_tickangle=45,
        barmode="group",
        height=500, width=1000,
        yaxis_title="%Female"
    )
    try:
        fig.show()
    except Exception:
        pass"""))

# ── Cell 10: Test 5 — Race enrichment per group (Fisher exact) ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 5: Race enrichment per group (Fisher exact)")
print("=" * 70)

race_sub = df[~df["RACE"].isin(["Unknown", "nan"])].copy()
race_counts = race_sub["RACE"].value_counts()
print(f"\\nRace categories (excluding Unknown):")
for v, c in race_counts.items():
    print(f"  {v:40s} N={c}")

race_results = run_enrichment(race_sub, "RACE")
print(f"\\nTotal enrichment tests: {len(race_results)}")

for r in race_results:
    results_ph2.append({
        "FDR_Family": "Race enrichment",
        "Phase": "Phase 2",
        "Comparison": "Race enrichment",
        "Test": "Fisher exact",
        "Group": str(r["value"]) + " -> " + str(r["group"]),
        "N": r["n_group"],
        "N_events": r["n_with_val_in_group"],
        "Statistic": round(r["odds_ratio"], 4),
        "p_value": round(r["p_value"], 6),
        "Effect_Size": round(r["odds_ratio"], 4),
        "Significant": sig_star(r["p_value"])
    })
    if r["p_value"] < 0.05:
        print(f"  {str(r['value'])[:40]:40s} x {str(r['group'])[:35]:35s} "
              f"N={r['n_group']:4d}  n_race={r['n_with_val_in_group']:3d}  "
              f"OR={r['odds_ratio']:.2f}  p={r['p_value']:.4f}  {sig_star(r['p_value'])}")

# Heatmap: race x cancer group with counts
race_cross = pd.crosstab(race_sub["RACE"], race_sub["CANCER_GROUP"])
n20 = df["CANCER_GROUP"].value_counts()
n20 = n20[n20 >= 20].index
race_cross = race_cross[[c for c in race_cross.columns if c in n20]]

fig = px.imshow(
    race_cross.values,
    x=race_cross.columns.tolist(),
    y=race_cross.index.tolist(),
    text_auto=True,
    aspect="auto",
    color_continuous_scale="PuBuGn",
    title="Race x Cancer Group: Counts",
    height=max(400, 25 * len(race_cross)),
    width=max(600, 100 * len(race_cross.columns))
)
fig.update_layout(xaxis_tickangle=45)
try:
    fig.show()
except Exception:
    pass"""))

# ── Cell 11: Phase 2 FDR + display ──
cells.append(new_code_cell("""df_ph2 = pd.DataFrame(results_ph2)

if len(df_ph2) > 0:
    df_ph2["FDR_WithinFamily"] = np.nan
    families = df_ph2["FDR_Family"].unique() if "FDR_Family" in df_ph2.columns else ["All"]
    for fam in families:
        mask = df_ph2["FDR_Family"] == fam
        pvals = df_ph2.loc[mask, "p_value"].values
        if len(pvals) > 1:
            reject, qvals, _, _ = smm.multipletests(pvals, method="fdr_bh")
            df_ph2.loc[mask, "FDR_WithinFamily"] = qvals
        elif len(pvals) == 1:
            df_ph2.loc[mask, "FDR_WithinFamily"] = pvals[0]

    print(f"\\nPhase 2: {len(df_ph2)} total tests")
    print(f"Families: {df_ph2['FDR_Family'].value_counts().to_dict()}")
    sig_raw = (df_ph2["p_value"] < 0.05).sum()
    sig_fdr = (df_ph2["FDR_WithinFamily"] < 0.05).sum()
    print(f"Significant (p<0.05): {sig_raw}")
    print(f"FDR-significant (q<0.05): {sig_fdr}")

    display_cols = ["Phase", "Comparison", "Test", "Group", "N", "N_events",
                    "Statistic", "p_value", "FDR_WithinFamily", "Significant", "Effect_Size"]
    avail = [c for c in display_cols if c in df_ph2.columns]
    display(df_ph2[avail].style.format({
        "p_value": "{:.6f}",
        "FDR_WithinFamily": "{:.4f}",
        "Effect_Size": "{:.4f}"
    }))
else:
    print("No Phase 2 results.")"""))

# ── Cell 12: Phase 3 header ──
cells.append(new_markdown_cell(
    '## Phase 3: Numeric Comparisons\n\n'
    'Tests: Kruskal-Wallis (global across all cancer groups), per-group '
    'Mann-Whitney (each group vs all others), and Spearman correlations '
    'between numeric variables. FDR is applied within each Comparison family.'
))

# ── Cell 13: Phase 3 helpers ──
cells.append(new_code_cell("""def kw_per_group(df, num_col, cat_col="CANCER_GROUP", min_n=20):
    results = []
    sub = df[[num_col, cat_col]].dropna()
    sub = sub[sub[cat_col].isin(
        sub[cat_col].value_counts()[sub[cat_col].value_counts() >= min_n].index)]

    # Global KW
    groups_data = [sub[sub[cat_col] == g][num_col].values for g in sub[cat_col].unique()]
    if len(groups_data) >= 2:
        h, p = kruskal(*groups_data)
        n_total = len(sub)
        k = len(groups_data)
        eps_sq = (h - k + 1) / (n_total - k) if n_total > k else 0
        results.append({
            "Comparison": f"{num_col} x {cat_col}",
            "Test": "Kruskal-Wallis", "Group": "global",
            "N": n_total, "N_events": n_total,
            "Statistic": f"H={h:.4f}",
            "p_value": p, "Effect_Size": round(eps_sq, 4)
        })

    # Per-group MW
    for g_name in sub[cat_col].unique():
        g_vals = sub[sub[cat_col] == g_name][num_col].values
        other_vals = sub[sub[cat_col] != g_name][num_col].values
        try:
            u, p = mannwhitneyu(g_vals, other_vals, alternative="two-sided")
        except ValueError:
            continue
        n1, n2 = len(g_vals), len(other_vals)
        z = (u - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12) if n1 * n2 > 0 else 0
        r = z / math.sqrt(n1 + n2) if (n1 + n2) > 0 else 0
        results.append({
            "Comparison": f"{num_col} x {cat_col}",
            "Test": "Mann-Whitney", "Group": g_name,
            "N": len(sub), "N_events": len(g_vals),
            "Statistic": f"U={u:.1f}",
            "p_value": p, "Effect_Size": round(r, 4)
        })
    return results

def boxplot_numeric(df, num_col, cat_col, title, ylabel, min_n=20):
    sub = df[[num_col, cat_col]].dropna()
    sub = sub[sub[cat_col].isin(
        sub[cat_col].value_counts()[sub[cat_col].value_counts() >= min_n].index)]
    medians = sub.groupby(cat_col)[num_col].median().sort_values()
    sub[cat_col] = pd.Categorical(sub[cat_col], categories=medians.index, ordered=True)
    fig = px.box(sub, x=cat_col, y=num_col, title=title,
                 points="outliers", color=cat_col,
                 color_discrete_sequence=px.colors.qualitative.Set3)
    fig.update_layout(xaxis_tickangle=45, height=500, width=1000,
                      yaxis_title=ylabel, showlegend=False)
    try:
        fig.show()
    except Exception:
        pass

results_ph3 = []
print("✅ Phase 3 helpers defined")"""))

# ── Cell 14: Test 6 — AGE x CANCER_GROUP ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 6: AGE x CANCER_GROUP")
print("=" * 70)

age_results = kw_per_group(df, "AGE", "CANCER_GROUP")
for r in age_results:
    results_ph3.append({
        "FDR_Family": "AGE x CG",
        "Phase": "Phase 3",
        "Comparison": r["Comparison"],
        "Test": r["Test"],
        "Group": r["Group"],
        "N": r["N"],
        "N_events": r["N_events"],
        "Statistic": r["Statistic"],
        "p_value": round(r["p_value"], 6),
        "Effect_Size": r["Effect_Size"],
        "Significant": sig_star(r["p_value"])
    })
    if r["Test"] == "Kruskal-Wallis":
        print(f"  KW: {r['Statistic']}, p={r['p_value']:.6f}, eps2={r['Effect_Size']:.4f}")

sig_mw = [r for r in age_results if r["Test"] == "Mann-Whitney" and r["p_value"] < 0.05]
n_mw = len([r for r in age_results if r["Test"] == "Mann-Whitney"])
print(f"  Mann-Whitney: {n_mw} tests, {len(sig_mw)} significant (p<0.05)")

n_age = df["AGE"].notna().sum()
boxplot_numeric(df, "AGE", "CANCER_GROUP",
                f"AGE distribution by CANCER_GROUP (N={n_age})",
                "AGE (years)")"""))

# ── Cell 15: Test 7 — TF x CANCER_GROUP ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 7: TUMOR_FRACTION x CANCER_GROUP")
print("=" * 70)

tf_results = kw_per_group(df, "TUMOR_FRACTION", "CANCER_GROUP")
for r in tf_results:
    results_ph3.append({
        "FDR_Family": "TF x CG",
        "Phase": "Phase 3",
        "Comparison": r["Comparison"],
        "Test": r["Test"],
        "Group": r["Group"],
        "N": r["N"],
        "N_events": r["N_events"],
        "Statistic": r["Statistic"],
        "p_value": round(r["p_value"], 6),
        "Effect_Size": r["Effect_Size"],
        "Significant": sig_star(r["p_value"])
    })
    if r["Test"] == "Kruskal-Wallis":
        print(f"  KW: {r['Statistic']}, p={r['p_value']:.6f}, eps2={r['Effect_Size']:.4f}")

sig_mw = [r for r in tf_results if r["Test"] == "Mann-Whitney" and r["p_value"] < 0.05]
n_mw = len([r for r in tf_results if r["Test"] == "Mann-Whitney"])
print(f"  Mann-Whitney: {n_mw} tests, {len(sig_mw)} significant (p<0.05)")

n_tf = df["TUMOR_FRACTION"].notna().sum()
boxplot_numeric(df, "TUMOR_FRACTION", "CANCER_GROUP",
                f"TF distribution by CANCER_GROUP (N={n_tf})",
                "TUMOR_FRACTION")"""))

# ── Cell 16: Test 8 — TP x CANCER_GROUP ──
cells.append(new_code_cell("""print("=" * 70)
print("Test 8: TUMOR_PLOIDY x CANCER_GROUP")
print("=" * 70)

tp_results = kw_per_group(df, "TUMOR_PLOIDY", "CANCER_GROUP")
for r in tp_results:
    results_ph3.append({
        "FDR_Family": "TP x CG",
        "Phase": "Phase 3",
        "Comparison": r["Comparison"],
        "Test": r["Test"],
        "Group": r["Group"],
        "N": r["N"],
        "N_events": r["N_events"],
        "Statistic": r["Statistic"],
        "p_value": round(r["p_value"], 6),
        "Effect_Size": r["Effect_Size"],
        "Significant": sig_star(r["p_value"])
    })
    if r["Test"] == "Kruskal-Wallis":
        print(f"  KW: {r['Statistic']}, p={r['p_value']:.6f}, eps2={r['Effect_Size']:.4f}")

sig_mw = [r for r in tp_results if r["Test"] == "Mann-Whitney" and r["p_value"] < 0.05]
n_mw = len([r for r in tp_results if r["Test"] == "Mann-Whitney"])
print(f"  Mann-Whitney: {n_mw} tests, {len(sig_mw)} significant (p<0.05)")

n_tp = df["TUMOR_PLOIDY"].notna().sum()
boxplot_numeric(df, "TUMOR_PLOIDY", "CANCER_GROUP",
                f"TP distribution by CANCER_GROUP (N={n_tp})",
                "TUMOR_PLOIDY")"""))

# ── Cell 17: Tests 9-11 — Correlations ──
cells.append(new_code_cell("""print("=" * 70)
print("Tests 9-11: Numeric correlations (Spearman + LOESS scatter)")
print("=" * 70)

def spearman_with_plot(df, col1, col2, label1, label2, min_n=20):
    sub = df[[col1, col2]].dropna()
    if len(sub) < min_n:
        print(f"  Skipping {label1}x{label2}: N={len(sub)} < {min_n}")
        return None
    rho, p = spearmanr(sub[col1], sub[col2])
    fig = px.scatter(
        sub, x=col1, y=col2,
        title=f"{label1} vs {label2} (rho={rho:.4f}, p={p:.6f}, N={len(sub)})",
        opacity=0.5,
        trendline="lowess",
        trendline_color_override="red"
    )
    fig.update_layout(height=500, width=600)
    try:
        fig.show()
    except Exception:
        pass
    print(f"  {label1} x {label2}: rho={rho:.4f}, p={p:.6f}, N={len(sub)}")
    return {
        "FDR_Family": "Numeric correlations",
        "Phase": "Phase 3",
        "Comparison": "Numeric correlations",
        "Test": "Spearman",
        "Group": "global",
        "N": len(sub),
        "N_events": len(sub),
        "Statistic": round(rho, 4),
        "p_value": round(p, 6),
        "Effect_Size": round(rho, 4),
        "Significant": sig_star(p)
    }

for c1, c2, l1, l2 in [
    ("AGE", "TUMOR_FRACTION", "AGE", "TF"),
    ("AGE", "TUMOR_PLOIDY", "AGE", "TP"),
    ("TUMOR_FRACTION", "TUMOR_PLOIDY", "TF", "TP"),
]:
    res = spearman_with_plot(df, c1, c2, l1, l2)
    if res:
        results_ph3.append(res)"""))

# ── Cell 18: Phase 3 FDR + display ──
cells.append(new_code_cell("""if len(results_ph3) > 0:
    df_ph3 = pd.DataFrame(results_ph3)
    df_ph3["FDR_WithinFamily"] = np.nan
    families = df_ph3["FDR_Family"].unique()
    for fam in families:
        mask = df_ph3["FDR_Family"] == fam
        pvals = df_ph3.loc[mask, "p_value"].values
        if len(pvals) > 1:
            reject, qvals, _, _ = smm.multipletests(pvals, method="fdr_bh")
            df_ph3.loc[mask, "FDR_WithinFamily"] = qvals
        elif len(pvals) == 1:
            df_ph3.loc[mask, "FDR_WithinFamily"] = pvals[0]

    print(f"\\nPhase 3: {len(df_ph3)} total tests")
    print(f"Families: {df_ph3['FDR_Family'].value_counts().to_dict()}")
    sig_raw = (df_ph3["p_value"] < 0.05).sum()
    sig_fdr = (df_ph3["FDR_WithinFamily"] < 0.05).sum()
    print(f"Significant (p<0.05): {sig_raw}")
    print(f"FDR-significant (q<0.05): {sig_fdr}")

    display_cols = ["Phase", "Comparison", "Test", "Group", "N", "N_events",
                    "Statistic", "p_value", "FDR_WithinFamily", "Significant", "Effect_Size"]
    avail = [c for c in display_cols if c in df_ph3.columns]
    display(df_ph3[avail].style.format({
        "p_value": "{:.6f}",
        "FDR_WithinFamily": "{:.4f}",
        "Effect_Size": "{:.4f}"
    }))
else:
    print("No Phase 3 results.")
    df_ph3 = pd.DataFrame()"""))

# ── Cell 19: Combined results + save ──
cells.append(new_code_cell("""# Combine Phase 2 and Phase 3
if len(df_ph2) > 0 and len(df_ph3) > 0:
    combined = pd.concat([df_ph2, df_ph3], ignore_index=True)
elif len(df_ph2) > 0:
    combined = df_ph2.copy()
elif len(df_ph3) > 0:
    combined = df_ph3.copy()
else:
    combined = pd.DataFrame()

# Recompute FDR within each family across merged
if len(combined) > 0 and "FDR_Family" in combined.columns:
    combined["FDR_WithinFamily"] = np.nan
    for fam in combined["FDR_Family"].unique():
        mask = combined["FDR_Family"] == fam
        pvals = combined.loc[mask, "p_value"].values
        if len(pvals) > 1:
            reject, qvals, _, _ = smm.multipletests(pvals, method="fdr_bh")
            combined.loc[mask, "FDR_WithinFamily"] = qvals
        elif len(pvals) == 1:
            combined.loc[mask, "FDR_WithinFamily"] = pvals[0]

# Standardise column order
final_cols = [
    "Phase", "Comparison", "Test", "Group", "N", "N_events",
    "Statistic", "p_value", "FDR_WithinFamily", "Significant", "Effect_Size"
]
combined = combined[final_cols]

out_dir = "/home/alon/menow_home_ass/notebooks/clinical_associations"
os.makedirs(out_dir, exist_ok=True)

# Save full results
combined.to_csv(f"{out_dir}/clinical_associations_results.csv", index=False)

# Save significant summary
sig = combined[combined["p_value"] < 0.05].copy()
sig.to_csv(f"{out_dir}/clinical_associations_summary.csv", index=False)

print(f"\\n{'=' * 60}")
print(f"Combined Results Summary")
print(f"{'=' * 60}")
print(f"Total tests: {len(combined)}")
print(f"  Phase 2: {len(df_ph2)} rows")
print(f"  Phase 3: {len(df_ph3)} rows")
sig_raw = (combined["p_value"] < 0.05).sum()
sig_fdr = (combined["FDR_WithinFamily"] < 0.05).sum()
print(f"Significant (p<0.05): {sig_raw} ({sig_raw/len(combined)*100:.1f}%)")
print(f"FDR-significant (q<0.05): {sig_fdr} ({sig_fdr/len(combined)*100:.1f}%)")

display(combined.style.format({
    "p_value": "{:.6f}",
    "FDR_WithinFamily": "{:.4f}",
    "Effect_Size": "{:.4f}"
}))

print(f"\\nSaved: {out_dir}/clinical_associations_results.csv")
print(f"Saved: {out_dir}/clinical_associations_summary.csv")"""))

# ── Cell 20: Summary markdown ──
cells.append(new_markdown_cell(
    '## Summary\n\n'
    'This analysis characterized relationships between clinical variables '
    'in the PBTA_RNA cohort, independent of outcome.\n\n'
    '**Phase 2 (Cross-Categorical)** identified per-group enrichments of SEX, '
    'cancer predisposition syndromes, molecular subtypes, sex-predisposition interactions, '
    'and race across cancer groups.\n\n'
    '**Phase 3 (Numeric Comparisons)** identified differences in AGE, TUMOR_FRACTION, '
    'and TUMOR_PLOIDY distributions across cancer groups, as well as correlations '
    'between these numeric variables.\n\n'
    '### Next Steps\n'
    '- Review significant findings in context of known biology\n'
    '- Cross-reference with Phase 1 outcome associations\n'
    '- Proceed to multivariate modeling (Phase 4)'
))

nb.cells = cells

with open(NB_PATH, 'w') as f:
    nbformat.write(nb, f)

print(f"✅ Notebook saved to {NB_PATH}")
print(f"   {len(cells)} cells created")
