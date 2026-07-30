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
    'cancer groups using binomial tests, chi-squared tests, and Fisher exact tests. '
    'Each cancer group gets its own p-value for every comparison. Results are '
    'visualized with stacked bar charts.\n\n'
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
from scipy.stats import chi2_contingency, fisher_exact, binom_test, kruskal, spearmanr, mannwhitneyu, chisquare
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

# ── Cell 3: Data loading + RACE_GROUP ──
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

# Exclude "nan" (missing) cancer group from all analyses
df = df[df["CANCER_GROUP"] != "nan"].copy()

# Create RACE_GROUP broad categories
race_group_map = {
    "White": "White",
    "Black or African American": "Black",
    "Asian": "Asian",
    "More Than One Race": "Other",
    "Other": "Other",
    "American Indian or Alaska Native": "Other",
    "Native Hawaiian or Other Pacific Islander": "Other",
    "Unknown": "Unknown"
}
df["RACE_GROUP"] = df["RACE"].map(race_group_map).fillna("Other")

print("\\nSample sizes per cancer group:")
print(df["CANCER_GROUP"].value_counts().to_string())
print(f"\\nTotal unique patients: {df['PATIENT_ID'].nunique()}")
print(f"Total samples: {len(df)}")"""))

# ── Cell 4: Phase 2 header ──
cells.append(new_markdown_cell(
    '## Phase 2: Cross-Categorical Associations\n\n'
    'Each test produces **per-group p-values** — one row per cancer group per comparison. '
    'Tests include binomial tests (Test 1), chi-squared / Fisher exact (Tests 2, 4), '
    'and descriptive heatmaps (Test 3). FDR is applied within each Comparison family.'
))

# ── Cell 5: Phase 2 helpers ──
cells.append(new_code_cell("""def cramers_v(chi2, n, k, r):
    \"\"\"Cramer's V effect size for chi-squared test. k=rows, r=cols.\"\"\"
    phi2 = chi2 / n
    phi2_corrected = max(0, phi2 - ((k-1)*(r-1))/(n-1))
    k_corrected = k - ((k-1)**2)/(n-1)
    r_corrected = r - ((r-1)**2)/(n-1)
    denom = min(k_corrected - 1, r_corrected - 1)
    if denom <= 0:
        return 0.0
    return math.sqrt(phi2_corrected / denom)

def sig_star(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

results_ph2 = []
print("✅ Phase 2 helpers defined")"""))

# ── Cell 6: Test 1 — SEX enrichment per group vs 50:50 ──
cells.append(new_code_cell("""# ── Test 1: SEX enrichment — Checks: "Is there a sex bias in this cancer group compared to 50:50?"
print("=" * 70)
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

# ── Cell 7: Test 2 — Predisposition profile per group (chi-squared / Fisher) ──
cells.append(new_code_cell("""# ── Test 2: Predisposition profile — Checks: "Does the predisposition makeup of this cancer group differ from all other cancers combined?"
print("=" * 70)
print("Test 2: Predisposition profile per group (chi-squared / Fisher)")
print("=" * 70)

# Filter out unwanted predisposition values
pred_filter = ["No predisposition", "Unknown", "nan", ""]
pred_sub = df[~df["CANCER_PREDISPOSITIONS"].isin(pred_filter)].copy()
print(f"Records after filtering out No predisposition/Unknown: {len(pred_sub)} (from {len(df)})")
print(f"Predisposition categories retained: {pred_sub['CANCER_PREDISPOSITIONS'].nunique()}")

pred_counts = pred_sub["CANCER_PREDISPOSITIONS"].value_counts()
print("\\nPredisposition distribution:")
for v, c in pred_counts.items():
    print(f"  {str(v)[:60]:60s} N={c}")

min_n = 20
pred_results = []
for group in sorted(pred_sub["CANCER_GROUP"].unique()):
    g = pred_sub[pred_sub["CANCER_GROUP"] == group]
    n_group = len(g)
    if n_group < min_n:
        continue
    other = pred_sub[pred_sub["CANCER_GROUP"] != group]
    # Build distributions
    g_dist = g["CANCER_PREDISPOSITIONS"].value_counts()
    other_dist = other["CANCER_PREDISPOSITIONS"].value_counts()
    # Union of categories
    all_cats = sorted(set(g_dist.index) | set(other_dist.index))
    g_counts = np.array([g_dist.get(c, 0) for c in all_cats])
    other_counts = np.array([other_dist.get(c, 0) for c in all_cats])
    # Contingency table: 2 x k
    table = np.array([g_counts, other_counts])
    # Chi-squared test
    try:
        expected_v = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
        has_small_expected = (expected_v < 5).any()
        if has_small_expected and table.shape[1] == 2:
            # 2x2 table with small expected → Fisher exact
            odds_ratio, p = fisher_exact(table)
            test_name = "Fisher exact"
        else:
            chi2, p, dof, exp = chi2_contingency(table)
            test_name = "Chi-squared"
            if has_small_expected:
                test_name += " *"
        # Compute effect size
        if test_name.startswith("Chi-squared"):
            cv = cramers_v(chi2, table.sum(), table.shape[0], table.shape[1])
            effect_val = round(cv, 4)
        else:
            effect_val = round(math.log(odds_ratio), 4) if odds_ratio > 0 else 0
    except Exception as e:
        p = 1.0
        test_name = f"Error: {e}"
        effect_val = 0

    top3 = sorted(g_dist.items(), key=lambda x: x[1], reverse=True)[:3]
    top3_str = "; ".join(f"{k} ({v})" for k, v in top3)
    pred_results.append({
        "group": group, "n": n_group,
        "n_pred_cats": len(g_dist), "p": p,
        "test": test_name, "top3": top3_str
    })
    results_ph2.append({
        "FDR_Family": "Predisposition profile",
        "Phase": "Phase 2",
        "Comparison": "Predisposition profile",
        "Test": test_name,
        "Group": group,
        "N": n_group,
        "N_events": n_group,
        "Statistic": f"{len(g_dist)} cats",
        "p_value": round(p, 6),
        "Effect_Size": effect_val,
        "Significant": sig_star(p)
    })
    print(f"  {group:40s} N={n_group:4d}  cats={len(g_dist):2d}  p={p:.6f}  {sig_star(p)}  Top3: {top3_str}")

# Stacked bar chart: predisposition composition per group
if pred_results:
    # Build a dataframe suitable for stacked bar
    groups_ok = [r["group"] for r in pred_results]
    pred_stack = pred_sub[pred_sub["CANCER_GROUP"].isin(groups_ok)]
    ct = pd.crosstab(pred_stack["CANCER_GROUP"], pred_stack["CANCER_PREDISPOSITIONS"])
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
    # Keep top 15 predispositions for readability
    top_cats = pred_counts.head(15).index.tolist()
    ct_plot = ct_pct[[c for c in top_cats if c in ct_pct.columns]]
    other_cols = [c for c in ct_pct.columns if c not in ct_plot.columns]
    if other_cols:
        ct_plot["Other"] = ct_pct[other_cols].sum(axis=1)
    ct_plot = ct_plot.fillna(0)
    fig = go.Figure()
    for cat in ct_plot.columns:
        fig.add_trace(go.Bar(
            name=str(cat)[:50], x=ct_plot.index, y=ct_plot[cat],
            hovertemplate=f"%{{y:.1f}}%<extra>{str(cat)[:50]}</extra>"
        ))
    fig.update_layout(
        barmode="stack", title="Predisposition composition per cancer group (%)",
        xaxis_title="Cancer Group", yaxis_title="Percentage",
        height=max(400, 20 * len(ct_plot)), width=1000,
        legend=dict(font=dict(size=8))
    )
    try:
        fig.show()
    except Exception:
        pass
    print(f"\\nTest 2 complete: {len(pred_results)} groups tested")"""))

# ── Cell 8: Test 3 — SUBTYPE descriptive heatmap ──
cells.append(new_code_cell("""# ── Test 3: SUBTYPE descriptive — Purpose: Show subtype distribution across cancer groups (descriptive only)
print("=" * 70)
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

# ── Cell 9: Test 4 — Race distribution per group (chi-squared GoF) ──
cells.append(new_code_cell("""# ── Test 4: Race distribution — Checks: "Does the racial composition of this cancer group differ from the overall cohort's racial composition?"
print("=" * 70)
print("Test 4: Race distribution per group (chi-squared goodness-of-fit)")
print("=" * 70)

race_order = ["White", "Black", "Asian", "Other", "Unknown"]
overall_race = df["RACE_GROUP"].value_counts()
# Ensure all categories present
for cat in race_order:
    if cat not in overall_race.index:
        overall_race[cat] = 0
overall_race = overall_race[race_order]
overall_props = overall_race / overall_race.sum()
print(f"Overall race distribution (N={overall_race.sum():.0f}):")
for cat in race_order:
    print(f"  {cat:10s} {overall_race[cat]:5d} ({overall_props[cat]*100:5.1f}%)")

race_results = []
race_plot_data = []
for group in sorted(df["CANCER_GROUP"].unique()):
    g = df[df["CANCER_GROUP"] == group]
    n = len(g)
    if n < 20:
        continue
    within = g["RACE_GROUP"].value_counts()
    for cat in race_order:
        if cat not in within.index:
            within[cat] = 0
    within = within[race_order]
    obs = within.values.astype(float)
    exp = (overall_props.values * n).astype(float)
    # Only include categories with expected > 0
    mask = exp > 0
    if mask.sum() < 2:
        continue
    try:
        chi2_stat, p = chisquare(f_obs=obs[mask], f_exp=exp[mask])
        test_name = "Chi2 GoF"
    except Exception as e:
        p = 1.0
        chi2_stat = 0
        test_name = f"Error: {e}"

    cohens_w = math.sqrt(chi2_stat / n) if n > 0 else 0
    race_results.append({
        "group": group, "n": n, "p": p, "test": test_name
    })
    results_ph2.append({
        "FDR_Family": "Race distribution",
        "Phase": "Phase 2",
        "Comparison": "Race distribution",
        "Test": test_name,
        "Group": group,
        "N": n,
        "N_events": n,
        "Statistic": f"chi2={chi2_stat:.2f}",
        "p_value": round(p, 6),
        "Effect_Size": round(cohens_w, 4),
        "Significant": sig_star(p)
    })
    print(f"  {group:40s} N={n:4d}  p={p:.6f}  {sig_star(p)}")
    race_plot_data.append({
        "Group": group, "Type": "Within-group",
        **{cat: within[cat] / n * 100 for cat in race_order}
    })
    race_plot_data.append({
        "Group": group, "Type": "Overall cohort",
        **{cat: overall_props[cat] * 100 for cat in race_order}
    })

# Stacked bar: per group vs overall
if race_plot_data:
    rdf = pd.DataFrame(race_plot_data)
    fig = go.Figure()
    for cat in race_order:
        for tp in ["Within-group", "Overall cohort"]:
            sub = rdf[(rdf["Type"] == tp)]
            fig.add_trace(go.Bar(
                name=f"{cat} ({tp})",
                x=sub["Group"] + " " + tp,
                y=sub[cat],
                legendgroup=cat,
                showlegend=(sub["Type"].iloc[0] == tp),
                marker=dict(
                    line=dict(width=1, color="black")
                ) if tp == "Overall cohort" else {}
            ))
    fig.update_layout(
        barmode="stack",
        title="Race composition: within-group vs overall cohort",
        xaxis_title="Cancer Group", yaxis_title="Percentage",
        height=max(400, 25 * len(race_results)), width=1000,
        legend=dict(font=dict(size=9))
    )
    try:
        fig.show()
    except Exception:
        pass
    print(f"\\nTest 4 complete: {len(race_results)} groups tested")"""))

# ── Cell 10: Phase 2 FDR + display ──
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

# ── Cell 11: Phase 3 header ──
cells.append(new_markdown_cell(
    '## Phase 3: Numeric Comparisons\n\n'
    'Tests: Kruskal-Wallis (global across all cancer groups), per-group '
    'Mann-Whitney (each group vs all others), and Spearman correlations '
    'between numeric variables. FDR is applied within each Comparison family.'
))

# ── Cell 12: Phase 3 helpers ──
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

# ── Cell 13: Test 5 — AGE × CANCER_GROUP ──
cells.append(new_code_cell("""# ── Test 5: AGE × CG — Checks: "Does AGE differ across cancer groups? Which groups are outliers?"
print("=" * 70)
print("Test 5: AGE x CANCER_GROUP")
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

# ── Cell 14: Test 6 — TF × CANCER_GROUP ──
cells.append(new_code_cell("""# ── Test 6: TF × CG — Checks: "Does TUMOR_FRACTION differ across cancer groups? Which groups are outliers?"
print("=" * 70)
print("Test 6: TUMOR_FRACTION x CANCER_GROUP")
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

# ── Cell 15: Test 7 — TP × CANCER_GROUP ──
cells.append(new_code_cell("""# ── Test 7: TP × CG — Checks: "Does TUMOR_PLOIDY differ across cancer groups? Which groups are outliers?"
print("=" * 70)
print("Test 7: TUMOR_PLOIDY x CANCER_GROUP")
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

# ── Cell 16: Tests 8-10 — Correlations ──
cells.append(new_code_cell("""# ── Test 8: AGE × TF — Checks: "Is AGE associated with TUMOR_FRACTION?"
# ── Test 9: AGE × TP — Checks: "Is AGE associated with TUMOR_PLOIDY?"
# ── Test 10: TF × TP — Checks: "Is TUMOR_FRACTION associated with TUMOR_PLOIDY?"
print("=" * 70)
print("Tests 8-10: Numeric correlations (Spearman + LOESS scatter)")
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

# ── Cell 17: Phase 3 FDR + display ──
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

# ── Cell 18: Combined results + save ──
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

# ── Cell 19: Summary markdown ──
cells.append(new_markdown_cell(
    '## Summary\n\n'
    'This analysis characterized relationships between clinical variables '
    'in the PBTA_RNA cohort, independent of outcome.\n\n'
    '**Phase 2 (Cross-Categorical)** identified per-group enrichments of SEX, '
    'cancer predisposition syndromes, molecular subtype profiles, '
    'and race distributions across cancer groups.\n\n'
    '**Phase 3 (Numeric Comparisons)** identified differences in AGE, TUMOR_FRACTION, '
    'and TUMOR_PLOIDY distributions across cancer groups, as well as correlations '
    'between these numeric variables.\n\n'
    '### Next Steps\n'
    '- Review significant findings in context of known biology\n'
    '- Cross-reference with Phase 1 outcome associations\n'
    '- Proceed to multivariate modeling (Phase 4)'
))

# ── Insert markdown cells before each test code cell ──

cells.insert(16, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Are AGE, TUMOR_FRACTION, and TUMOR_PLOIDY correlated with each other?**\n\n'
    'Spearman rank correlations between all three numeric variable pairs (AGE×TF, AGE×TP, TF×TP). Scatter plots with LOESS trendline. Report rho and p-value.'
))

cells.insert(15, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Does TUMOR_PLOIDY differ across cancer groups? Which groups are outliers?**\n\n'
    'Global Kruskal-Wallis test across all cancer groups with n ≥ 20, followed by per-group Mann-Whitney (each group vs all others combined). Boxplot with jittered points.'
))

cells.insert(14, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Does TUMOR_FRACTION differ across cancer groups? Which groups are outliers?**\n\n'
    'Global Kruskal-Wallis test across all cancer groups with n ≥ 20, followed by per-group Mann-Whitney (each group vs all others combined). Boxplot with jittered points.'
))

cells.insert(13, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Does AGE differ across cancer groups? Which groups are outliers?**\n\n'
    'Global Kruskal-Wallis test across all cancer groups with n ≥ 20, followed by per-group Mann-Whitney (each group vs all others combined). Report epsilon-squared for KW, rank-biserial r for MW. Boxplot with jittered points.'
))

cells.insert(9, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Does the racial composition of this cancer group differ from the overall cohort\'s racial composition?**\n\n'
    'For each cancer group with n ≥ 20, perform a chi-squared goodness-of-fit test comparing the within-group race distribution (White/Black/Asian/Other/Unknown) against the overall cohort proportions. Report Cohen\'s w as effect size.'
))

cells.insert(8, new_markdown_cell(
    '## 📌 Purpose\n'
    '**Show the molecular subtype distribution across cancer groups (descriptive only, no tests performed).**\n\n'
    'Display a heatmap of MOLECULAR_SUBTYPE × CANCER_GROUP counts. This is purely descriptive to visualize which subtypes appear in which cancer groups.'
))

cells.insert(7, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Does the predisposition makeup of this cancer group differ from all other cancers combined?**\n\n'
    'Exclude "No predisposition" and "Unknown" patients. For each cancer group with n ≥ 20 (after filtering), compare the distribution of predisposition categories within that group vs all other groups combined using a chi-squared test (Fisher exact if expected < 5). Report Cramer\'s V as effect size.'
))

cells.insert(6, new_markdown_cell(
    '## 📌 What We\'re Checking\n'
    '**Is there a sex bias in this cancer group compared to 50:50?**\n\n'
    'For each cancer group with n ≥ 20, perform a two-sided binomial test comparing the observed female proportion to the expected 50%. Report direction of enrichment (Male/Female) and effect size (difference from 50%).'
))

print(f"  Inserted 8 markdown cells → total {len(cells)} cells")

nb.cells = cells

with open(NB_PATH, 'w') as f:
    nbformat.write(nb, f)

print(f"✅ Notebook saved to {NB_PATH}")
print(f"   {len(cells)} cells created")
