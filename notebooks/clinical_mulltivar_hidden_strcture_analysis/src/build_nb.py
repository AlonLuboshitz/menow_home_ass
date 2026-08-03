#!/usr/bin/env python3
"""Build the Phase 4+5 multivariate & unsupervised analysis notebook."""

import os, sys
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT_DIR = '/home/alon/menow_home_ass/notebooks/clinical_mulltivar_hidden_strcture_analysis'
NB_PATH = os.path.join(OUT_DIR, 'clinical_mulltivar_hidden_strcture_analysis.ipynb')
os.makedirs(OUT_DIR, exist_ok=True)

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
    '# Phase 4+5: Multivariate Models & Unsupervised Subgroup Discovery'
))

# ── Cell 1: Summary ──
cells.append(new_markdown_cell(
    '## Summary\n\n'
    '**Phase 4: Multivariate Models** — We fit stratified Cox proportional hazards models '
    '(stratified by CANCER_GROUP) to assess whether AGE, SEX, TUMOR_FRACTION, and TUMOR_PLOIDY '
    'independently predict overall survival (OS) and event-free survival (EFS). We also test '
    'whether CANCER_PREDISPOSITIONS adds explanatory power. Forest plots visualise hazard ratios '
    'with 95% confidence intervals, and Schoenfeld residuals test the proportional hazards '
    'assumption.\n\n'
    '**Phase 5: Unsupervised Subgroup Discovery** — We apply PCA, t-SNE, and '
    'FAMD (Factor Analysis of Mixed Data) to clinical features to uncover hidden patient '
    'subgroups. K-means clustering on PCA and FAMD components identifies discrete clusters, '
    'validated via Kaplan-Meier survival curves with log-rank tests. Per-cancer-group analyses '
    'reveal whether distinct subpopulations exist within individual cancer types.'
))

# ── Cell 2: Imports ──
cells.append(new_code_cell("""import sys, os, warnings, math
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly_white"

sys.path.insert(0, "/home/alon/menow_home_ass/notebooks/clinical_analysis/src")
from imports import (
    read_patients, read_samples,
    clean_os, clean_efs, clean_race_eth,
    clean_pred, clean_subtype, clean_tf_tp,
    kaplan_meier, add_km, logrank_multi
)

DATA_DIR = "/home/alon/menow_home_ass/PBTA_RNA"

print("Imports loaded.")"""))

# ── Cell 3: Data loading ──
cells.append(new_code_cell("""pat = read_patients()
smp = read_samples()
df = pat.merge(smp, on="PATIENT_ID", how="inner")
print(f"Merged dataset: {df.shape[0]} rows, {df.shape[1]} columns")

df = clean_os(df)
df = clean_efs(df)
df = clean_race_eth(df)
df = clean_pred(df)
df = clean_subtype(df)
df = clean_tf_tp(df)

for col in ["SEX", "CANCER_GROUP", "CANCER_PREDISPOSITIONS", "RACE", "MOLECULAR_SUBTYPE"]:
    if col in df.columns:
        df[col] = df[col].astype(str)

df = df[df["CANCER_GROUP"] != "nan"].copy()

print("\\nSample sizes per cancer group:")
cg_counts = df["CANCER_GROUP"].value_counts()
for g, n in cg_counts.items():
    print(f"  {g}: {n}")
print(f"\\nTotal unique patients: {df['PATIENT_ID'].nunique()}")
print(f"Total samples: {len(df)}")
print(f"OS events: {df['os_event'].sum():.0f}")
print(f"EFS events: {df['efs_event'].sum():.0f}")"""))

# ── Cell 4: Phase 4 header ──
cells.append(new_markdown_cell(
    '## Phase 4: Multivariate Models'
))

# ── Cell 5: What We're Checking ──
cells.append(new_markdown_cell(
    '## 📌 What We\'re Checking\\n\\n'
    '**Does each variable independently predict OS after adjusting for the others?**'
))

# ── Cell 6: OS Stratified Cox ──
cells.append(new_code_cell("""# ── OS Stratified Cox PH ──
def print_cox_summary(cph, outcome_name):
    hr = np.exp(cph.params_)
    ci = np.exp(cph.confidence_intervals_)
    pvals = cph.summary['p']
    print(f"{'='*70}")
    print(f"  {outcome_name} — Stratified Cox PH")
    print(f"{'='*70}")
    print(f"  {'Variable':<30s} {'HR':>8s}  {'CI.lower':>10s}  {'CI.upper':>10s}  {'p':>8s}")
    print(f"  {'-'*70}")
    for var in hr.index:
        if var.startswith('CANCER_GROUP'):
            continue
        lo = ci.loc[var, ci.columns[0]]
        hi = ci.loc[var, ci.columns[1]]
        pv = pvals[var]
        print(f"  {var:<30s} {hr[var]:>8.4f}  {lo:>10.4f}  {hi:>10.4f}  {pv:>8.4f}")
    print(f"  {'-'*70}")
    print(f"  Concordance: {cph.concordance_index_:.4f}")
    print(f"  AIC: {cph.AIC_partial_:.2f}")
    print(f"  Log-likelihood: {cph.log_likelihood_:.2f}")
    print()

os_df = df.dropna(subset=["OS_MONTHS", "os_event", "AGE", "SEX", "TUMOR_FRACTION", "TUMOR_PLOIDY"]).copy()
print(f"OS model: {len(os_df)} samples, {int(os_df['os_event'].sum())} events")

cph_os = CoxPHFitter(strata=["CANCER_GROUP"])
cph_os.fit(os_df, duration_col="OS_MONTHS", event_col="os_event",
           formula="AGE + C(SEX) + TUMOR_FRACTION + TUMOR_PLOIDY")
print_cox_summary(cph_os, "Overall Survival")

# Test PREDISPOSITIONS
try:
    cph_pred = CoxPHFitter(strata=["CANCER_GROUP"])
    cph_pred.fit(os_df, duration_col="OS_MONTHS", event_col="os_event",
                 formula="AGE + C(SEX) + TUMOR_FRACTION + TUMOR_PLOIDY + C(CANCER_PREDISPOSITIONS)")
    has_any_sig = any(cph_pred.summary['p'] < 0.05)
    print(f"CANCER_PREDISPOSITIONS included. Any category significant? {has_any_sig}")
    if has_any_sig:
        print_cox_summary(cph_pred, "OS with PREDISPOSITIONS")
    else:
        print("No significant PREDISPOSITION categories. Using base model.")
except Exception as e:
    print(f"Could not fit with PREDISPOSITIONS: {e}")
    print("Using model without PREDISPOSITIONS.")"""))

# ── Cell 7: What We're Checking EFS ──
cells.append(new_markdown_cell(
    '## 📌 What We\'re Checking\\n\\n'
    '**Does each variable independently predict EFS after adjusting for the others?**'
))

# ── Cell 8: EFS Stratified Cox ──
cells.append(new_code_cell("""efs_df = df.dropna(subset=["EFS_MONTHS", "efs_event", "AGE", "SEX", "TUMOR_FRACTION", "TUMOR_PLOIDY"]).copy()
print(f"EFS model: {len(efs_df)} samples, {int(efs_df['efs_event'].sum())} events")

cph_efs = CoxPHFitter(strata=["CANCER_GROUP"])
cph_efs.fit(efs_df, duration_col="EFS_MONTHS", event_col="efs_event",
            formula="AGE + C(SEX) + TUMOR_FRACTION + TUMOR_PLOIDY")
print_cox_summary(cph_efs, "Event-Free Survival")"""))

# ── Cell 9: Forest Plot OS ──
cells.append(new_code_cell("""# ── Forest Plot Helper ──
def forest_plot(df_hr, title, subtitle="", x_range=(0.1, 10), ref_line=1.0, ref_label="HR=1"):
    n = len(df_hr)
    fig = go.Figure()
    fig.add_vline(x=ref_line, line=dict(color="gray", width=1, dash="dash"))
    colors = []
    for _, r in df_hr.iterrows():
        if r['p_value'] < 0.05:
            if r['hr'] < 1:
                colors.append('forestgreen')
            else:
                colors.append('crimson')
        else:
            colors.append('lightgray')
    for idx, (_, r) in enumerate(df_hr.iterrows()):
        fig.add_trace(go.Scatter(
            x=[r['hr']], y=[r['label']],
            mode='markers',
            marker=dict(size=10, color=colors[idx], symbol='diamond'),
            error_x=dict(
                type='data', symmetric=False,
                array=[r['ci_upper'] - r['hr']],
                arrayminus=[r['hr'] - r['ci_lower']],
                visible=True, color=colors[idx], thickness=2, width=6,
            ),
            hovertemplate=f'<b>{r["label"]}</b><br>HR: {r["hr"]:.3f}<br>95% CI: ({r["ci_lower"]:.3f}, {r["ci_upper"]:.3f})<extra></extra>',
            showlegend=False,
        ))
    fig.update_xaxes(type='log', range=[np.log10(x_range[0]), np.log10(x_range[1])],
                     title='Hazard Ratio (95% CI)',
                     tickvals=[0.1, 0.2, 0.5, 1, 2, 5, 10],
                     ticktext=['0.1', '0.2', '0.5', '1', '2', '5', '10'])
    fig.update_yaxes(autorange='reversed')
    fig.update_layout(title=title + (f'<br><sup>{subtitle}</sup>' if subtitle else ''),
                      height=max(200, 60 + 40 * n),
                      margin=dict(l=200, r=50, t=80, b=60))
    return fig

# ── OS Forest Plot ──
df_hr_os = []
for var in cph_os.params_.index:
    if var.startswith("CANCER_GROUP"):
        continue
    hr = np.exp(cph_os.params_[var])
    ci = np.exp(cph_os.confidence_intervals_.loc[var])
    pv = cph_os.summary.loc[var, 'p']
    df_hr_os.append({
        'label': var,
        'hr': hr,
        'ci_lower': ci.iloc[0],
        'ci_upper': ci.iloc[1],
        'p_value': pv
    })
df_hr_os = pd.DataFrame(df_hr_os)
fig_os = forest_plot(df_hr_os, "OS: Stratified Cox PH", subtitle="Stratified by CANCER_GROUP")
try:
    fig_os.show()
except:
    pass"""))

# ── Cell 10: Forest Plot EFS ──
cells.append(new_code_cell("""df_hr_efs = []
for var in cph_efs.params_.index:
    if var.startswith("CANCER_GROUP"):
        continue
    hr = np.exp(cph_efs.params_[var])
    ci = np.exp(cph_efs.confidence_intervals_.loc[var])
    pv = cph_efs.summary.loc[var, 'p']
    df_hr_efs.append({
        'label': var,
        'hr': hr,
        'ci_lower': ci.iloc[0],
        'ci_upper': ci.iloc[1],
        'p_value': pv
    })
df_hr_efs = pd.DataFrame(df_hr_efs)
fig_efs = forest_plot(df_hr_efs, "EFS: Stratified Cox PH", subtitle="Stratified by CANCER_GROUP")
try:
    fig_efs.show()
except:
    pass"""))

# ── Cell 11: PH Check ──
cells.append(new_code_cell("""# ── Proportional Hazards Check: Schoenfeld Residuals ──
print("=" * 70)
print("  SCHOENFELD RESIDUAL TEST — OS")
print("=" * 70)
try:
    schoenfeld_os = proportional_hazard_test(cph_os, os_df, time_transform="rank")
    print(schoenfeld_os.summary)
    flagged = [v for v in schoenfeld_os.summary.index if schoenfeld_os.summary.loc[v, 'p'] < 0.05]
    if flagged:
        print(f"\\n⚠️  Variables violating PH assumption (p<0.05): {flagged}")
    else:
        print("\\n✅ No significant PH violation detected.")
except Exception as e:
    print(f"Could not compute Schoenfeld test: {e}")

print("\\n" + "=" * 70)
print("  SCHOENFELD RESIDUAL TEST — EFS")
print("=" * 70)
try:
    schoenfeld_efs = proportional_hazard_test(cph_efs, efs_df, time_transform="rank")
    print(schoenfeld_efs.summary)
    flagged = [v for v in schoenfeld_efs.summary.index if schoenfeld_efs.summary.loc[v, 'p'] < 0.05]
    if flagged:
        print(f"\\n⚠️  Variables violating PH assumption (p<0.05): {flagged}")
    else:
        print("\\n✅ No significant PH violation detected.")
except Exception as e:
    print(f"Could not compute Schoenfeld test: {e}")"""))

# ── Cell 12: What We're Checking subgroup ──
cells.append(new_markdown_cell(
    '## 📌 What We\'re Checking\\n\\n'
    '**Is the effect of each predictor consistent across cancer groups?**'
))

# ── Cell 13: Subgroup forest OS ──
cells.append(new_code_cell("""# ── Subgroup Forest Plot Helper ──
def subgroup_forest(df, outcome_col, event_col, predictors, min_events=5):
    cgs = df.groupby("CANCER_GROUP").filter(lambda g: g[event_col].sum() >= min_events)
    cgs_valid = cgs["CANCER_GROUP"].unique()
    n_predictors = len(predictors)
    fig = make_subplots(
        rows=1, cols=n_predictors,
        subplot_titles=[f"{p} — OS" for p in predictors],
        horizontal_spacing=0.15
    )
    for j, pred in enumerate(predictors):
        results = []
        for cg in sorted(cgs_valid):
            sub = df[df["CANCER_GROUP"] == cg].dropna(subset=[outcome_col, event_col, pred])
            if sub[event_col].sum() < min_events:
                continue
            if pred == "SEX":
                sexes = sub["SEX"].dropna().unique()
                if len(sexes) < 2:
                    continue
                sub = sub[sub["SEX"].isin(sexes)]
                try:
                    cph = CoxPHFitter()
                    cph.fit(sub, duration_col=outcome_col, event_col=event_col, formula=f"C({pred})")
                except Exception:
                    continue
            else:
                try:
                    cph = CoxPHFitter()
                    cph.fit(sub, duration_col=outcome_col, event_col=event_col, formula=pred)
                except Exception:
                    continue
            hr_val = np.exp(cph.params_.iloc[0])
            ci = np.exp(cph.confidence_intervals_.iloc[0])
            pv = cph.summary['p'].iloc[0]
            results.append({
                'group': cg,
                'hr': hr_val,
                'ci_lower': ci.iloc[0],
                'ci_upper': ci.iloc[1],
                'p_value': pv,
                'n': len(sub),
                'events': int(sub[event_col].sum())
            })
        if not results:
            continue
        res_df = pd.DataFrame(results).sort_values('hr')
        colors = ['forestgreen' if r['p_value'] < 0.05 and r['hr'] < 1
                  else 'crimson' if r['p_value'] < 0.05 and r['hr'] > 1
                  else 'lightgray' for _, r in res_df.iterrows()]
        for idx, (_, r) in enumerate(res_df.iterrows()):
            fig.add_trace(go.Scatter(
                x=[r['hr']], y=[f"{r['group'][:20]}"],
                mode='markers',
                marker=dict(size=9, color=colors[idx], symbol='diamond'),
                error_x=dict(
                    type='data', symmetric=False,
                    array=[r['ci_upper'] - r['hr']],
                    arrayminus=[r['hr'] - r['ci_lower']],
                    visible=True, color=colors[idx], thickness=2, width=5,
                ),
                hovertemplate=f'<b>{r["group"][:20]}</b><br>HR: {r["hr"]:.3f}<br>n={r["n"]}, events={r["events"]}<extra></extra>',
                showlegend=False,
            ), row=1, col=j + 1)
        fig.add_vline(x=1.0, line=dict(color="gray", width=1, dash="dash"), row=1, col=j + 1)
        fig.update_xaxes(type='log', range=[np.log10(0.05), np.log10(20)],
                         title='HR (95% CI)', row=1, col=j + 1,
                         tickvals=[0.1, 0.5, 1, 2, 5, 10],
                         ticktext=['0.1', '0.5', '1', '2', '5', '10'])
    fig.update_layout(
        title="Subgroup Forest Plot — OS per Predictor per Cancer Group",
        height=400, width=300 * n_predictors,
        margin=dict(l=50, r=30, t=60, b=40)
    )
    return fig

cox_df_os = df.dropna(subset=["OS_MONTHS", "os_event", "AGE", "SEX", "TUMOR_FRACTION", "TUMOR_PLOIDY"]).copy()
fig_sub_os = subgroup_forest(cox_df_os, "OS_MONTHS", "os_event",
                              predictors=["AGE", "SEX", "TUMOR_FRACTION", "TUMOR_PLOIDY"])
try:
    fig_sub_os.show()
except:
    pass"""))

# ── Cell 14: Subgroup forest EFS ──
cells.append(new_code_cell("""cox_df_efs = df.dropna(subset=["EFS_MONTHS", "efs_event", "AGE", "SEX", "TUMOR_FRACTION", "TUMOR_PLOIDY"]).copy()
fig_sub_efs = subgroup_forest(cox_df_efs, "EFS_MONTHS", "efs_event",
                               predictors=["AGE", "SEX", "TUMOR_FRACTION", "TUMOR_PLOIDY"])
try:
    fig_sub_efs.show()
except:
    pass"""))

# ── Cell 15: Phase 5 header ──
cells.append(new_markdown_cell(
    '## Phase 5: Unsupervised Subgroup Discovery'
))

# ── Cell 16: What We're Checking Phase 5 ──
cells.append(new_markdown_cell(
    '## 📌 What We\'re Checking\\n\\n'
    '**Can we discover hidden patient subgroups using clinical features?**'
))

# ── Cell 17: PCA Numeric ──
cells.append(new_code_cell("""# ── Scatter with Dropdown Helper ──
def scatter_with_dropdown(df, x, y, color_cols, color_labels, title, hover_data=None):
    hover = ['PATIENT_ID'] + (hover_data or [])
    color_maps = {}
    for col in color_cols:
        unique_vals = df[col].dropna().unique()
        palette = px.colors.qualitative.Plotly + px.colors.qualitative.Set2 + px.colors.qualitative.Set3
        color_maps[col] = {v: palette[i % len(palette)] for i, v in enumerate(sorted(unique_vals))}
    default_col = color_cols[0]
    fig = go.Figure()
    hover_text = []
    for _, row in df.iterrows():
        parts = [f"ID: {row.get('PATIENT_ID', '')}"]
        for h in hover:
            if h in row and h != 'PATIENT_ID':
                parts.append(f"{h}: {row[h]}")
        hover_text.append('<br>'.join(parts))
    fig.add_trace(go.Scatter(
        x=df[x], y=df[y],
        mode='markers',
        marker=dict(
            color=[color_maps[default_col].get(v, 'gray') for v in df[default_col]],
            size=6, line=dict(width=0.5, color='gray')
        ),
        text=hover_text,
        hoverinfo='text',
        showlegend=False,
    ))
    buttons = []
    for col, label in zip(color_cols, color_labels):
        colors = [color_maps[col].get(v, 'gray') for v in df[col]]
        buttons.append(dict(
            label=label, method='restyle',
            args=[{'marker.color': [colors]}]
        ))
    fig.update_layout(
        title=title,
        updatemenus=[{
            'buttons': buttons,
            'direction': 'down',
            'showactive': True,
            'x': 1.1, 'y': 1.0,
            'xanchor': 'left', 'yanchor': 'top',
        }],
        width=800, height=600,
    )
    fig.update_xaxes(title=x)
    fig.update_yaxes(title=y)
    return fig

# ── PCA on Numeric Features ──
pca_data = df[["PATIENT_ID", "AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY",
                "CANCER_GROUP", "OS_STATUS", "MOLECULAR_SUBTYPE"]].dropna().copy()
scaler = StandardScaler()
pca_features = scaler.fit_transform(pca_data[["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY"]])
pca = PCA()
pca_components = pca.fit_transform(pca_features)
pca_data["PC1"] = pca_components[:, 0]
pca_data["PC2"] = pca_components[:, 1]
var_exp = pca.explained_variance_ratio_ * 100

# Scree plot
fig_scree = go.Figure()
fig_scree.add_trace(go.Bar(x=[f"PC{i+1}" for i in range(len(var_exp))], y=var_exp,
                            marker_color="steelblue"))
fig_scree.add_trace(go.Scatter(x=[f"PC{i+1}" for i in range(len(var_exp))], y=np.cumsum(var_exp),
                                mode="lines+markers", name="Cumulative",
                                marker=dict(color="crimson"), yaxis="y2"))
fig_scree.update_layout(
    title="PCA Scree Plot — Numeric Features",
    xaxis_title="Principal Component",
    yaxis_title="Variance Explained (%)",
    yaxis2=dict(overlaying="y", side="right", title="Cumulative Variance (%)",
                range=[0, 110]),
    width=600, height=400
)
try:
    fig_scree.show()
except:
    pass

print("Variance explained:")
for i, v in enumerate(var_exp):
    print(f"  PC{i+1}: {v:.2f}% (cumulative: {np.sum(var_exp[:i+1]):.2f}%)")

# 2D scatter
pca_data["OS_EVENT_LABEL"] = pca_data["OS_STATUS"].str.replace(r"^\d+:", "", regex=True)
fig_pca = scatter_with_dropdown(
    pca_data, "PC1", "PC2",
    color_cols=["CANCER_GROUP", "OS_EVENT_LABEL", "MOLECULAR_SUBTYPE"],
    color_labels=["Cancer Group", "OS Status", "Molecular Subtype"],
    title="PCA: 2-Component Projection of Clinical Features",
    hover_data=["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY", "CANCER_GROUP", "OS_EVENT_LABEL"]
)
fig_pca.update_xaxes(title=f"PC1 ({var_exp[0]:.1f}%)")
fig_pca.update_yaxes(title=f"PC2 ({var_exp[1]:.1f}%)")
try:
    fig_pca.show()
except:
    pass"""))

# ── Cell 18: t-SNE Numeric ──
cells.append(new_code_cell("""# ── t-SNE on Numeric Features ──
tsne = TSNE(n_components=2, perplexity=30, random_state=42)
tsne_components = tsne.fit_transform(pca_features)
pca_data["tSNE1"] = tsne_components[:, 0]
pca_data["tSNE2"] = tsne_components[:, 1]

fig_tsne = scatter_with_dropdown(
    pca_data, "tSNE1", "tSNE2",
    color_cols=["CANCER_GROUP", "OS_EVENT_LABEL", "MOLECULAR_SUBTYPE"],
    color_labels=["Cancer Group", "OS Status", "Molecular Subtype"],
    title="t-SNE: 2D Embedding of Clinical Features (perplexity=30)",
    hover_data=["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY", "CANCER_GROUP"]
)
try:
    fig_tsne.show()
except:
    pass"""))

# ── Cell 19: FAMD Mixed ──
cells.append(new_code_cell("""# ── FAMD (Manual): Standardize Numeric + One-Hot Categorical → PCA ──
famd_data = df[["PATIENT_ID", "AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY",
                 "SEX", "CANCER_PREDISPOSITIONS", "RACE", "MOLECULAR_SUBTYPE",
                 "CANCER_GROUP", "OS_STATUS"]].dropna().copy()

# Standardize numeric
scaler_famd = StandardScaler()
numeric_famd = scaler_famd.fit_transform(famd_data[["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY"]])
numeric_cols = ["AGE_scaled", "TF_scaled", "TP_scaled"]
numeric_df = pd.DataFrame(numeric_famd, columns=numeric_cols, index=famd_data.index)

# One-hot encode categorical
cat_cols = ["SEX", "CANCER_PREDISPOSITIONS", "RACE", "MOLECULAR_SUBTYPE"]
cat_df = pd.get_dummies(famd_data[cat_cols], prefix=cat_cols, drop_first=True).astype(float)

# Combine
famd_matrix = pd.concat([numeric_df, cat_df], axis=1)
print(f"FAMD input matrix: {famd_matrix.shape[0]} rows, {famd_matrix.shape[1]} columns")

# PCA on combined matrix
pca_famd = PCA(n_components=2)
famd_pcs = pca_famd.fit_transform(famd_matrix)
famd_data["FAMD1"] = famd_pcs[:, 0]
famd_data["FAMD2"] = famd_pcs[:, 1]
famd_var = pca_famd.explained_variance_ratio_ * 100

famd_data["OS_EVENT_LABEL"] = famd_data["OS_STATUS"].str.replace(r"^\d+:", "", regex=True)

fig_famd = scatter_with_dropdown(
    famd_data, "FAMD1", "FAMD2",
    color_cols=["CANCER_GROUP", "OS_EVENT_LABEL", "MOLECULAR_SUBTYPE"],
    color_labels=["Cancer Group", "OS Status", "Molecular Subtype"],
    title="FAMD: PCA on Standardized Numeric + One-Hot Categorical Features",
    hover_data=["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY", "CANCER_GROUP"]
)
fig_famd.update_xaxes(title=f"FAMD1 ({famd_var[0]:.1f}%)")
fig_famd.update_yaxes(title=f"FAMD2 ({famd_var[1]:.1f}%)")
try:
    fig_famd.show()
except:
    pass"""))

# ── Cell 20: t-SNE on FAMD ──
cells.append(new_code_cell("""# ── t-SNE on FAMD Components ──
# First get enough FAMD components for good t-SNE
pca_famd_full = PCA(n_components=min(20, famd_matrix.shape[1]))
famd_full_pcs = pca_famd_full.fit_transform(famd_matrix)

tsne_famd = TSNE(n_components=2, perplexity=30, random_state=42)
tsne_famd_components = tsne_famd.fit_transform(famd_full_pcs)
famd_data["tSNE_FAMD1"] = tsne_famd_components[:, 0]
famd_data["tSNE_FAMD2"] = tsne_famd_components[:, 1]

fig_tsne_famd = scatter_with_dropdown(
    famd_data, "tSNE_FAMD1", "tSNE_FAMD2",
    color_cols=["CANCER_GROUP", "OS_EVENT_LABEL", "MOLECULAR_SUBTYPE"],
    color_labels=["Cancer Group", "OS Status", "Molecular Subtype"],
    title="t-SNE on First 20 FAMD Components",
    hover_data=["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY", "CANCER_GROUP"]
)
try:
    fig_tsne_famd.show()
except:
    pass"""))

# ── Cell 21: K-means ──
cells.append(new_code_cell("""# ── K-means Clustering: PCA Numeric + FAMD ──
def elbow_silhouette_plot(data, prefix, k_range=range(2, 11)):
    inertias = []
    sil_scores = []
    k_list = []
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = km.fit_predict(data)
        inertias.append(km.inertia_)
        sil = silhouette_score(data, labels)
        sil_scores.append(sil)
        k_list.append(k)
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=(f"Elbow Method — {prefix}", f"Silhouette Score — {prefix}"),
                        horizontal_spacing=0.15)
    fig.add_trace(go.Scatter(x=list(k_range), y=inertias, mode="lines+markers",
                             marker=dict(size=8, color="steelblue")), row=1, col=1)
    fig.add_trace(go.Scatter(x=list(k_range), y=sil_scores, mode="lines+markers",
                             marker=dict(size=8, color="crimson")), row=1, col=2)
    fig.update_xaxes(title="k", row=1, col=1)
    fig.update_yaxes(title="Inertia", row=1, col=1)
    fig.update_xaxes(title="k", row=1, col=2)
    fig.update_yaxes(title="Silhouette Score", row=1, col=2)
    fig.update_layout(title=f"K-means on {prefix} Components", height=400, width=800)
    best_k = k_list[np.argmax(sil_scores)]
    best_sil = max(sil_scores)
    return fig, best_k, best_sil, sil_scores

# PCA numeric
fig_elbow_pca, best_k_pca, best_sil_pca, sil_pca = elbow_silhouette_plot(
    pca_features, "PCA Numeric (3 features)")
try:
    fig_elbow_pca.show()
except:
    pass
print(f"PCA Numeric: Best k={best_k_pca}, Silhouette={best_sil_pca:.4f}")

# FAMD
fig_elbow_famd, best_k_famd, best_sil_famd, sil_famd = elbow_silhouette_plot(
    famd_matrix.values, "FAMD Mixed")
try:
    fig_elbow_famd.show()
except:
    pass
print(f"FAMD Mixed: Best k={best_k_famd}, Silhouette={best_sil_famd:.4f}")

# Assign cluster labels
km_pca = KMeans(n_clusters=best_k_pca, random_state=42, n_init=10)
pca_data["cluster_pca"] = km_pca.fit_predict(pca_features).astype(str)

km_famd = KMeans(n_clusters=best_k_famd, random_state=42, n_init=10)
famd_data["cluster_famd"] = km_famd.fit_predict(famd_matrix.values).astype(str)

print(f"\\nPCA cluster sizes:")
print(pca_data["cluster_pca"].value_counts().to_string())
print(f"\\nFAMD cluster sizes:")
print(famd_data["cluster_famd"].value_counts().to_string())"""))

# ── Cell 22: Survival Validation ──
cells.append(new_code_cell("""# ── Survival Validation: KM per Cluster ──
def km_cluster_validation(data, cluster_col, time_col, event_col, outcome_label, prefix):
    fig = go.Figure()
    colors = px.colors.qualitative.Plotly + px.colors.qualitative.Set2
    clusters = sorted(data[cluster_col].unique())
    groups = []
    for i, cl in enumerate(clusters):
        sub = data[data[cluster_col] == cl]
        km = kaplan_meier(sub[time_col], sub[event_col])
        fig = add_km(fig, km, f"Cluster {cl} (n={len(sub)})", colors[i % len(colors)])
        groups.append((sub[time_col], sub[event_col]))
    p = logrank_multi(groups) if len(groups) >= 2 else 1.0
    fig.update_layout(
        title=f"{outcome_label} by {prefix} Clusters — log-rank p={p:.6f}",
        xaxis_title="Months", yaxis_title="Survival Probability",
        height=500, width=700, template="plotly_white"
    )
    try:
        fig.show()
    except:
        pass
    return p

# Merge cluster labels for validation
pca_merge = pca_data[["PATIENT_ID", "cluster_pca"]].drop_duplicates(subset=["PATIENT_ID"])
famd_merge = famd_data[["PATIENT_ID", "cluster_famd"]].drop_duplicates(subset=["PATIENT_ID"])

val_df = df.merge(pca_merge, on="PATIENT_ID", how="inner").merge(famd_merge, on="PATIENT_ID", how="inner")
print(f"Survival validation dataset: {len(val_df)} rows, {val_df['PATIENT_ID'].nunique()} patients")

print("\\n--- PCA Clusters ---")
p_os_pca = km_cluster_validation(val_df, "cluster_pca", "OS_MONTHS", "os_event", "OS", "PCA")
p_efs_pca = km_cluster_validation(val_df, "cluster_pca", "EFS_MONTHS", "efs_event", "EFS", "PCA")
print(f"  OS log-rank p={p_os_pca:.6f}")
print(f"  EFS log-rank p={p_efs_pca:.6f}")

print("\\n--- FAMD Clusters ---")
p_os_famd = km_cluster_validation(val_df, "cluster_famd", "OS_MONTHS", "os_event", "OS", "FAMD")
p_efs_famd = km_cluster_validation(val_df, "cluster_famd", "EFS_MONTHS", "efs_event", "EFS", "FAMD")
print(f"  OS log-rank p={p_os_famd:.6f}")
print(f"  EFS log-rank p={p_efs_famd:.6f}")"""))

# ── Cell 23: Cluster Profiles ──
cells.append(new_code_cell("""# ── Cluster Profiles ──
def cluster_profiles(data, cluster_col, features_cat, prefix):
    clusters = sorted(data[cluster_col].unique())
    rows = []
    for cl in clusters:
        sub = data[data[cluster_col] == cl]
        n = len(sub)
        age_mean = sub["AGE"].mean()
        age_sd = sub["AGE"].std()
        pct_female = (sub["SEX"] == "Female").mean() * 100 if "SEX" in sub else np.nan
        os_rate = sub["os_event"].mean() * 100 if "os_event" in sub else np.nan
        # Top 3 predispositions
        top_preds = sub["CANCER_PREDISPOSITIONS"].value_counts().head(3)
        top_pred_str = "; ".join([f"{v}({c})" for v, c in top_preds.items()])
        # Race distribution
        race_dist = sub["RACE"].value_counts().to_dict()
        race_str = "; ".join([f"{v}:{c}" for v, c in sorted(race_dist.items(), key=lambda x: -x[1])[:3]])
        # CG composition
        cg_dist = sub["CANCER_GROUP"].value_counts().to_dict()
        cg_str = "; ".join([f"{v}:{c}" for v, c in sorted(cg_dist.items(), key=lambda x: -x[1])[:3]])
        # Top subtypes
        top_subtypes = sub["MOLECULAR_SUBTYPE"].value_counts().head(3)
        subtype_str = "; ".join([f"{v}({c})" for v, c in top_subtypes.items()])
        rows.append({
            "Cluster": cl, "N": n,
            "AGE_mean": f"{age_mean:.1f} ({age_sd:.1f})",
            "%Female": f"{pct_female:.1f}",
            "Top_Predispositions": top_pred_str,
            "Race_Top3": race_str,
            "CG_Top3": cg_str,
            "OS_rate": f"{os_rate:.1f}%",
            "Top_Subtypes": subtype_str
        })
    profile_df = pd.DataFrame(rows)
    print(f"\\n{'='*70}")
    print(f"  CLUSTER PROFILES — {prefix} (k={len(clusters)})")
    print(f"{'='*70}")
    print(profile_df.to_string(index=False))
    return profile_df

# Build merged profile dataset
profile_data = val_df.drop_duplicates(subset=["PATIENT_ID"]).copy()

profile_pca = cluster_profiles(profile_data, "cluster_pca", None, "PCA")
profile_famd = cluster_profiles(profile_data, "cluster_famd", None, "FAMD")"""))

# ── Cell 24: What We're Checking per-group ──
cells.append(new_markdown_cell(
    '## 📌 What We\'re Checking\\n\\n'
    '**Within each cancer type, do distinct subpopulations exist?**'
))

# ── Cell 25: Per-group loop ──
cells.append(new_code_cell("""# ── Per-Group Unsupervised Analysis ──
MIN_CG_SIZE = 50
per_group_results = {}

cgs_large = df["CANCER_GROUP"].value_counts()
cgs_large = cgs_large[cgs_large >= MIN_CG_SIZE].index.tolist()
print(f"Cancer groups with ≥{MIN_CG_SIZE} samples: {len(cgs_large)}")
for cg in cgs_large:
    print(f"  {cg}")

for cg in cgs_large:
    print(f"\\n{'='*60}")
    print(f"  Analyzing: {cg}")
    print(f"{'='*60}")
    
    gdata = df[df["CANCER_GROUP"] == cg].dropna(
        subset=["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY",
                "OS_MONTHS", "os_event", "EFS_MONTHS", "efs_event"]).copy()
    
    if len(gdata) < MIN_CG_SIZE:
        print(f"  Skipping {cg}: only {len(gdata)} samples after dropna")
        continue
    
    # Standardize numeric features
    scaler_g = StandardScaler()
    g_features = scaler_g.fit_transform(gdata[["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY"]])
    
    # PCA
    pca_g = PCA(n_components=2)
    pca_g_comp = pca_g.fit_transform(g_features)
    gdata["PC1"] = pca_g_comp[:, 0]
    gdata["PC2"] = pca_g_comp[:, 1]
    pca_g_var = pca_g.explained_variance_ratio_ * 100
    
    # t-SNE
    tsne_g = TSNE(n_components=2, perplexity=30, random_state=42)
    tsne_g_comp = tsne_g.fit_transform(g_features)
    gdata["tSNE1"] = tsne_g_comp[:, 0]
    gdata["tSNE2"] = tsne_g_comp[:, 1]
    
    # Elbow + Silhouette
    inertias_g = []
    sil_scores_g = []
    k_range = range(2, 11)
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = km.fit_predict(g_features)
        inertias_g.append(km.inertia_)
        sil_scores_g.append(silhouette_score(g_features, labels))
    
    best_k_g = k_range[np.argmax(sil_scores_g)]
    best_sil_g = max(sil_scores_g)
    
    # Plot elbow + silhouette
    fig_eg = make_subplots(rows=1, cols=2,
                           subplot_titles=(f"Elbow — {cg[:25]}", f"Silhouette — {cg[:25]}"),
                           horizontal_spacing=0.2)
    fig_eg.add_trace(go.Scatter(x=list(k_range), y=inertias_g, mode="lines+markers",
                                marker=dict(size=8, color="steelblue")), row=1, col=1)
    fig_eg.add_trace(go.Scatter(x=list(k_range), y=sil_scores_g, mode="lines+markers",
                                marker=dict(size=8, color="crimson")), row=1, col=2)
    fig_eg.update_xaxes(title="k", row=1, col=1)
    fig_eg.update_yaxes(title="Inertia", row=1, col=1)
    fig_eg.update_xaxes(title="k", row=1, col=2)
    fig_eg.update_yaxes(title="Silhouette", row=1, col=2)
    fig_eg.update_layout(title=f"K-means: {cg}", height=400, width=800)
    try:
        fig_eg.show()
    except:
        pass
    
    print(f"  Best k={best_k_g}, Silhouette={best_sil_g:.4f}")
    
    # K-means with best k
    km_g = KMeans(n_clusters=best_k_g, random_state=42, n_init=10)
    gdata["cluster"] = km_g.fit_predict(g_features).astype(str)
    
    # Scatter with dropdown
    gdata["OS_EVENT_LABEL"] = gdata["OS_STATUS"].str.replace(r"^\d+:", "", regex=True)
    
    for proj, xc, yc, ptitle in [("PCA", "PC1", "PC2", f"PCA: {cg}"),
                                   ("t-SNE", "tSNE1", "tSNE2", f"t-SNE: {cg}")]:
        fig_g = scatter_with_dropdown(
            gdata, xc, yc,
            color_cols=["cluster", "OS_EVENT_LABEL", "MOLECULAR_SUBTYPE"],
            color_labels=["Cluster", "OS Status", "Molecular Subtype"],
            title=ptitle,
            hover_data=["AGE", "TUMOR_FRACTION", "TUMOR_PLOIDY"]
        )
        fig_g.update_xaxes(title=f"{xc} ({pca_g_var[0]:.1f}%)" if xc == "PC1" else xc)
        fig_g.update_yaxes(title=f"{yc} ({pca_g_var[1]:.1f}%)" if yc == "PC2" else yc)
        try:
            fig_g.show()
        except:
            pass
    
    per_group_results[cg] = {
        "data": gdata,
        "k": best_k_g,
        "silhouette": best_sil_g,
        "n": len(gdata)
    }

print(f"\\nPer-group analysis complete. {len(per_group_results)} groups analyzed.")"""))

# ── Cell 26: Per-group survival ──
cells.append(new_code_cell("""# ── Per-Group Survival Validation ──
def per_group_survival(gdata, cg, time_col, event_col, outcome_label):
    fig = go.Figure()
    colors = px.colors.qualitative.Plotly + px.colors.qualitative.Set2
    clusters = sorted(gdata["cluster"].unique())
    groups = []
    for i, cl in enumerate(clusters):
        sub = gdata[gdata["cluster"] == cl]
        km = kaplan_meier(sub[time_col], sub[event_col])
        fig = add_km(fig, km, f"Cluster {cl} (n={len(sub)})", colors[i % len(colors)])
        groups.append((sub[time_col], sub[event_col]))
    p = logrank_multi(groups) if len(groups) >= 2 else 1.0
    fig.update_layout(
        title=f"{outcome_label} by Cluster: {cg} — log-rank p={p:.6f}",
        xaxis_title="Months", yaxis_title="Survival Probability",
        height=450, width=650, template="plotly_white"
    )
    try:
        fig.show()
    except:
        pass
    return p

per_group_surv_results = []

for cg, res in per_group_results.items():
    gdata = res["data"]
    print(f"\\n--- {cg} (k={res['k']}, n={res['n']}) ---")
    p_os = per_group_survival(gdata, cg, "OS_MONTHS", "os_event", "OS")
    p_efs = per_group_survival(gdata, cg, "EFS_MONTHS", "efs_event", "EFS")
    print(f"  OS log-rank p={p_os:.6f}")
    print(f"  EFS log-rank p={p_efs:.6f}")
    per_group_surv_results.append({
        "CANCER_GROUP": cg,
        "N": res["n"],
        "k": res["k"],
        "Silhouette": round(res["silhouette"], 4),
        "OS_logrank_p": round(p_os, 6),
        "EFS_logrank_p": round(p_efs, 6)
    })"""))

# ── Cell 27: Per-group summary table ──
cells.append(new_code_cell("""# ── Per-Group Summary Table ──
if per_group_surv_results:
    surv_summary = pd.DataFrame(per_group_surv_results)
    print("\\n" + "=" * 90)
    print("  PER-GROUP CLUSTERING SUMMARY")
    print("=" * 90)
    print(f"  {'Cancer Group':<35s} {'N':>6s} {'k':>4s} {'Silhouette':>11s} {'OS p':>12s} {'EFS p':>12s}")
    print(f"  {'-'*85}")
    for _, r in surv_summary.sort_values("OS_logrank_p").iterrows():
        os_sig = " *" if r["OS_logrank_p"] < 0.05 else ""
        efs_sig = " *" if r["EFS_logrank_p"] < 0.05 else ""
        print(f"  {r['CANCER_GROUP']:<35s} {r['N']:>6d} {r['k']:>4d} {r['Silhouette']:>11.4f} {r['OS_logrank_p']:>10.6f}{os_sig} {r['EFS_logrank_p']:>10.6f}{efs_sig}")
    
    print(f"\\n  Significant (p<0.05): OS={sum(r['OS_logrank_p']<0.05 for r in per_group_surv_results)}, EFS={sum(r['EFS_logrank_p']<0.05 for r in per_group_surv_results)}")
    print()
    
    # Save
    out_dir = '/home/alon/menow_home_ass/notebooks/clinical_mulltivar_hidden_strcture_analysis'
    os.makedirs(out_dir, exist_ok=True)
    surv_summary.to_csv(os.path.join(out_dir, "per_group_clustering_results.csv"), index=False)
    print(f"Saved: {os.path.join(out_dir, 'per_group_clustering_results.csv')}")
else:
    print("No per-group results to summarize.")"""))

# ── Cell 28: Summary markdown ──
cells.append(new_markdown_cell(
    '## Summary\n\n'
    '### Phase 4: Multivariate Models\n'
    '- Stratified Cox PH models were fit for OS and EFS, adjusting for AGE, SEX, TUMOR_FRACTION, and TUMOR_PLOIDY.\n'
    '- Forest plots display hazard ratios with 95% confidence intervals.\n'
    '- Schoenfeld residual tests evaluated the proportional hazards assumption.\n'
    '- Subgroup forest plots per cancer group assessed consistency of predictor effects.\n\n'
    '### Phase 5: Unsupervised Subgroup Discovery\n'
    '- PCA and t-SNE projections of clinical features (AGE, TF, TP) revealed data structure.\n'
    '- FAMD extended the analysis to mixed numeric-categorical data.\n'
    '- K-means clustering identified patient subgroups, validated via survival analysis.\n'
    '- Per-cancer-group analysis highlighted cancer-specific subpopulations.\n\n'
    'Results are saved in the notebook output directory.'
))

# ── Build notebook ──
nb.cells = cells

with open(NB_PATH, "w") as f:
    import json
    json.dump(nb, f, indent=1)

n_md = sum(1 for c in cells if c['cell_type'] == 'markdown')
n_code = sum(1 for c in cells if c['cell_type'] == 'code')
print(f"Written {len(cells)} cells ({n_md} markdown + {n_code} code) to {NB_PATH}")
