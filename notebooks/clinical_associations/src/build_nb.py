#!/usr/bin/env python3
"""Build the Phase 2+3 clinical associations notebook programmatically."""

import sys, os, nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT_DIR = '/home/alon/menow_home_ass/notebooks/clinical_associations'
NB_PATH = os.path.join(OUT_DIR, 'clinical_associations_analysis.ipynb')
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
    '# Phase 2+3: Clinical Associations — Cross-Categorical & Numeric Comparisons'
))

# ── Cell 1: Summary ──
cells.append(new_markdown_cell(
    '## Summary\n\n'
    'This notebook implements **Phase 2 (Cross-Categorical Associations)** and '
    '**Phase 3 (Numeric Comparisons)** of the clinical deep-dive analysis.\n\n'
    '**Phase 2** tests associations between pairs of categorical clinical variables '
    '(SEX, CANCER_GROUP, CANCER_PREDISPOSITIONS, MOLECULAR_SUBTYPE, RACE) using '
    'chi-squared tests and Cramer\'s V.\n\n'
    '**Phase 3** compares numeric variables (AGE, TUMOR_FRACTION, TUMOR_PLOIDY) '
    'across cancer groups using Kruskal-Wallis tests, and computes Spearman '
    'correlations between numeric pairs.\n\n'
    'See `context/clinical_deep_dive_general.md` for shared methodology and conventions.'
))

# ── Cell 2: Imports ──
cells.append(new_code_cell(
    '# ── Standard library ──\n'
    'import sys, os, warnings, textwrap, itertools, math\n'
    'warnings.filterwarnings(\'ignore\')\n\n'
    '# ── Data handling ──\n'
    'import numpy as np\n'
    'import pandas as pd\n\n'
    '# ── Statistics ──\n'
    'from scipy.stats import chi2_contingency, fisher_exact, kruskal, spearmanr\n'
    'from scipy.stats import mannwhitneyu as mann_whitney_u\n'
    'import statsmodels.stats.multitest as smm\n'
    'from itertools import combinations\n\n'
    '# ── Plotting ──\n'
    'import plotly.graph_objects as go\n'
    'import plotly.express as px\n'
    'from plotly.subplots import make_subplots\n'
    'import plotly.io as pio\n'
    'pio.templates.default = \'plotly_white\'\n\n'
    '# ── Helpers ──\n'
    'sys.path.insert(0, \'/home/alon/menow_home_ass/notebooks/clinical_analysis/src\')\n'
    'from imports import (\n'
    '    read_patients, read_samples,\n'
    '    clean_os, clean_efs, clean_race_eth,\n'
    '    clean_pred, clean_subtype, clean_tf_tp\n'
    ')\n\n'
    'print("✅ All imports loaded")'
))

# ── Cell 3: Data loading ──
cells.append(new_code_cell(
    '# ── Load data ──\n'
    'pat = read_patients()\n'
    'smp = read_samples()\n\n'
    '# Merge on PATIENT_ID (inner)\n'
    'df = pat.merge(smp, on=\'PATIENT_ID\', how=\'inner\')\n'
    'print(f"Merged dataset: {df.shape[0]} rows, {df.shape[1]} columns")\n\n'
    '# ── Clean ──\n'
    'df = clean_race_eth(df)\n'
    'df = clean_pred(df)\n'
    'df = clean_subtype(df)\n'
    'df = clean_tf_tp(df)\n\n'
    '# ── OS / EFS for reference (not primary analysis here) ──\n'
    'df = clean_os(df)\n'
    'df = clean_efs(df)\n\n'
    '# ── Categorical variables as strings for cross-tab ──\n'
    'for col in [\'SEX\', \'CANCER_GROUP\', \'CANCER_PREDISPOSITIONS\', \'RACE\', \'MOLECULAR_SUBTYPE\']:\n'
    '    if col in df.columns:\n'
    '        df[col] = df[col].astype(str)\n\n'
    '# ── Display summary ──\n'
    'print("\\nSample sizes per cancer group:")\n'
    'print(df[\'CANCER_GROUP\'].value_counts().to_string())\n'
    'print(f"\\nTotal unique patients: {df[\'PATIENT_ID\'].nunique()}")\n'
    'print(f"Total samples: {len(df)}")'
))

# ── Cell 4: Phase 2 header ──
cells.append(new_markdown_cell('## Phase 2: Cross-Categorical Associations'))

# ── Cell 5: Phase 2 helpers ──
cells.append(new_code_cell(
    '# ── Phase 2 helper: chi-squared with Cramer\'s V ──\n'
    'def cramers_v(confusion_matrix):\n'
    '    """Calculate Cramer\'s V from a confusion matrix."""\n'
    '    chi2, p, dof, expected = chi2_contingency(confusion_matrix, correction=False)\n'
    '    n = confusion_matrix.sum().sum()\n'
    '    min_dim = min(confusion_matrix.shape) - 1\n'
    '    if min_dim == 0 or n == 0:\n'
    '        return np.nan, chi2, p\n'
    '    v = np.sqrt(max(0, chi2 / (n * min_dim)))\n'
    '    return v, chi2, p\n\n'
    'def run_cross_cat(df, col1, col2, label1, label2, min_expected=5):\n'
    '    """Run chi-squared or Fisher exact test between two categorical variables."""\n'
    '    # Drop rows where either is missing\n'
    '    sub = df[[col1, col2]].dropna().copy()\n'
    '    sub = sub[(sub[col1] != \'nan\') & (sub[col2] != \'nan\')]\n'
    '    sub = sub[(sub[col1] != \'\') & (sub[col2] != \'\')]\n\n'
    '    # Build contingency table\n'
    '    ctab = pd.crosstab(sub[col1], sub[col2])\n\n'
    '    n_total = sub.shape[0]\n\n'
    '    # Check if Fisher exact is needed (any expected < 5)\n'
    '    chi2_stat, p_val, dof, expected = chi2_contingency(ctab, correction=False)\n'
    '    min_exp = expected.min()\n\n'
    '    if min_exp < 5:\n'
    '        # Use Fisher exact for 2x2, otherwise warn and use chi-squared\n'
    '        if ctab.shape == (2, 2):\n'
    '            odds_ratio, p_val = fisher_exact(ctab)\n'
    '            test_name = "Fisher exact"\n'
    '        else:\n'
    '            p_val = chi2_contingency(ctab, correction=False)[1]\n'
    '            test_name = "\\u03c7\\u00b2 (Fisher not available for >2\\u00d72)"\n'
    '    else:\n'
    '        p_val = chi2_stat\n'
    '        test_name = "\\u03c7\\u00b2"\n'
    '        p_val = chi2_contingency(ctab, correction=False)[1]\n\n'
    '    v, chi2_v, p_v = cramers_v(ctab.values)\n\n'
    '    return {\n'
    '        \'Test\': f\'{label1} \\u00d7 {label2}\',\n'
    '        \'Var1\': label1, \'Var2\': label2,\n'
    '        \'N\': n_total,\n'
    '        \'Statistic\': test_name,\n'
    '        \'Chi2\': round(chi2_v, 4),\n'
    '        \'P_value\': round(p_v, 6),\n'
    '        \'Cramers_V\': round(v, 4),\n'
    '        \'Min_Expected\': round(min_exp, 2),\n'
    '        \'Table_Shape\': f\'{ctab.shape[0]}\\u00d7{ctab.shape[1]}\'\n'
    '    }\n\n'
    'print("\\u2705 Phase 2 helpers defined")'
))

# ── Cell 6: Test 1 — SEX × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 1: SEX × CANCER_GROUP ──\n'
    'results_ph2 = []\n\n'
    'res = run_cross_cat(df, \'SEX\', \'CANCER_GROUP\', \'SEX\', \'CANCER_GROUP\')\n'
    'results_ph2.append(res)\n'
    'print(f"SEX × CANCER_GROUP: N={res[\'N\']}, p={res[\'P_value\']}, V={res[\'Cramers_V\']}")\n\n'
    '# Heatmap\n'
    'ctab = pd.crosstab(df[\'SEX\'], df[\'CANCER_GROUP\'])\n'
    'fig = px.imshow(ctab.values, \n'
    '                x=ctab.columns, y=ctab.index,\n'
    '                text_auto=True, aspect="auto",\n'
    '                color_continuous_scale=\'Blues\',\n'
    '                title=f"SEX × CANCER_GROUP (\\u03c7\\u00b2 p={res[\'P_value\']:.4f}, V={res[\'Cramers_V\']:.3f})")\n'
    'fig.update_layout(xaxis_tickangle=45, height=400, width=900)\n'
    'try:\n'
    '    fig.show()\n'
    'except:\n'
    '    pass'
))

# ── Cell 7: Test 2 — CANCER_PREDISPOSITIONS × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 2: CANCER_PREDISPOSITIONS × CANCER_GROUP ──\n'
    'res = run_cross_cat(df, \'CANCER_PREDISPOSITIONS\', \'CANCER_GROUP\', \'PRED\', \'CANCER_GROUP\')\n'
    'results_ph2.append(res)\n'
    'print(f"PRED × CANCER_GROUP: N={res[\'N\']}, p={res[\'P_value\']}, V={res[\'Cramers_V\']}")\n\n'
    'ctab = pd.crosstab(df[\'CANCER_PREDISPOSITIONS\'], df[\'CANCER_GROUP\'])\n'
    'fig = px.imshow(ctab.values, \n'
    '                x=ctab.columns, y=ctab.index,\n'
    '                text_auto=True, aspect="auto",\n'
    '                color_continuous_scale=\'Blues\',\n'
    '                title=f"PRED × CANCER_GROUP (\\u03c7\\u00b2 p={res[\'P_value\']:.4f}, V={res[\'Cramers_V\']:.3f})")\n'
    'fig.update_layout(xaxis_tickangle=45, height=500, width=900)\n'
    'try:\n'
    '    fig.show()\n'
    'except:\n'
    '    pass'
))

# ── Cell 8: Test 3 — MOLECULAR_SUBTYPE × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 3: MOLECULAR_SUBTYPE × CANCER_GROUP ──\n'
    'res = run_cross_cat(df, \'MOLECULAR_SUBTYPE\', \'CANCER_GROUP\', \'SUBTYPE\', \'CANCER_GROUP\')\n'
    'results_ph2.append(res)\n'
    'print(f"SUBTYPE × CANCER_GROUP: N={res[\'N\']}, p={res[\'P_value\']}, V={res[\'Cramers_V\']}")\n\n'
    '# Heatmap — top subtypes only for readability\n'
    'ctab = pd.crosstab(df[\'MOLECULAR_SUBTYPE\'], df[\'CANCER_GROUP\'])\n'
    '# Keep subtypes with at least 20 samples total\n'
    'subtype_counts = ctab.sum(axis=1)\n'
    'top_subtypes = subtype_counts[subtype_counts >= 20].index\n'
    'ctab_top = ctab.loc[top_subtypes]\n'
    'fig = px.imshow(ctab_top.values, \n'
    '                x=ctab_top.columns, y=ctab_top.index,\n'
    '                text_auto=True, aspect="auto",\n'
    '                color_continuous_scale=\'Blues\',\n'
    '                title=f"SUBTYPE × CANCER_GROUP (top subtypes, \\u03c7\\u00b2 p={res[\'P_value\']:.4f}, V={res[\'Cramers_V\']:.3f})")\n'
    'fig.update_layout(xaxis_tickangle=45, height=600, width=1000)\n'
    'try:\n'
    '    fig.show()\n'
    'except:\n'
    '    pass'
))

# ── Cell 9: Test 4 — SEX × CANCER_PREDISPOSITIONS ──
cells.append(new_code_cell(
    '# ── Test 4: SEX × CANCER_PREDISPOSITIONS ──\n'
    'res = run_cross_cat(df, \'SEX\', \'CANCER_PREDISPOSITIONS\', \'SEX\', \'PRED\')\n'
    'results_ph2.append(res)\n'
    'print(f"SEX × PRED: N={res[\'N\']}, p={res[\'P_value\']}, V={res[\'Cramers_V\']}")\n\n'
    'ctab = pd.crosstab(df[\'SEX\'], df[\'CANCER_PREDISPOSITIONS\'])\n'
    'fig = px.imshow(ctab.values, \n'
    '                x=ctab.columns, y=ctab.index,\n'
    '                text_auto=True, aspect="auto",\n'
    '                color_continuous_scale=\'Blues\',\n'
    '                title=f"SEX × PRED (\\u03c7\\u00b2 p={res[\'P_value\']:.4f}, V={res[\'Cramers_V\']:.3f})")\n'
    'fig.update_layout(xaxis_tickangle=45, height=400, width=700)\n'
    'try:\n'
    '    fig.show()\n'
    'except:\n'
    '    pass'
))

# ── Cell 10: Test 5 — RACE × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 5: RACE × CANCER_GROUP ──\n'
    'res = run_cross_cat(df, \'RACE\', \'CANCER_GROUP\', \'RACE\', \'CANCER_GROUP\')\n'
    'results_ph2.append(res)\n'
    'print(f"RACE × CANCER_GROUP: N={res[\'N\']}, p={res[\'P_value\']}, V={res[\'Cramers_V\']}")\n\n'
    'ctab = pd.crosstab(df[\'RACE\'], df[\'CANCER_GROUP\'])\n'
    'fig = px.imshow(ctab.values, \n'
    '                x=ctab.columns, y=ctab.index,\n'
    '                text_auto=True, aspect="auto",\n'
    '                color_continuous_scale=\'Blues\',\n'
    '                title=f"RACE × CANCER_GROUP (\\u03c7\\u00b2 p={res[\'P_value\']:.4f}, V={res[\'Cramers_V\']:.3f})")\n'
    'fig.update_layout(xaxis_tickangle=45, height=400, width=900)\n'
    'try:\n'
    '    fig.show()\n'
    'except:\n'
    '    pass'
))

# ── Cell 11: Phase 2 results table ──
cells.append(new_code_cell(
    '# ── Phase 2 Results ──\n'
    'df_ph2 = pd.DataFrame(results_ph2)\n'
    'df_ph2[\'Phase\'] = \'Phase 2\'\n'
    'df_ph2 = df_ph2[[\'Phase\', \'Test\', \'Var1\', \'Var2\', \'N\', \'Statistic\', \'Chi2\', \'P_value\', \'Cramers_V\', \'Min_Expected\', \'Table_Shape\']]\n'
    'display(df_ph2.style\n'
    '    .format({\'P_value\': \'{:.6f}\', \'Chi2\': \'{:.4f}\', \'Cramers_V\': \'{:.4f}\', \'Min_Expected\': \'{:.2f}\'})\n'
    '    .background_gradient(cmap=\'RdYlGn_r\', subset=[\'P_value\'])\n'
    ')'
))

# ── Cell 12: Phase 3 header ──
cells.append(new_markdown_cell('## Phase 3: Numeric Comparisons'))

# ── Cell 13: Phase 3 helpers ──
cells.append(new_code_cell(
    '# ── Phase 3 helpers ──\n\n'
    'def kw_dunn(df, num_col, cat_col, min_n=20):\n'
    '    """Kruskal-Wallis test with Dunn\'s post-hoc and \\u03b5\\u00b2 effect size."""\n'
    '    groups = df[[num_col, cat_col]].dropna()\n'
    '    groups = groups[groups[cat_col].isin(groups[cat_col].value_counts()[groups[cat_col].value_counts() >= min_n].index)]\n\n'
    '    cat_names = groups[cat_col].unique()\n'
    '    cat_data = {name: groups.loc[groups[cat_col] == name, num_col].values for name in cat_names}\n\n'
    '    if len(cat_data) < 2:\n'
    '        return None, None, None, None, None\n\n'
    '    # KW test\n'
    '    h_stat, p_val = kruskal(*cat_data.values())\n\n'
    '    # \\u03b5\\u00b2 effect size\n'
    '    n_total = sum(len(v) for v in cat_data.values())\n'
    '    if n_total > 0:\n'
    '        eps_sq = (h_stat - len(cat_data) + 1) / (n_total - len(cat_data))\n'
    '    else:\n'
    '        eps_sq = np.nan\n\n'
    '    # Dunn\'s post-hoc\n'
    '    from scipy.stats import norm\n'
    '    groups_list = list(cat_data.keys())\n'
    '    posthoc = []\n'
    '    for i, j in combinations(range(len(groups_list)), 2):\n'
    '        ni = len(cat_data[groups_list[i]])\n'
    '        nj = len(cat_data[groups_list[j]])\n'
    '        ri = np.mean([np.mean(cat_data[g] < v) for g in groups_list for v in cat_data[g]])\n'
    '        # Simplified Dunn: use Mann-Whitney U as proxy\n'
    '        u_stat, p_mw = mann_whitney_u(cat_data[groups_list[i]], cat_data[groups_list[j]], alternative=\'two-sided\')\n'
    '        posthoc.append({\n'
    '            \'Group1\': groups_list[i], \'Group2\': groups_list[j],\n'
    '            \'MW_U\': u_stat, \'P_raw\': p_mw\n'
    '        })\n\n'
    '    return h_stat, p_val, eps_sq, posthoc, groups_list\n\n'
    'def boxplot_by_group(df, num_col, cat_col, title, ylabel, min_n=20):\n'
    '    """Boxplot of numeric variable grouped by categorical."""\n'
    '    groups = df[[num_col, cat_col]].dropna()\n'
    '    groups = groups[groups[cat_col].isin(groups[cat_col].value_counts()[groups[cat_col].value_counts() >= min_n].index)]\n\n'
    '    # Order by median\n'
    '    medians = groups.groupby(cat_col)[num_col].median().sort_values()\n'
    '    groups[cat_col] = pd.Categorical(groups[cat_col], categories=medians.index, ordered=True)\n\n'
    '    fig = px.box(groups, x=cat_col, y=num_col, \n'
    '                 title=title,\n'
    '                 points="outliers",\n'
    '                 color=cat_col, color_discrete_sequence=px.colors.qualitative.Set3)\n'
    '    fig.update_layout(xaxis_tickangle=45, height=500, width=1000,\n'
    '                      yaxis_title=ylabel, showlegend=False)\n'
    '    try:\n'
    '        fig.show()\n'
    '    except:\n'
    '        pass\n\n'
    'print("\\u2705 Phase 3 helpers defined")'
))

# ── Cell 14: Test 6 — AGE × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 6: AGE × CANCER_GROUP ──\n'
    'print("=== AGE × CANCER_GROUP ===")\n'
    'h_stat, p_val, eps_sq, posthoc, groups = kw_dunn(df, \'AGE\', \'CANCER_GROUP\')\n'
    'print(f"Kruskal-Wallis: H={h_stat:.4f}, p={p_val:.6f}, \\u03b5\\u00b2={eps_sq:.4f}, groups={len(groups)}")\n'
    'if posthoc:\n'
    '    posthoc_df = pd.DataFrame(posthoc)\n'
    '    sig_posthoc = posthoc_df[posthoc_df[\'P_raw\'] < 0.05].sort_values(\'P_raw\')\n'
    '    print(f"Significant pairwise comparisons (MW, p<0.05): {len(sig_posthoc)}")\n'
    '    if len(sig_posthoc) > 0:\n'
    '        display(sig_posthoc.head(10))\n\n'
    'boxplot_by_group(df, \'AGE\', \'CANCER_GROUP\', \n'
    '                 f\'AGE distribution by CANCER_GROUP (KW p={p_val:.4f}, \\u03b5\\u00b2={eps_sq:.3f})\',\n'
    '                 \'AGE (years)\')\n\n'
    'results_ph3 = [{\n'
    '    \'Phase\': \'Phase 3\', \'Test\': \'AGE × CANCER_GROUP\',\n'
    '    \'Num_Var\': \'AGE\', \'Cat_Var\': \'CANCER_GROUP\',\n'
    '    \'N\': df[\'AGE\'].notna().sum(),\n'
    '    \'Test\': \'Kruskal-Wallis\',\n'
    '    \'Statistic\': round(h_stat, 4),\n'
    '    \'P_value\': round(p_val, 6),\n'
    '    \'Effect_Size\': round(eps_sq, 4),\n'
    '    \'Effect_Type\': \'\\u03b5\\u00b2\',\n'
    '    \'N_Groups\': len(groups)\n'
    '}]'
))

# ── Cell 15: Test 7 — TF × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 7: TF × CANCER_GROUP ──\n'
    'print("=== TF × CANCER_GROUP ===")\n'
    'h_stat, p_val, eps_sq, posthoc, groups = kw_dunn(df, \'TUMOR_FRACTION\', \'CANCER_GROUP\')\n'
    'print(f"Kruskal-Wallis: H={h_stat:.4f}, p={p_val:.6f}, \\u03b5\\u00b2={eps_sq:.4f}, groups={len(groups)}")\n'
    'if posthoc:\n'
    '    posthoc_df = pd.DataFrame(posthoc)\n'
    '    sig_posthoc = posthoc_df[posthoc_df[\'P_raw\'] < 0.05].sort_values(\'P_raw\')\n'
    '    print(f"Significant pairwise comparisons (MW, p<0.05): {len(sig_posthoc)}")\n\n'
    'boxplot_by_group(df, \'TUMOR_FRACTION\', \'CANCER_GROUP\',\n'
    '                 f\'TF distribution by CANCER_GROUP (KW p={p_val:.4f}, \\u03b5\\u00b2={eps_sq:.3f})\',\n'
    '                 \'TUMOR_FRACTION\')\n\n'
    'results_ph3.append({\n'
    '    \'Phase\': \'Phase 3\', \'Test\': \'TF × CANCER_GROUP\',\n'
    '    \'Num_Var\': \'TUMOR_FRACTION\', \'Cat_Var\': \'CANCER_GROUP\',\n'
    '    \'N\': df[\'TUMOR_FRACTION\'].notna().sum(),\n'
    '    \'Test\': \'Kruskal-Wallis\',\n'
    '    \'Statistic\': round(h_stat, 4),\n'
    '    \'P_value\': round(p_val, 6),\n'
    '    \'Effect_Size\': round(eps_sq, 4),\n'
    '    \'Effect_Type\': \'\\u03b5\\u00b2\',\n'
    '    \'N_Groups\': len(groups)\n'
    '})'
))

# ── Cell 16: Test 8 — TP × CANCER_GROUP ──
cells.append(new_code_cell(
    '# ── Test 8: TP × CANCER_GROUP ──\n'
    'print("=== TP × CANCER_GROUP ===")\n'
    'h_stat, p_val, eps_sq, posthoc, groups = kw_dunn(df, \'TUMOR_PLOIDY\', \'CANCER_GROUP\')\n'
    'print(f"Kruskal-Wallis: H={h_stat:.4f}, p={p_val:.6f}, \\u03b5\\u00b2={eps_sq:.4f}, groups={len(groups)}")\n'
    'if posthoc:\n'
    '    posthoc_df = pd.DataFrame(posthoc)\n'
    '    sig_posthoc = posthoc_df[posthoc_df[\'P_raw\'] < 0.05].sort_values(\'P_raw\')\n'
    '    print(f"Significant pairwise comparisons (MW, p<0.05): {len(sig_posthoc)}")\n\n'
    'boxplot_by_group(df, \'TUMOR_PLOIDY\', \'CANCER_GROUP\',\n'
    '                 f\'TP distribution by CANCER_GROUP (KW p={p_val:.4f}, \\u03b5\\u00b2={eps_sq:.3f})\',\n'
    '                 \'TUMOR_PLOIDY\')\n\n'
    'results_ph3.append({\n'
    '    \'Phase\': \'Phase 3\', \'Test\': \'TP × CANCER_GROUP\',\n'
    '    \'Num_Var\': \'TUMOR_PLOIDY\', \'Cat_Var\': \'CANCER_GROUP\',\n'
    '    \'N\': df[\'TUMOR_PLOIDY\'].notna().sum(),\n'
    '    \'Test\': \'Kruskal-Wallis\',\n'
    '    \'Statistic\': round(h_stat, 4),\n'
    '    \'P_value\': round(p_val, 6),\n'
    '    \'Effect_Size\': round(eps_sq, 4),\n'
    '    \'Effect_Type\': \'\\u03b5\\u00b2\',\n'
    '    \'N_Groups\': len(groups)\n'
    '})'
))

# ── Cell 17: Tests 9-11 — Numeric × Numeric correlations ──
cells.append(new_code_cell(
    '# ── Tests 9-11: Numeric × Numeric correlations ──\n\n'
    'def spearman_test(df, col1, col2, label1, label2):\n'
    '    """Spearman rank correlation with scatter plot."""\n'
    '    sub = df[[col1, col2]].dropna()\n'
    '    if len(sub) < 20:\n'
    '        print(f"  Not enough data: {len(sub)} < 20")\n'
    '        return None\n\n'
    '    rho, p_val = spearmanr(sub[col1], sub[col2])\n\n'
    '    fig = px.scatter(sub, x=col1, y=col2, \n'
    '                     title=f\'{label1} vs {label2} (\\u03c1={rho:.4f}, p={p_val:.4f})\',\n'
    '                     opacity=0.6, trendline=\'lowess\')\n'
    '    fig.update_layout(height=500, width=600)\n'
    '    try:\n'
    '        fig.show()\n'
    '    except:\n'
    '        pass\n\n'
    '    print(f"  Spearman \\u03c1={rho:.4f}, p={p_val:.6f}, N={len(sub)}")\n'
    '    return {\n'
    '        \'Phase\': \'Phase 3\',\n'
    '        \'Test\': f\'{label1} \\u00d7 {label2}\',\n'
    '        \'Num_Var\': label1,\n'
    '        \'Cat_Var\': label2,\n'
    '        \'N\': len(sub),\n'
    '        \'Test\': \'Spearman\',\n'
    '        \'Statistic\': round(rho, 4),\n'
    '        \'P_value\': round(p_val, 6),\n'
    '        \'Effect_Size\': round(rho, 4),\n'
    '        \'Effect_Type\': \'\\u03c1\',\n'
    '        \'N_Groups\': 2\n'
    '    }\n\n'
    'print("=== Numeric × Numeric Correlations ===")\n\n'
    '# Test 9: AGE × TF\n'
    'print("Test 9: AGE × TF")\n'
    'res9 = spearman_test(df, \'AGE\', \'TUMOR_FRACTION\', \'AGE\', \'TF\')\n'
    'if res9: results_ph3.append(res9)\n\n'
    '# Test 10: AGE × TP\n'
    'print("Test 10: AGE × TP")\n'
    'res10 = spearman_test(df, \'AGE\', \'TUMOR_PLOIDY\', \'AGE\', \'TP\')\n'
    'if res10: results_ph3.append(res10)\n\n'
    '# Test 11: TF × TP\n'
    'print("Test 11: TF × TP")\n'
    'res11 = spearman_test(df, \'TUMOR_FRACTION\', \'TUMOR_PLOIDY\', \'TF\', \'TP\')\n'
    'if res11: results_ph3.append(res11)'
))

# ── Cell 18: Phase 3 results table ──
cells.append(new_code_cell(
    '# ── Phase 3 Results ──\n'
    'df_ph3 = pd.DataFrame(results_ph3)\n'
    'cols = [\'Phase\', \'Test\', \'Num_Var\', \'Cat_Var\', \'N\', \'Test\', \'Statistic\', \'P_value\', \'Effect_Size\', \'Effect_Type\', \'N_Groups\']\n'
    '# Remove duplicate \'Test\' column issue\n'
    'df_ph3 = df_ph3.loc[:, ~df_ph3.columns.duplicated()]\n'
    'display(df_ph3.style\n'
    '    .format({\'P_value\': \'{:.6f}\', \'Statistic\': \'{:.4f}\', \'Effect_Size\': \'{:.4f}\'})\n'
    '    .background_gradient(cmap=\'RdYlGn_r\', subset=[\'P_value\'])\n'
    ')'
))

# ── Cell 19: Combined results + save ──
cells.append(new_code_cell(
    '# ── Combine Phase 2 + Phase 3 results ──\n'
    '# Phase 2 results\n'
    'df2 = pd.DataFrame(results_ph2)\n'
    'df2[\'Phase\'] = \'Phase 2\'\n'
    'df2 = df2.rename(columns={\n'
    '    \'Chi2\': \'Statistic_Val\', \'Cramers_V\': \'Effect_Size\'\n'
    '})\n'
    '# Already has P_value\n\n'
    '# Phase 3 results\n'
    'df3 = pd.DataFrame(results_ph3)\n'
    'df3[\'Phase\'] = \'Phase 3\'\n'
    '# Already has Statistic, P_value, Effect_Size\n\n'
    '# Combine\n'
    'all_cols = [\'Phase\', \'Test\', \'N\', \'P_value\', \'Effect_Size\']\n'
    'combined = pd.concat([df2, df3], ignore_index=True)\n\n'
    '# FDR correction within Phase\n'
    'for phase in combined[\'Phase\'].unique():\n'
    '    mask = combined[\'Phase\'] == phase\n'
    '    pvals = combined.loc[mask, \'P_value\'].values\n'
    '    reject, qvals, _, _ = smm.multipletests(pvals, method=\'fdr_bh\')\n'
    '    combined.loc[mask, \'FDR_q\'] = qvals\n\n'
    'display(combined.style\n'
    '    .format({\'P_value\': \'{:.6f}\', \'Effect_Size\': \'{:.4f}\', \'FDR_q\': \'{:.4f}\'})\n'
    '    .background_gradient(cmap=\'RdYlGn_r\', subset=[\'P_value\'])\n'
    ')\n\n'
    '# ── Save ──\n'
    'out_dir = \'/home/alon/menow_home_ass/notebooks/clinical_associations\'\n'
    'os.makedirs(out_dir, exist_ok=True)\n'
    'combined.to_csv(f\'{out_dir}/clinical_associations_results.csv\', index=False)\n\n'
    '# Significant summary\n'
    'sig = combined[combined[\'FDR_q\'] < 0.05].copy()\n'
    'sig.to_csv(f\'{out_dir}/clinical_associations_summary.csv\', index=False)\n\n'
    'print(f"\\n\\u2705 Results saved: {len(combined)} total tests, {len(sig)} FDR-significant (q<0.05)")'
))

# ── Cell 20: Summary markdown ──
cells.append(new_markdown_cell(
    '## Summary & Next Steps\n\n'
    'This analysis has characterized the relationships between clinical variables '
    'in the PBTA_RNA cohort, independent of outcome.\n\n'
    '**Phase 2** identified associations between categorical clinical variables.\n\n'
    '**Phase 3** identified differences in numeric variables across cancer groups '
    'and correlations between numeric pairs.\n\n'
    '### Next Steps\n'
    '- **Phase 4**: Multivariate models (Cox PH) incorporating all variables\n'
    '- **Phase 5**: Unsupervised subgroup discovery\n'
    '- **Phase 6**: Consolidated summary across all phases'
))

# ── Cell 21: Placeholder for detailed findings ──
cells.append(new_markdown_cell(
    '### Detailed Findings\n\n'
    '*(This section to be populated after reviewing the output CSV files.)*\n\n'
    'Key observations:\n'
    '- Significant associations from Phase 2 indicate which categorical variables are related\n'
    '- Significant KW results from Phase 3 show which numeric variables differ by cancer group\n'
    '- Spearman correlations reveal monotonic relationships between AGE, TF, and TP\n\n'
    'See `clinical_associations_results.csv` for the complete table of all tests '
    'with p-values and effect sizes.'
))

nb.cells = cells

with open(NB_PATH, 'w') as f:
    nbformat.write(nb, f)

print(f"✅ Notebook saved to {NB_PATH}")
print(f"   {len(cells)} cells created")
