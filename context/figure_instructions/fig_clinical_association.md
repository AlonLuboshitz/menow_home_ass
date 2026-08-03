# Figure 3: Cross-dataset clinical associations (SEX, predisposition, race, AGE/TF/TP per group, correlations)

## Figure-level info
- Report subsection: "Cross datasets analysis"
- LaTeX label / reference: referenced in the report text as literal `Figure 3A` … `Figure 3H` (no `\ref`); the figure environment has `\label{fig:clinical_association}`. Resolves to "Figure 3" in report order.
- Dummy file to replace: `Figures/dummy_3.png`
- Output file: `Figures/fig_clinical_association.png` (300 dpi) [+ `.pdf`]
- Layout: 2-column grid, 8 panels (A–H). Reading order left→right, top→bottom: A, B, C, D, E, F, G, H. Wide stacked-bar panels (A, B, C) may each span both columns; D, E, F share rows; G and H (3 scatter subplots each) span both columns on the last two rows.
- Notebooks relied upon:
  - `notebooks/clinical_associations_analysis_executed.ipynb` — Phase 2 (Test 1 SEX, Test 2 predispositions, Test 4 race) and Phase 3 (Test 5 AGE, Test 6 TF, Test 7 TP, Tests 8–10 correlations)
  - src: `notebooks/clinical_associations/src/build_nb.py` (cells for Tests 1–10)
  - saved results: `notebooks/clinical_associations/clinical_associations_results.csv` (has `FDR_WithinFamily`, `N`, `N_events`, `Effect_Size`) and `clinical_associations_summary.csv`
- General styling (apply to all panels, per general instruction §4):
  - Panel letter: **bold, uppercase, top-left corner of the axes area** (A–H).
  - Significance for ALL panels: **FDR Benjamini-Hochberg, `FDR_WithinFamily` from the results CSV, q < 0.05** (user decision 2026-08-02). Star notation per `context/clinical_deep_dive_general.md` §3.6 applied to q (`*` q<0.05, `**` q<0.01, `***` q<0.001); groups with q≥0.05 are "ns". Always annotate raw p and FDR q per group.
  - All tests are **sample-level** (merged patient×sample table; 4074 rows after excluding samples with missing CANCER_GROUP; 2701 unique patients). Per-group N = sample count, except SEX (non-missing SEX) and predispositions (samples with a known predisposition).
  - Informative axis labels with units; legend whenever >1 series; "Unknown"/missing as explicit category where the notebook does; colors per §4.2 (`px.colors.qualitative.Set3` / `Plotly`).
  - No global figure title needed (caption lives in the .tex), but give each panel an informative title.

## Panel overview
| Panel | Title | Plot type | Notebook step | N | Statistics to annotate |
|-------|-------|-----------|---------------|---|------------------------|
| A | SEX composition per cancer group (binomial vs 50:50) | Horizontal stacked %Male/%Female bar, 50% vline, significance markers | clinical_associations Test 1 (build_nb.py cell 7) | 4051 samples with SEX (26 groups, n≥20) | Binomial p & FDR q per group; 6 FDR-sig CGs |
| B | Predisposition composition per cancer group | Stacked bar (top-15 categories + Other) | Test 2 (cell 9) | 368 samples with pred (9.0% of 4074); 5 CGs tested | 5/5 FDR-sig; p, q, Cramer's V per group in legend |
| C | Race composition per group vs overall cohort | Stacked bar per CG, groups sorted by composition similarity | Test 4 (cell 13) | 4074 samples (26 groups, n≥20) | GoF p & FDR q per group; 9 FDR-sig CGs; q written vertically on the bars; overall race distribution |
| D | AGE distribution by cancer group | Box + outlier points | Test 5 (cell 18) | 3800 samples with AGE | KW H=826.74, p<0.001, ε²=0.2124; 14 FDR-sig CGs |
| E | TF distribution by cancer group | Box + outlier points | Test 6 (cell 20) | 2538 samples with TF | KW H=100.09, p<0.001, ε²=0.0318; 6 FDR-sig CGs |
| F | TP distribution by cancer group | Box + outlier points | Test 7 (cell 22) | 2591 samples with TP | KW H=54.50, p=4.9e-05, ε²=0.0134; 4 FDR-sig CGs |
| G | Spearman correlations (AGE×TF, AGE×TP, TF×TP) | 3 scatter subplots + LOESS trendline | Tests 8–10 (cell 24) | 2694 / 2751 / 2744 | Spearman ρ, p, N per pair (all significant, |ρ|<0.1) |
| H | Pearson correlations (AGE×TF, AGE×TP, TF×TP) | 3 scatter subplots + LOESS trendline | ad-hoc from merged table (same pairwise dropna) | 2694 / 2751 / 2744 | Pearson r, p, N per pair (2/3 significant, |r|<0.1) |

## Panel A: SEX composition per cancer group (binomial vs 50:50)
- **Shows:** For 6 of 26 cancer groups the male:female split deviates significantly from the 50:50 baseline (FDR q<0.05).
- **Data:** merged patient×sample table; `SEX` filtered to `Female`/`Male` (overall 2192 M / 1859 F = 4051 samples, 45.9%F). Per-group two-sided `binom_test(n_f, n, 0.5)`; groups with n≥20 retained (26 groups, sample-level).
- **Plot type:** Horizontal stacked bar, `%Male` (steelblue) + `%Female` (lightcoral), `barmode="overlay"`, dashed vertical line at 50%, significance markers in a side column — mirrors the Test 1 plot (build_nb.py cell 7). Groups on the y-axis sorted by effect size (difference from 50%).
- **X axis:** Percent (0–100; significance marker column to ~105).
- **Y axis:** Cancer group.
- **Grouping/colors:** 2 series (Male/Female); emphasize the 6 FDR-significant groups (bold / marker), mute the "ns" groups.
- **Statistics to annotate:** binomial p and FDR q per group; star notation from q.
- **N to show:** N per group = samples with SEX in that group (sample-level, NOT patient-level); total N = 4051.
- **Exact values that must appear (FDR-significant, q<0.05):** Medulloblastoma **40.7%F** q=0.0028; Ewing sarcoma **26.1%F** q=0.0142; Adamantinomatous Craniopharyngioma **34.8%F** q=0.0142; Ependymoma **41.7%F** q=0.0177; Low-grade glioma **45.1%F** q=0.0217; Ganglioglioma **38.4%F** q=0.0238. Fall below the FDR threshold (raw p<0.05 but ns): Chordoma 72.4%F q=0.0896; Diffuse midline glioma 55.1%F q=0.1373. Overall cohort 45.9%F.
- **Special annotations:** Bold uppercase "A" top-left. 50% dashed reference line. Annotate "%F / q" on or beside the significant bars.

## Panel B: Predisposition composition per cancer group
- **Shows:** 9.0% of samples (368/4074) carry a documented cancer predisposition across 23 categories; 5 cancer groups (n≥20 predisposed samples) each have a predisposition composition that differs from all other CGs combined (5/5 FDR-significant).
- **Data:** merged table; `CANCER_PREDISPOSITIONS` cleaned via `clean_pred()`; exclude `No predisposition` / `Unknown` → 368 records, 23 categories. Per-group chi-squared (Fisher exact if expected <5) vs all-other-CGs combined; groups with n≥20 retained (5 groups). Cramer's V as effect size.
- **Plot type:** Stacked bar chart — % of predisposed samples per category per CG (mirrors the Test 2 plot: top-15 categories by prevalence + pooled "Other").
- **X axis:** Cancer group.
- **Y axis:** % of predisposed samples (within-CG).
- **Grouping/colors:** one stacked series per category (legend); "Other" pools the remaining categories.
- **Statistics to annotate:** per-group p, FDR q, Cramer's V; all 5 groups significant (q<0.05).
- **N to show:** N = 368 predisposed samples (9.0% of 4074). Per-test N (predisposed samples in the CG): Schwannoma 22, ATRT 29, HGG 80, Neurofibroma/Plexiform 43, LGG 72.
- **Exact values that must appear:** Top predisposition categories (sample counts): NF-1 **130**, Other inherited conditions NOS **64**, Li-Fraumeni (TP53) **58**, NF-2 **40**, Tuberous Sclerosis **14**, Lynch **10**, VHL **8**, Rhabdoid predisposition (SMARCB1/SMARCA4) **7**, CMMRD **6**, Li-Fraumeni+NF-1 **6**, all others ≤3. Per-group tests: Schwannoma N=22 p=0 q=0 V=0.6176; ATRT N=29 p=0 q=0 V=0.5592; HGG N=80 p=0 q=0 V=0.4736; Neurofibroma/Plexiform N=43 p=0 q=0 V=0.3890; LGG N=72 p=1.1e-05 q=1.1e-05 V=0.3304. The report's "Y CGs" = 5.
- **Special annotations:** Bold uppercase "B" top-left. Annotate "9.0% of samples (368/4074) have a documented predisposition". Note per-CG top categories (e.g. LGG NF-1 48/72; Neurofibroma/Plexiform NF-1 41/43; Schwannoma NF-2 16/22; HGG Li-Fraumeni 21; ATRT Other inherited conditions NOS 12). Per-group test stats (p, q, V) are shown in the legend (label-only entries below the category swatches), NOT as text above the bars — user decision 2026-08-02 (the above-bar text overlapped the tall stacked bars).

## Panel C: Race composition per group vs overall cohort
- **Shows:** The racial composition of 9 of 26 cancer groups differs from the overall cohort's race distribution (FDR q<0.05); groups with similar composition are placed next to each other.
- **Data:** merged table; `RACE_GROUP` (White / Black / Asian / Other / Unknown, via the notebook's `race_group_map`); per-group chi-squared goodness-of-fit vs the overall proportions, groups with n≥20 retained (26 groups, sample-level). Overall N=4074.
- **Plot type:** Stacked bar per cancer group (5 race categories); groups ordered by composition similarity, NOT alphabetical (user decision) — hierarchical clustering ordering (e.g. `scipy.cluster.hierarchy.linkage(per-CG race proportion vectors, method="complete", metric="euclidean")`). Overall-cohort composition shown as a reference.
- **X axis:** Cancer group (clustering order).
- **Y axis:** % of samples.
- **Grouping/colors:** 5 stacked categories (White/Black/Asian/Other/Unknown); overall cohort reference distinguishable (e.g. black bar outline).
- **Statistics to annotate:** per-group GoF p and FDR q; mark the 9 FDR-significant CGs.
- **N to show:** per-CG sample counts; overall N = 4074 samples.
- **Exact values that must appear:** Overall cohort: White 60.8% (2475), Black 8.0% (327), Asian 3.4% (140), Other 1.7% (71), Unknown 26.0% (1061). 9 FDR-significant CGs — q and within-group W/B/A/O/U%: Diffuse midline glioma q=0 (45.8/5.2/5.5/2.1/41.3); Chordoma q=0.0008 (34.5/0.0/0.0/0.0/65.5); Low-grade glioma q=0.0016 (67.7/7.9/2.1/1.5/20.8); Meningioma q=0.0016 (48.0/20.0/5.0/1.0/26.0); Neurofibroma/Plexiform q=0.0016 (58.3/25.0/0.0/2.1/14.6); Schwannoma q=0.0036 (77.9/4.4/2.9/5.9/8.8); ATRT q=0.0195 (51.3/15.8/3.9/2.6/26.3); Ganglioglioma q=0.0195 (74.2/6.6/4.0/2.0/13.2); Ewing sarcoma q=0.0374 (84.8/0.0/0.0/2.2/13.0). Fall below FDR: Oligodendroglioma q=0.0749; Medulloblastoma q=0.0883.
- **Special annotations:** Bold uppercase "C" top-left. Per user decision 2026-08-02 the report's "sort the groups that has similar composition" instruction is implemented via hierarchical clustering on the per-CG race proportion vectors (specified above); do NOT use alphabetical order. q-values of FDR-significant groups are written VERTICALLY on their bars (rotation 90°, centered in each bar), NOT horizontally above them — user decision 2026-08-02 (the above-bar text overlapped the tall stacked bars).

## Panel D: AGE distribution by cancer group
- **Shows:** AGE at diagnosis differs across cancer groups; 14 of the 26 tested CGs differ significantly from all others combined (FDR q<0.05).
- **Data:** merged table, `AGE` dropna (N=3800 samples); global Kruskal-Wallis + per-group Mann-Whitney (each CG vs all others combined), groups with n≥20.
- **Plot type:** Box plot with jittered outlier points, x-axis sorted by median AGE (`px.box(..., points="outliers", color=cat_col)`, `Set3`) — mirrors Test 5 (cell 18).
- **X axis:** Cancer group (by median AGE); tick angle −45°.
- **Y axis:** AGE (years).
- **Grouping/colors:** one box per CG (Set3); significant CGs emphasized.
- **Statistics to annotate:** KW in title: H=826.74, p<0.001, ε²=0.2124; per-group MW p & q; mark the 14 FDR-significant CGs.
- **N to show:** N = 3800 samples with AGE (of 4074).
- **Exact values that must appear:** 14 FDR-significant CGs (q<0.05): Embryonal tumor with multilayer rosettes, High-grade glioma, Meningioma, Oligodendroglioma, Ependymoma, Glial-neuronal tumor NOS, Choroid plexus tumor, Neurofibroma/Plexiform, Schwannoma, Atypical Teratoid Rhabdoid Tumor, Dysembryoplastic neuroepithelial tumor, Hemangioblastoma, Diffuse hemispheric glioma, Neuroblastoma. 18 raw-significant (p<0.05); 4 drop after FDR: CNS Embryonal tumor q=0.0539, Chordoma q=0.0519, Pilocytic astrocytoma q=0.0680, Sarcoma q=0.0507.
- **Special annotations:** Bold uppercase "D" top-left.

## Panel E: TF distribution by cancer group
- **Shows:** Tumor fraction differs across cancer groups; 6 of the tested CGs differ significantly from all others (FDR q<0.05).
- **Data:** merged table, `TUMOR_FRACTION` dropna (N=2538 samples); KW + per-group MW, groups with n≥20.
- **Plot type:** Box plot with jittered outlier points, x sorted by median TF — mirrors Test 6 (cell 20).
- **X axis:** Cancer group (by median TF); tick angle −45°.
- **Y axis:** Tumor fraction (0–1).
- **Grouping/colors:** one box per CG (Set3); significant CGs emphasized.
- **Statistics to annotate:** KW in title: H=100.09, p<0.001, ε²=0.0318; per-group MW p & q; mark the 6 FDR-significant CGs.
- **N to show:** N = 2538 samples with TF.
- **Exact values that must appear:** 6 FDR-significant CGs (q<0.05): High-grade glioma q=2.2e-05, Diffuse midline glioma q=0, Ependymoma q=0.0005, Adamantinomatous Craniopharyngioma q=0.0002, Low-grade glioma q=0.0414, Neurofibroma/Plexiform q=0.0123. 7 raw-significant; 1 drops after FDR: Schwannoma q=0.1203.
- **Special annotations:** Bold uppercase "E" top-left.

## Panel F: TP distribution by cancer group
- **Shows:** Tumor ploidy differs across cancer groups; 4 of the tested CGs differ significantly from all others (FDR q<0.05).
- **Data:** merged table, `TUMOR_PLOIDY` dropna (N=2591 samples); KW + per-group MW, groups with n≥20.
- **Plot type:** Box plot with jittered outlier points, x sorted by median TP — mirrors Test 7 (cell 22).
- **X axis:** Cancer group (by median TP); tick angle −45°.
- **Y axis:** Tumor ploidy (2–4).
- **Grouping/colors:** one box per CG (Set3); significant CGs emphasized.
- **Statistics to annotate:** KW in title: H=54.50, p=4.9e-05, ε²=0.0134; per-group MW p & q; mark the 4 FDR-significant CGs.
- **N to show:** N = 2591 samples with TP.
- **Exact values that must appear:** 4 FDR-significant CGs (q<0.05): Choroid plexus tumor q=0.0265, Chordoma q=0.0265, Low-grade glioma q=0.0488, Dysembryoplastic neuroepithelial tumor q=0.0488. 6 raw-significant; 2 drop after FDR: Meningioma q=0.1425, Diffuse midline glioma q=0.1425.
- **Special annotations:** Bold uppercase "F" top-left.

## Panel G: Spearman correlations (AGE×TF, AGE×TP, TF×TP)
- **Shows:** The numeric variables are significantly but negligibly correlated (all Spearman p<0.05 with |ρ|<0.1). Significance is nominal — driven by the large N (≈2700); the |ρ|<0.1 effects are not practically meaningful.
- **Data:** merged table, pairwise dropna per pair; Spearman rank correlation (scipy `spearmanr`). Pearson r is computed ad-hoc and shown in Panel H.
- **Plot type:** 3 scatter subplots (AGE×TF, AGE×TP, TF×TP) with LOESS trendline in red (`px.scatter(..., trendline="lowess", trendline_color_override="red", opacity=0.5)`) — mirrors Tests 8–10 (cell 24).
- **X axis:** AGE (years) / TF (0–1) / TF (0–1) per subplot.
- **Y axis:** TF (0–1) / TP (2–4) / TP (2–4) per subplot.
- **Grouping/colors:** single point series per subplot + red LOESS line.
- **Statistics to annotate:** Spearman ρ, p, N in each subplot title/annotation.
- **N to show:** AGE×TF N=2694, AGE×TP N=2751, TF×TP N=2744 (samples with both values).
- **Exact values that must appear:** AGE×TF ρ=−0.098 (−0.0978), p<0.001, N=2694; AGE×TP ρ=+0.044 (0.0437), p=0.0219, N=2751; TF×TP ρ=−0.086 (−0.0855), p=7e-06, N=2744.
- **Special annotations:** Bold uppercase "G" top-left. Presented as significant p-values with no observed trends (|ρ|<0.1), possibly due to the high sample amount (user decision 2026-08-02).

## Panel H: Pearson correlations (AGE×TF, AGE×TP, TF×TP)
- **Shows:** Linear correlations between the same numeric variables as Panel G. Significant p-values with no observed trends (|r|<0.1), possibly due to the high sample amount. 2 of 3 pairs are significant (AGE×TF, TF×TP); AGE×TP is ns (Pearson p=0.0986, unlike its significant Spearman p=0.0219 in Panel G). Significance is nominal — driven by the large N (≈2700); the |r|<0.1 effects are not practically meaningful.
- **Data:** merged table, pairwise dropna per pair (identical to Panel G, so N matches); Pearson correlation via scipy `pearsonr`. This is the "pearson" of the report's "pearson/spearman". Ad-hoc computation (not in the analysis notebook, so no FDR CSV check).
- **Plot type:** 3 scatter subplots (AGE×TF, AGE×TP, TF×TP) with red LOESS trendline — same style as Panel G, so both panels look consistent.
- **X axis:** AGE (years) / TF (0–1) / TF (0–1) per subplot.
- **Y axis:** TF (0–1) / TP (2–4) / TP (2–4) per subplot.
- **Grouping/colors:** single point series per subplot + red LOESS line.
- **Statistics to annotate:** Pearson r, p, N in each subplot title/annotation; mark AGE×TP as ns.
- **N to show:** AGE×TF N=2694, AGE×TP N=2751, TF×TP N=2744 (samples with both values; same as Panel G).
- **Exact values that must appear:** AGE×TF r=−0.053 (−0.0527), p=0.0063, N=2694; AGE×TP r=+0.032 (0.0315), p=0.0986 (ns), N=2751; TF×TP r=−0.065 (−0.0645), p=0.0007, N=2744.
- **Special annotations:** Bold uppercase "H" top-left. Same conclusion wording as G: significant p-values with no observed trends (|r|<0.1), possibly due to the high sample amount; annotate AGE×TP as not significant (ns).

## Open questions / flags
- **Significance convention (FDR vs raw) — RESOLVED (2026-08-02):** All panels annotate FDR q (`FDR_WithinFamily`, Benjamini-Hochberg, results CSV), significant at q<0.05. The report previously quoted raw-count results; the user chose FDR. Counts: SEX **6** FDR-sig (Medulloblastoma q=0.0028, Ewing sarcoma q=0.0142, Adamantinomatous Craniopharyngioma q=0.0142, Ependymoma q=0.0177, Low-grade glioma q=0.0217, Ganglioglioma q=0.0238; Chordoma q=0.0896 and Diffuse midline glioma q=0.1373 drop); Race **9** FDR-sig (Diffuse midline glioma q=0, Chordoma q=0.0008, Low-grade glioma q=0.0016, Meningioma q=0.0016, Neurofibroma/Plexiform q=0.0016, Schwannoma q=0.0036, ATRT q=0.0195, Ganglioglioma q=0.0195, Ewing sarcoma q=0.0374; Oligodendroglioma q=0.0749 and Medulloblastoma q=0.0883 drop); AGE **14** FDR-sig of 18 raw (drops: CNS Embryonal tumor q=0.0539, Chordoma q=0.0519, Pilocytic astrocytoma q=0.0680, Sarcoma q=0.0507); TF **6** FDR-sig of 7 raw (Schwannoma q=0.1203 drops); TP **4** FDR-sig of 6 raw (Meningioma q=0.1425, Diffuse midline glioma q=0.1425 drop).
- **Predisposition prevalence "~4%" vs 9.0% — RESOLVED (2026-08-02):** Report text said "~9%" in one sentence and the figure brief previously said "~4%". The executed value is **9.0%** (368/4074 samples, 23 categories). Panel B uses 9.0% (368/4074).
- **Panel C group ordering — RESOLVED (2026-08-02):** Report instruction "sort the groups that has similar composition" is implemented as hierarchical clustering ordering (scipy `linkage`, method="complete", Euclidean metric) on the per-CG race proportion vectors — NOT alphabetical (user decision).
- **Correlation panels G/H renumbering — RESOLVED (2026-08-02):** Correlations are split into **Panel G (Spearman)** and **Panel H (Pearson)** (user decision 2026-08-02); the report text references "Figure 3G, H". There is no subtype-heatmap panel in Figure 3.
- **Correlations significance conclusion — RESOLVED (2026-08-02):** The report text now reads "did yield significant p values with no observed trends possibly due to high sample amount". Both panels annotate the significant p-values and note the negligible effects (|r|,|ρ|<0.1) and the large N (≈2700). Exception: Panel H AGE×TP Pearson is ns (p=0.0986) — flagged below.
- **Correlations "pearson/spearman" — RESOLVED (2026-08-02):** The report says "pearson/spearman"; Pearson is now computed and shown in **Panel H**, Spearman in **Panel G** (user decision 2026-08-02). Pearson is ad-hoc (not in the analysis notebook); Spearman reproduces Tests 8–10 (cell 24).
- **Panel H AGE×TP Pearson not significant — FLAG (2026-08-02):** AGE×TP Pearson r=+0.032 (0.0315), p=0.0986 is **not** significant (p≥0.05), unlike its Spearman p=0.0219 (Panel G). Panel H annotates it as "ns". The conclusion wording (significant p-values, no observed trends, possibly due to high sample amount) therefore applies to 2 of the 3 Pearson pairs; AGE×TP shows a weak positive linear trend that is not significant.
- **Notebook comment name — flag (no decision needed):** The figure environment comment says "% relay on clinical_association_excuted"; the actual notebook is `clinical_associations_analysis_executed.ipynb`. Mapping by content (general instruction §8).
- **Dummy file name — flag (no decision needed):** The report references `Figures/dummy_3.png`, which is not present in `Figures/` (only `dummy_1.png` and misspelled `dummpy_2.png`). Rendering agent should write `Figures/fig_clinical_association.png` regardless.
- **Panel B/C statistic placement — RESOLVED (2026-08-02):** Panel B stats moved into the legend; Panel C q-values written vertically on the bars. Both previously sat at y=102 above the stacked bars and overlapped them.
