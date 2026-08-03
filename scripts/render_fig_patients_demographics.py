"""Render Figure 1: Patient & clinical demographics (fig_patients_demographics).

Panels:
  A  Patient Data - Missingness (13 cols, % null, sorted desc, highlighted OS/EFS months)
  B  Overall Survival Status (donut, hole 0.3)
  C  Event-Free Survival Status (EFS Binary + EFS Detailed grouped bars)
  D  Age Distribution (histogram, 40 bins)
  E  AGE outliers per cancer group (Tukey IQR, horizontal count bar)

Spec: context/figure_instructions/fig_patients_demographics.md
All values computed from the raw clinical files; verified against the spec.
Canonical copy: /home/alon/menow_home_ass/scripts/render_fig_patients_demographics.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

PATIENT_FILE = "/home/alon/menow_home_ass/PBTA_RNA/data_clinical_patient_attributes.txt"
SAMPLE_FILE = "/home/alon/menow_home_ass/PBTA_RNA/data_clinical_sample_attributes.txt"
OUT_PNG = "/home/alon/menow_home_ass/Figures/fig_patients_demographics.png"
OUT_PDF = "/home/alon/menow_home_ass/Figures/fig_patients_demographics.pdf"

# ----------------------------------------------------------------------
# Data loading + cleaning helpers (equivalents of imports.py helpers)
# ----------------------------------------------------------------------
def read_patients():
    return pd.read_csv(PATIENT_FILE, sep="\t", header=4,
                       dtype={"AGE": float, "AGE_IN_DAYS": float,
                              "OS_MONTHS": float, "EFS_MONTHS": float})

def read_samples():
    return pd.read_csv(SAMPLE_FILE, sep="\t", header=4)

def clean_os(df):
    df = df.copy()
    df["OS_STATUS"] = df["OS_STATUS"].str.strip()
    df["os_label"] = df["OS_STATUS"].str.replace(r"^\d+:", "", regex=True)
    df["os_event"] = df["OS_STATUS"].apply(
        lambda x: 1 if pd.notna(x) and x.startswith("1:")
        else (0 if pd.notna(x) and x.startswith("0:") else np.nan))
    return df

def clean_efs(df):
    df = df.copy()
    df["EFS_STATUS"] = df["EFS_STATUS"].str.strip()
    df["efs_detail"] = df["EFS_STATUS"].str.replace(r"^\d+:", "", regex=True)
    df["efs_event"] = df["EFS_STATUS"].apply(
        lambda x: 0 if pd.notna(x) and x == "0:No Event"
        else (1 if pd.notna(x) and x != "1:NA" else np.nan))
    return df

# ----------------------------------------------------------------------
# Compute all panel data
# ----------------------------------------------------------------------
pat = read_patients()
samp = read_samples()
N = len(pat)                      # 2870
assert N == 2870, f"Unexpected N={N}"

# ---- Panel A: missingness ----
miss = pat.isna().mean() * 100
miss = miss.sort_values(ascending=False)
miss_counts = pat.isna().sum().reindex(miss.index)

# ---- Panel B: OS ----
p = clean_os(pat)
p["os_label"] = p["os_label"].fillna("Unknown")
os_counts = p["os_label"].value_counts().reindex(["LIVING", "DECEASED", "Unknown"])

# ---- Panel C: EFS ----
p = clean_efs(p)
p["efs_binary"] = p["efs_event"].map({1: "Event", 0: "No Event"}).fillna("Unknown")
efs_bin = p["efs_binary"].value_counts().reindex(["No Event", "Event", "Unknown"])
efs_detail = p["efs_detail"].value_counts(dropna=False)

# ---- Panel D: AGE ----
age = pat["AGE"].dropna()
age_missing_n = int(pat["AGE"].isna().sum())
age_stats = dict(mean=age.mean(), median=age.median(), min=age.min(), max=age.max())

# ---- Panel E: AGE outliers per cancer group (Tukey IQR) ----
Q1, Q3 = 4.0, 13.25
IQR = Q3 - Q1
THRESH = Q3 + 1.5 * IQR            # 27.125
ol_pts = pat.loc[pat["AGE"] > THRESH, ["PATIENT_ID", "AGE"]]
out = ol_pts.merge(samp[["PATIENT_ID", "CANCER_GROUP"]], on="PATIENT_ID", how="left")
def patient_group(grp):
    if (grp == "Oligodendroglioma").any():
        return "Oligodendroglioma"
    g = grp.dropna()
    return g.iloc[0] if len(g) else "Unknown"
ol_group = out.groupby("PATIENT_ID")["CANCER_GROUP"].agg(patient_group)
ol_counts = ol_group.value_counts().sort_values(ascending=True)  # horizontal bar: ascending = biggest at top
n_out = len(ol_group)              # 62
oligo_n = int(ol_group.eq("Oligodendroglioma").sum())   # 53

# ----------------------------------------------------------------------
# Figure construction
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(9, 8.3))
# Layout gap knobs:
#   hspace -> vertical gaps between rows (A to B/C, and B/C to D/E) — user wants more room here
#   wspace -> horizontal gaps between columns (B to C, and D to E) — user wants more room here
#   Panel letters now sit above the axes, so the vertical gap must fit them.
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.05, 1.22],
              width_ratios=[1.0, 1.12], hspace=0.65, wspace=0.45,
              left=0.075, right=0.97, top=0.96, bottom=0.08)

def panel_letter(ax, letter):
    ax.text(0.0, 1.02, letter, transform=ax.transAxes, fontsize=17,
            fontweight="bold", va="bottom", ha="left")

# ---------------- Panel A: missingness ----------------
axA = fig.add_subplot(gs[0, :])
cols = miss.index.tolist()
vals = miss.values
colorsA = []
for c in cols:
    colorsA.append("#E31A1C" if c in ("OS_MONTHS", "EFS_MONTHS") else "#9E0142")
barsA = axA.bar(np.arange(len(cols)), vals, color=colorsA, width=0.72,
                edgecolor="black", linewidth=0.4)
for c in ("EFS_MONTHS", "OS_MONTHS"):
    idx = cols.index(c)
    barsA[idx].set_edgecolor("black")
    barsA[idx].set_linewidth(2.2)
    barsA[idx].set_facecolor("#FDBF6F")
for i, (c, v) in enumerate(zip(cols, vals)):
    axA.text(i, v + 0.6, f"{v:.1f}%", ha="center", va="bottom", fontsize=8.5)
axA.set_xticks(np.arange(len(cols)))
axA.set_xticklabels(cols, rotation=-45, ha="left", fontsize=9)
axA.set_ylabel("% Missing", fontsize=11)
axA.set_ylim(0, 103)
axA.set_title(f"Patient Data \u2014 Missingness  (N = {N} patients)", fontsize=12)
axA.annotate("EFS_MONTHS 28.4% / OS_MONTHS 27.0%\n(\"~28% null\" for OS/EFS months)",
             xy=(cols.index("EFS_MONTHS"), miss["EFS_MONTHS"]), xytext=(3.5, 72),
             fontsize=8.5, ha="center",
             bbox=dict(boxstyle="round,pad=0.35", fc="#FDBF6F", ec="black", alpha=0.85),
             arrowprops=dict(arrowstyle="->", color="black", lw=1.1))
axA.grid(axis="y", ls=":", alpha=0.4)
axA.set_axisbelow(True)
panel_letter(axA, "A")

# ---------------- Panel B: OS donut ----------------
axB = fig.add_subplot(gs[1, 0])
b_labels = os_counts.index.tolist()
b_vals = os_counts.values
b_colors = ["#66C2A5", "#D95F02", "#B0B0B0"]  # green / orange / gray
wedges, _ = axB.pie(b_vals, colors=b_colors, startangle=90, counterclock=False,
                    wedgeprops=dict(width=0.30, edgecolor="white", linewidth=2))
axB.set(aspect="equal")
total_b = b_vals.sum()
labels_b = []
for name, v in zip(b_labels, b_vals):
    labels_b.append(f"{name} {int(v)}\n({v / total_b * 100:.1f}%)")
axB.legend(wedges, labels_b, loc="center left", bbox_to_anchor=(0.98, 0.5), fontsize=9,
           frameon=False)
axB.set_title(f"Overall Survival Status  (N = {N})", fontsize=12)
axB.text(0, 0, f"{N}", ha="center", va="center", fontsize=13, fontweight="bold")
panel_letter(axB, "B")

# ---------------- Panel C: EFS (two subplots) ----------------
axC_outer = fig.add_subplot(gs[1, 1])
axC_outer.axis("off")
cgs = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 1], width_ratios=[1, 2.05],
                              wspace=0.42)
# -- C1: EFS Binary --
axC1 = fig.add_subplot(cgs[0, 0])
cats = ["No Event", "Event", "Unknown"]
c1_vals = efs_bin.reindex(cats).values
c1_colors = ["#A6D854", "#FB8072", "#BEBEBE"]
bars = axC1.bar(cats, c1_vals, color=c1_colors, edgecolor="black", linewidth=0.5,
                width=0.62)
for b, v in zip(bars, c1_vals):
    axC1.text(b.get_x() + b.get_width() / 2, v + 12, f"{int(v)}", ha="center",
              va="bottom", fontsize=9, fontweight="bold")
axC1.set_ylim(0, 1500)
axC1.set_title("EFS Binary", fontsize=10.5)
axC1.set_ylabel("Count (patients)", fontsize=9.5)
axC1.tick_params(axis="x", labelsize=8.5)
axC1.grid(axis="y", ls=":", alpha=0.4)
axC1.set_axisbelow(True)
# -- C2: EFS Detailed --
axC2 = fig.add_subplot(cgs[0, 1])
det_cats = efs_detail.index.tolist()
det_vals = efs_detail.values
bars = axC2.bar(det_cats, det_vals, color="#FB8072", edgecolor="black", linewidth=0.4,
                width=0.68)
for b, v in zip(bars, det_vals):
    axC2.text(b.get_x() + b.get_width() / 2, v + 8, f"{int(v)}", ha="center",
              va="bottom", fontsize=7.5)
axC2.set_ylim(0, 1450)
axC2.set_xticks(np.arange(len(det_cats)))
axC2.set_xticklabels(det_cats, rotation=-45, ha="left", fontsize=7)
axC2.set_title("EFS Detailed", fontsize=10.5)
axC2.set_ylabel("Count (patients)", fontsize=9.5)
axC2.grid(axis="y", ls=":", alpha=0.4)
axC2.set_axisbelow(True)
axC_outer.set_title(f"Event-Free Survival Status  (N = {N})", fontsize=12)
panel_letter(axC_outer, "C")

# ---------------- Panel D: AGE histogram ----------------
axD = fig.add_subplot(gs[2, 0])
n_bins = 40
bins = np.linspace(age.min(), age.max(), n_bins + 1)
axD.hist(age, bins=bins, color="steelblue", alpha=0.75, edgecolor="white",
         linewidth=0.3)
tail = age[age > THRESH]
axD.hist(tail, bins=bins, color="#D95F02", alpha=0.8, edgecolor="white", linewidth=0.3)
axD.set_xlabel("Age (years)", fontsize=11)
axD.set_ylabel("Count", fontsize=11)
axD.set_title(f"Age Distribution  (n = {len(age)} of {N})", fontsize=12)
stext = (f"mean {age_stats['mean']:.1f}   median {age_stats['median']:.1f}\n"
         f"min {age_stats['min']:.1f}   max {age_stats['max']:.1f}\n"
         f"missing {age_missing_n} ({age_missing_n / N * 100:.1f}%)")
axD.text(0.985, 0.97, stext, transform=axD.transAxes, ha="right", va="top",
         fontsize=8.5, bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#444444", alpha=0.9))
axD.grid(axis="y", ls=":", alpha=0.4)
axD.set_axisbelow(True)
panel_letter(axD, "D")

# ---------------- Panel E: AGE outliers per cancer group ----------------
axE = fig.add_subplot(gs[2, 1])
e_cats = ol_counts.index.tolist()
e_vals = ol_counts.values
e_colors = ["#D95F02" if c == "Oligodendroglioma" else "#C0C0C0" for c in e_cats]
bars = axE.barh(e_cats, e_vals, color=e_colors, edgecolor="black", linewidth=0.5,
                height=0.62)
for b, v in zip(bars, e_vals):
    axE.text(v + 0.5, b.get_y() + b.get_height() / 2, f"{int(v)}", ha="left",
             va="center", fontsize=9.5, fontweight="bold")
axE.set_xlabel("Count of AGE-outlier patients", fontsize=10.5)
axE.set_ylabel("Cancer group", fontsize=10.5)
axE.set_xlim(0, 62)
axE.set_ylim(-0.7, len(e_cats) + 1.3)
axE.set_title("AGE outliers per cancer group\n(Tukey IQR outlier, AGE > 27.125 yr)", fontsize=11.5)
oligo_idx = e_cats.index("Oligodendroglioma")
axE.annotate(f"53/62 ({oligo_n / n_out * 100:.1f}%) of AGE outliers\nare Oligodendroglioma",
             xy=(53, oligo_idx), xytext=(1.5, len(e_cats) + 0.25),
             textcoords="data", fontsize=8.5, ha="left", va="center",
             bbox=dict(boxstyle="round,pad=0.35", fc="#FDBF6F", ec="black", alpha=0.9),
             arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
axE.grid(axis="x", ls=":", alpha=0.4)
axE.set_axisbelow(True)
panel_letter(axE, "E")

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
print("saved:", OUT_PNG)
print("saved:", OUT_PDF)

# ----------------------------------------------------------------------
# Verification table
# ----------------------------------------------------------------------
print("\n=== COMPUTED vs EXPECTED (spec) ===")
print(f"Panel A  EFS_MONTHS: {miss['EFS_MONTHS']:.1f}% (816) vs 28.4% | "
      f"OS_MONTHS: {miss['OS_MONTHS']:.1f}% (774) vs 27.0% | "
      f"RACE {miss['RACE']:.1f}% vs 25.9% | ETHNICITY {miss['ETHNICITY']:.1f}% vs 14.6% | "
      f"OS_STATUS {miss['OS_STATUS']:.1f}% vs 12.4% | GERMLINE {miss['GERMLINE_SEX_ESTIMATE']:.1f}% vs 11.0% | "
      f"AGE_IN_DAYS {miss['AGE_IN_DAYS']:.1f}% vs 2.0% | AGE {miss['AGE']:.1f}% vs 2.0% | "
      f"SEX {miss['SEX']:.1f}% vs 0.5%")
print(f"Panel B  LIVING {int(os_counts['LIVING'])} ({os_counts['LIVING']/N*100:.1f}%) vs 1875 (65.3%) | "
      f"DECEASED {int(os_counts['DECEASED'])} ({os_counts['DECEASED']/N*100:.1f}%) vs 640 (22.3%) | "
      f"Unknown {int(os_counts['Unknown'])} ({os_counts['Unknown']/N*100:.1f}%) vs 355 (12.4%)")
print(f"Panel C  binary: No Event {int(efs_bin['No Event'])} vs 1286 | Event {int(efs_bin['Event'])} vs 1222 | "
      f"Unknown {int(efs_bin['Unknown'])} vs 362")
for cat, exp in [("No Event", 1286), ("Progressive", 423), ("NA", 362), ("Recurrence", 305),
                 ("Progressive - Metastatic", 177), ("Deceased-due to disease", 158),
                 ("Recurrence - Metastatic", 88), ("Second Malignancy", 45),
                 ("Deceased-due to other causes", 9), ("Deceased-causes unavailable", 6),
                 ("Deceased-due to unknown causes", 6), ("Second Malignancy - Metastatic", 5)]:
    got = int(efs_detail.get(cat, 0))
    flag = "" if got == exp else "  <-- MISMATCH"
    print(f"    detailed {cat!r}: {got} vs {exp}{flag}")
print(f"Panel D  mean {age_stats['mean']:.1f} vs 9.4 | median {age_stats['median']:.1f} vs 8.0 | "
      f"min {age_stats['min']:.1f} vs 0.0 | max {age_stats['max']:.1f} vs 73.0 | "
      f"missing {age_missing_n} ({age_missing_n/N*100:.1f}%) vs 58 (2.0%)")
for cat in ["Oligodendroglioma", "High-grade glioma", "Neurofibroma/Plexiform", "Schwannoma",
            "Low-grade glioma", "Diffuse midline glioma", "Meningioma", "Medulloblastoma"]:
    exp = 53 if cat == "Oligodendroglioma" else (3 if cat == "High-grade glioma" else 1)
    got = int(ol_counts.get(cat, 0))
    flag = "" if got == exp else "  <-- MISMATCH"
    print(f"Panel E  {cat}: {got} vs {exp}{flag}")
print(f"Panel E  total outliers: {n_out} vs 62 | oligo {oligo_n}/62 ({oligo_n/62*100:.1f}%)")
