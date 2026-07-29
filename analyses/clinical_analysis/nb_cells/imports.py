import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import kruskal, mannwhitneyu, chi2_contingency
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = "/home/alon/menow_home_ass/PBTA_RNA"
PATIENT_FILE = f"{DATA_DIR}/data_clinical_patient_attributes.txt"
SAMPLE_FILE = f"{DATA_DIR}/data_clinical_sample_attributes.txt"

def read_patients():
    return pd.read_csv(PATIENT_FILE, sep="\t", header=4,
                       dtype={"AGE": float, "AGE_IN_DAYS": float,
                              "OS_MONTHS": float, "EFS_MONTHS": float})

def read_samples():
    return pd.read_csv(SAMPLE_FILE, sep="\t", header=4)

# OS/EFS cleaners
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

# Categorical cleaners
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

# Kaplan-Meier helper
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

# Log-rank tests
def logrank2(t1,e1,t2,e2):
    from scipy.stats import chi2
    all_t = sorted(set(pd.concat([pd.Series(t1.dropna()),pd.Series(t2.dropna())]).dropna()))
    if len(all_t) < 2: return 1.0
    o1e=0; v=0
    d1=pd.DataFrame({"t":t1,"e":e1}).dropna()
    d2=pd.DataFrame({"t":t2,"e":e2}).dropna()
    for t in all_t:
        r1=(d1["t"]>=t).sum(); r2=(d2["t"]>=t).sum(); nr=r1+r2
        if nr==0: continue
        o1=int(((d1["t"]==t)&(d1["e"]==1)).sum())
        o2=int(((d2["t"]==t)&(d2["e"]==1)).sum())
        ot=o1+o2
        if ot==0: continue
        e1=ot*r1/nr; o1e+=(o1-e1)
        if nr>1: v+=ot*(r1/nr)*(r2/nr)*(nr-ot)/(nr-1)
    if v<=0: return 1.0
    return 1-chi2.cdf(o1e**2/v,1)

def logrank_multi(groups):
    from scipy.stats import chi2
    import numpy as np
    ng=len(groups)
    if ng<2: return 1.0
    all_t=sorted(set(pd.concat([pd.Series(g[0].dropna()) for g in groups]).dropna()))
    if len(all_t)<2: return 1.0
    O=np.zeros(ng);E=np.zeros(ng);V=np.zeros((ng,ng))
    for t in all_t:
        ar=np.array([(g[0]>=t).sum() for g in groups]); nr=ar.sum()
        if nr==0: continue
        ev=np.array([((g[0]==t)&(g[1]==1)).sum() for g in groups]); ot=ev.sum()
        if ot==0: continue
        O+=ev;E+=ot*ar/nr
        if nr>1:
            for i in range(ng):
                for j in range(ng):
                    if i==j: V[i,j]+=ot*ar[i]/nr*(1-ar[i]/nr)*(nr-ot)/(nr-1)
                    else: V[i,j]-=ot*ar[i]/nr*ar[j]/nr*(nr-ot)/(nr-1)
    try: return 1-chi2.cdf((O-E)@np.linalg.pinv(V)@(O-E),ng-1)
    except: return 1.0

print("Imports and helpers loaded.")
