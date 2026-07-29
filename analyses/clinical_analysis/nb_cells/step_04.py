import re
df = read_patients()
df = clean_pred(df)
df = clean_race_eth(df)
def explode(df):
    rows = []
    for _,r in df.iterrows():
        p = r['CANCER_PREDISPOSITIONS']
        if p in ('Unknown','No predisposition'): rows.append({**r,'pred':p})
        elif '),' in str(p):
            parts = re.split(r'\),\s*',str(p))
            for i,pr in enumerate(parts):
                if i < len(parts)-1: pr = pr + ')'
                rows.append({**r,'pred':pr.strip()})
        else: rows.append({**r,'pred':p})
    return pd.DataFrame(rows)
dex = explode(df)
total_p = df['PATIENT_ID'].nunique()
print(f'Patients: {total_p}, Exploded rows: {len(dex)}')
multi = df[df['CANCER_PREDISPOSITIONS'].str.contains(r'\),',na=False,regex=True)]
print(f'Multi-syndrome patients: {len(multi)}')
