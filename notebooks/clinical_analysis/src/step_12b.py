# Step 12: Per-group summary table for ALL cancer groups (sorted by N samples).
def fmt_p(p): return f'{p:.4g}' if p is not None else 'n/a'
rows = []
for cg in merged['CANCER_GROUP'].value_counts().index.tolist():
    s = merged[merged['CANCER_GROUP']==cg]
    os_s = s[['OS_MONTHS','os_event']].dropna()
    ef_s = s[['EFS_MONTHS','efs_event']].dropna()
    n_os = len(os_s); n_efs = len(ef_s)
    po = os_p.get(cg) if n_os>=MIN_N else None
    qo = os_q.get(cg) if n_os>=MIN_N else None
    pe = ef_p.get(cg) if n_efs>=MIN_N else None
    qe = ef_q.get(cg) if n_efs>=MIN_N else None
    rows.append({'Cancer Group':cg,'N samples':len(s),'N patients':s['PATIENT_ID'].nunique(),
        'N_OS':n_os,'Med OS':f'{os_s["OS_MONTHS"].median():.1f}' if n_os>0 else 'N/A',
        'OS p(vs others)':fmt_p(po),'OS q(FDR)':fmt_p(qo),'OS sig':star(po,qo) if (po is not None and qo is not None) else 'n/a',
        'N_EFS':n_efs,'Med EFS':f'{ef_s["EFS_MONTHS"].median():.1f}' if n_efs>0 else 'N/A',
        'EFS p(vs others)':fmt_p(pe),'EFS q(FDR)':fmt_p(qe),'EFS sig':star(pe,qe) if (pe is not None and qe is not None) else 'n/a'})
tbl = pd.DataFrame(rows)
print(tbl.to_string(index=False))
print('Significance: * p<0.05, ** FDR<0.05 (BH) -- log-rank of this group vs all other groups pooled.')
