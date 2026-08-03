# Step 12: Survival by Cancer Group -- OS/EFS KM curves, per-group log-rank (each group vs all others).
patients = read_patients()
patients = clean_os(patients); patients = clean_efs(patients)
samples = read_samples()
merged = samples.merge(patients,on='PATIENT_ID',how='left',suffixes=('','_p'))
merged = merged[merged['CANCER_GROUP'].notna() & (merged['CANCER_GROUP'].astype(str).str.strip()!='')]

# BH-FDR (Benjamini-Hochberg) + significance labels
def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0: return np.array([])
    order = np.argsort(pvals)
    adj = pvals[order]*n/np.arange(1,n+1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(adj,0,1)
    return out

def star(p,q):
    """Table label: ** FDR<0.05, * p<0.05, ns otherwise."""
    if q<0.05: return '**'
    if p<0.05: return '*'
    return 'ns'

def legend_sig(p,q):
    """Legend suffix: only marks significant differences."""
    if q<0.05: return '**'
    if p<0.05: return '*'
    return ''

# Complete (time,event) pairs per group
osd = merged[['CANCER_GROUP','OS_MONTHS','os_event']].dropna()
efd = merged[['CANCER_GROUP','EFS_MONTHS','efs_event']].dropna()
os_n  = osd.groupby('CANCER_GROUP').size().sort_values(ascending=False)
efs_n = efd.groupby('CANCER_GROUP').size().sort_values(ascending=False)

MIN_N = 20      # minimum complete records to enter a per-group test
MAX_PLOT = 15   # maximum groups drawn
test_groups_os  = os_n[os_n>=MIN_N].index.tolist()
test_groups_efs = efs_n[efs_n>=MIN_N].index.tolist()
plot_groups = os_n[os_n>=MIN_N].head(MAX_PLOT).index.tolist()

# Per-group log-rank: group vs ALL other groups pooled (OS and EFS)
os_p = {g: logrank2(osd.loc[osd['CANCER_GROUP']==g,'OS_MONTHS'], osd.loc[osd['CANCER_GROUP']==g,'os_event'],
                    osd.loc[osd['CANCER_GROUP']!=g,'OS_MONTHS'], osd.loc[osd['CANCER_GROUP']!=g,'os_event'])
        for g in test_groups_os}
ef_p = {g: logrank2(efd.loc[efd['CANCER_GROUP']==g,'EFS_MONTHS'], efd.loc[efd['CANCER_GROUP']==g,'efs_event'],
                    efd.loc[efd['CANCER_GROUP']!=g,'EFS_MONTHS'], efd.loc[efd['CANCER_GROUP']!=g,'efs_event'])
        for g in test_groups_efs}
os_q = dict(zip(test_groups_os, bh_fdr([os_p[g] for g in test_groups_os])))
ef_q = dict(zip(test_groups_efs, bh_fdr([ef_p[g] for g in test_groups_efs])))

print(f'Cancer groups with complete OS: {len(os_n)} | EFS: {len(efs_n)}')
print(f'Tested vs others (N>={MIN_N}): {len(test_groups_os)} OS, {len(test_groups_efs)} EFS | plotted: {len(plot_groups)}')
print(f'FDR<0.05 vs all others: {sum(1 for g in test_groups_os if os_q[g]<0.05)}/{len(test_groups_os)} OS, '
      f'{sum(1 for g in test_groups_efs if ef_q[g]<0.05)}/{len(test_groups_efs)} EFS')

# KM plot: OS (top) + EFS (bottom); legend driven by OS, n = complete OS records
fig = make_subplots(rows=2,cols=1,subplot_titles=('OS by Cancer Group','EFS by Cancer Group'),vertical_spacing=0.15)
colors = (px.colors.qualitative.Set1 + px.colors.qualitative.Dark2 + px.colors.qualitative.Set3)
os_d=[]; ef_d=[]
for i,cg in enumerate(plot_groups):
    sub = merged[merged['CANCER_GROUP']==cg]
    os_s = sub[['OS_MONTHS','os_event']].dropna()
    label = f'{cg} (n={len(os_s)}){legend_sig(os_p.get(cg,1.0),os_q.get(cg,1.0))}'
    fig = add_km(fig,kaplan_meier(os_s['OS_MONTHS'],os_s['os_event']),label,colors[i%len(colors)])
    os_d.append((os_s['OS_MONTHS'],os_s['os_event']))
    ef_s = sub[['EFS_MONTHS','efs_event']].dropna()
    km = kaplan_meier(ef_s['EFS_MONTHS'],ef_s['efs_event'])
    fig.add_trace(go.Scatter(x=km['t'],y=km['s'],mode='lines',name=label,
        line=dict(color=colors[i%len(colors)],width=2,shape='hv'),legendgroup=label,showlegend=False),row=2,col=1)
    ef_d.append((ef_s['EFS_MONTHS'],ef_s['efs_event']))
if len(os_d)>=2:
    po=logrank_multi(os_d); pe=logrank_multi(ef_d)
    fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'global OS log-rank p={po:.4f}',showarrow=False,font=dict(size=11,color='darkblue'),row=1,col=1)
    fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'global EFS log-rank p={pe:.4f}',showarrow=False,font=dict(size=11,color='darkred'),row=2,col=1)
fig.update_layout(height=800,title='Survival by Cancer Group (n = complete OS records; * p<0.05, ** FDR<0.05 vs all other groups)',
                  legend=dict(font=dict(size=10)))
fig.update_yaxes(range=[-0.05,1.05])
fig.show()
