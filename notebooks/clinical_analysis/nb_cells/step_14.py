patients = read_patients()
samples = read_samples()
merged = samples.merge(patients,on='PATIENT_ID',how='left',suffixes=('','_p'))
top8 = merged['CANCER_GROUP'].value_counts().head(8).index.tolist()
pd8 = merged[merged['CANCER_GROUP'].isin(top8)].dropna(subset=['AGE'])
fig = go.Figure()
for i,cg in enumerate(top8):
    sub = pd8[pd8['CANCER_GROUP']==cg]['AGE']
    fig.add_trace(go.Box(y=sub,name=cg,boxmean='sd',marker_color=px.colors.qualitative.Plotly[i]))
groups = [merged[merged['CANCER_GROUP']==cg]['AGE'].dropna() for cg in top8]
stat,p_kw = kruskal(*groups)
fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.05,text=f'KW: H={stat:.2f}, p={p_kw:.4f}',showarrow=False,font=dict(size=12,color='darkred'))
fig.update_layout(title='Age at Diagnosis by Cancer Group',yaxis_title='Age (years)',height=500,xaxis_tickangle=-45)
fig.show()
if p_kw<0.05:
    print('Post-hoc (MW p<0.01):')
    for i,c1 in enumerate(top8):
        for j,c2 in enumerate(top8):
            if i>=j: continue
            g1 = merged[merged['CANCER_GROUP']==c1]['AGE'].dropna()
            g2 = merged[merged['CANCER_GROUP']==c2]['AGE'].dropna()
            if len(g1)>5 and len(g2)>5:
                _,p = mannwhitneyu(g1,g2,alternative='two-sided')
                if p<0.01: print(f'  {c1:35s} vs {c2:35s}: p={p:.6f}')
