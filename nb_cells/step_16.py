patients = read_patients()
samples = read_samples(); samples = clean_tf_tp(samples)
merged = samples.merge(patients,on='PATIENT_ID',how='left',suffixes=('','_p'))
fig = make_subplots(rows=1,cols=2,subplot_titles=('TF by Cancer Group','TF by Tumor Type'))
top8 = merged['CANCER_GROUP'].value_counts().head(8).index.tolist()
for i,cg in enumerate(top8):
    sub = merged[merged['CANCER_GROUP']==cg]['TUMOR_FRACTION'].dropna()
    if len(sub)>0: fig.add_trace(go.Box(y=sub,name=cg,boxmean='sd',marker_color=px.colors.qualitative.Plotly[i]),row=1,col=1)
top_tt = ['primary','metastatic','progression','recurrence']
for tt in top_tt:
    sub = merged[merged['TUMOR_TYPE']==tt]['TUMOR_FRACTION'].dropna()
    if len(sub)>0: fig.add_trace(go.Box(y=sub,name=tt,boxmean='sd',marker_color=px.colors.qualitative.Set2[top_tt.index(tt)]),row=1,col=2)
fig.update_layout(height=500,title='Tumor Fraction Analysis',yaxis_title='Fraction',xaxis_tickangle=-45,xaxis2_tickangle=-45)
fig.show()
for name,grps in [('Cancer Group',top8),('Tumor Type',top_tt)]:
    g = [merged[merged['CANCER_GROUP']==cg]['TUMOR_FRACTION'].dropna() for cg in grps if len(merged[merged['CANCER_GROUP']==cg]['TUMOR_FRACTION'].dropna())>5]
    if len(g)>=2:
        s,p = kruskal(*g); print(f'KW ({name}): H={s:.2f}, p={p:.4f}')
