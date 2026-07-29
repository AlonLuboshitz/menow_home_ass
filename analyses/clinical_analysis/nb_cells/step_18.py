patients = read_patients(); samples = read_samples()
merged = samples.merge(patients,on='PATIENT_ID',how='left',suffixes=('','_p'))
ct = pd.crosstab(merged['CNS_REGION'],merged['CANCER_GROUP'])
tr = merged['CNS_REGION'].value_counts().head(8).index
tc = merged['CANCER_GROUP'].value_counts().head(10).index
cf = ct.loc[ct.index.intersection(tr),ct.columns.intersection(tc)]
cn = cf.div(cf.sum(axis=1),axis=0).fillna(0)
fig = go.Figure(data=go.Heatmap(z=cn.values,x=cn.columns,y=cn.index,text=cf.values,texttemplate='%{text}',
    textfont=dict(size=10),colorscale='Purples',colorbar=dict(title='Proportion')))
fig.update_layout(title='CNS Region x Cancer Group (Row-Normalized)',xaxis_tickangle=-45,height=550)
fig.show()
