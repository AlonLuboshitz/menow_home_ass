patients = read_patients()
patients = clean_os(patients); patients = clean_efs(patients)
samples = read_samples()
merged = samples.merge(patients,on='PATIENT_ID',how='left',suffixes=('','_p'))
top6 = merged['CANCER_GROUP'].value_counts().head(6).index.tolist()
fig = make_subplots(rows=2,cols=1,subplot_titles=('OS by Cancer Group','EFS by Cancer Group'),vertical_spacing=0.15)
colors = px.colors.qualitative.Set1
os_d=[]; ef_d=[]
for i,cg in enumerate(top6):
    sub = merged[merged['CANCER_GROUP']==cg]
    os_s = sub[['OS_MONTHS','os_event']].dropna()
    if len(os_s)>5:
        fig = add_km(fig,kaplan_meier(os_s['OS_MONTHS'],os_s['os_event']),cg,colors[i%len(colors)])
        os_d.append((os_s['OS_MONTHS'],os_s['os_event']))
    ef_s = sub[['EFS_MONTHS','efs_event']].dropna()
    if len(ef_s)>5:
        km = kaplan_meier(ef_s['EFS_MONTHS'],ef_s['efs_event'])
        fig.add_trace(go.Scatter(x=km['t'],y=km['s'],mode='lines',name=cg,
            line=dict(color=colors[i%len(colors)],width=2,shape='hv'),legendgroup=cg,showlegend=False),row=2,col=1)
        ef_d.append((ef_s['EFS_MONTHS'],ef_s['efs_event']))
if len(os_d)>=2:
    po=logrank_multi(os_d); pe=logrank_multi(ef_d)
    fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'OS p={po:.4f}',showarrow=False,font=dict(size=11,color='darkblue'),row=1,col=1)
    fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'EFS p={pe:.4f}',showarrow=False,font=dict(size=11,color='darkred'),row=2,col=1)
fig.update_layout(height=650,title='Survival by Cancer Group')
fig.update_yaxes(range=[-0.05,1.05])
fig.show()
