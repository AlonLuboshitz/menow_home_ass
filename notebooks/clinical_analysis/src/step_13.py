patients = read_patients()
patients = clean_os(patients); patients = clean_efs(patients)
samples = read_samples(); samples = clean_subtype(samples)
merged = samples.merge(patients,on='PATIENT_ID',how='left',suffixes=('','_p'))
top6 = merged['MOLECULAR_SUBTYPE'].value_counts().head(6).index.tolist()
fig = make_subplots(rows=2,cols=1,subplot_titles=('OS by Subtype','EFS by Subtype'),vertical_spacing=0.15)
colors = px.colors.qualitative.Set1 + px.colors.qualitative.Set2
os_d=[]; ef_d=[]
for i,st in enumerate(top6):
    sub = merged[merged['MOLECULAR_SUBTYPE']==st]
    os_s = sub[['OS_MONTHS','os_event']].dropna()
    if len(os_s)>5:
        fig = add_km(fig,kaplan_meier(os_s['OS_MONTHS'],os_s['os_event']),st,colors[i%len(colors)])
        os_d.append((os_s['OS_MONTHS'],os_s['os_event']))
    ef_s = sub[['EFS_MONTHS','efs_event']].dropna()
    if len(ef_s)>5:
        km = kaplan_meier(ef_s['EFS_MONTHS'],ef_s['efs_event'])
        fig.add_trace(go.Scatter(x=km['t'],y=km['s'],mode='lines',name=st,
            line=dict(color=colors[i%len(colors)],width=2,shape='hv'),legendgroup=st,showlegend=False),row=2,col=1)
        ef_d.append((ef_s['EFS_MONTHS'],ef_s['efs_event']))
if len(os_d)>=2:
    fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'OS p={logrank_multi(os_d):.4f}',showarrow=False,font=dict(size=11,color='darkblue'),row=1,col=1)
    fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'EFS p={logrank_multi(ef_d):.4f}',showarrow=False,font=dict(size=11,color='darkred'),row=2,col=1)
fig.update_layout(height=650,title='Survival by Subtype (Global)')
fig.update_yaxes(range=[-0.05,1.05])
fig.show()
