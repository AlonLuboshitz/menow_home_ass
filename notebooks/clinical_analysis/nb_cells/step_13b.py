# Per-cancer-group subtype analysis
div = merged.groupby('CANCER_GROUP')['MOLECULAR_SUBTYPE'].nunique().sort_values(ascending=False)
print('Subtype diversity:\n',div)
rich = div[div>=5].index.tolist()
print(f'\nGroups with >=5 subtypes: {rich}')
for cg in rich[:3]:
    cgm = merged[merged['CANCER_GROUP']==cg]
    top = cgm['MOLECULAR_SUBTYPE'].value_counts().head(5).index.tolist()
    if len(top)<2: continue
    fig = make_subplots(rows=2,cols=1,subplot_titles=(f'OS -- {cg}',f'EFS -- {cg}'),vertical_spacing=0.15)
    os_d2=[]; ef_d2=[]
    for j,st in enumerate(top):
        sub = cgm[cgm['MOLECULAR_SUBTYPE']==st]
        os_s = sub[['OS_MONTHS','os_event']].dropna()
        if len(os_s)>3:
            fig = add_km(fig,kaplan_meier(os_s['OS_MONTHS'],os_s['os_event']),st,colors[j%len(colors)])
            os_d2.append((os_s['OS_MONTHS'],os_s['os_event']))
        ef_s = sub[['EFS_MONTHS','efs_event']].dropna()
        if len(ef_s)>3:
            km = kaplan_meier(ef_s['EFS_MONTHS'],ef_s['efs_event'])
            fig.add_trace(go.Scatter(x=km['t'],y=km['s'],mode='lines',name=st,
                line=dict(color=colors[j%len(colors)],width=2,shape='hv'),legendgroup=st,showlegend=False),row=2,col=1)
            ef_d2.append((ef_s['EFS_MONTHS'],ef_s['efs_event']))
    if len(os_d2)>=2:
        fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'OS p={logrank_multi(os_d2):.4f}',showarrow=False,font=dict(size=10,color='darkblue'),row=1,col=1)
        fig.add_annotation(xref='paper',yref='paper',x=0.5,y=1.02,text=f'EFS p={logrank_multi(ef_d2):.4f}',showarrow=False,font=dict(size=10,color='darkred'),row=2,col=1)
    fig.update_layout(height=550,title=f'Survival by Subtype -- {cg}')
    fig.update_yaxes(range=[-0.05,1.05])
    fig.show()
