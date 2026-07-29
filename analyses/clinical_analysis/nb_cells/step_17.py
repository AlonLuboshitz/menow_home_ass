patients = read_patients()
patients = clean_pred(patients); patients = clean_os(patients); patients = clean_efs(patients)
patients['hp'] = ~patients['CANCER_PREDISPOSITIONS'].isin(['No predisposition','Unknown'])
print(f'With predisposition: {patients["hp"].sum()} / {len(patients)}')
fig = make_subplots(rows=2,cols=2,subplot_titles=('OS by Pred','EFS by Pred','Age by Pred','OS Event'),
    specs=[[{'type':'scatter'},{'type':'scatter'}],[{'type':'box'},{'type':'bar'}]])
for i,(lab,has) in enumerate([('No',False),('Yes',True)]):
    sub = patients[patients['hp']==has]
    os_s = sub[['OS_MONTHS','os_event']].dropna()
    if len(os_s)>3:
        fig = add_km(fig,kaplan_meier(os_s['OS_MONTHS'],os_s['os_event']),lab,px.colors.qualitative.Set1[i])
    ef_s = sub[['EFS_MONTHS','efs_event']].dropna()
    if len(ef_s)>3:
        km = kaplan_meier(ef_s['EFS_MONTHS'],ef_s['efs_event'])
        fig.add_trace(go.Scatter(x=km['t'],y=km['s'],mode='lines',name=lab,
            line=dict(color=px.colors.qualitative.Set1[i],width=2,shape='hv'),legendgroup=lab,showlegend=False),row=1,col=2)
    fig.add_trace(go.Box(y=sub['AGE'].dropna(),name=lab,boxmean='sd',marker_color=px.colors.qualitative.Set1[i]),row=2,col=1)
fig.update_layout(height=600,title='Predisposition vs Outcome')
fig.show()
