df = read_samples()
pg = df.groupby('PATIENT_ID')['CANCER_GROUP'].apply(set).reset_index()
pg['n'] = pg['CANCER_GROUP'].apply(len)
fig = make_subplots(rows=1,cols=2,subplot_titles=('Groups per Patient','Multi-Group Patients'),
    specs=[[{'type':'bar'},{'type':'table'}]])
gc = pg['n'].value_counts().sort_index()
fig.add_trace(go.Bar(x=gc.index.astype(str),y=gc.values,marker_color='steelblue',text=gc.values,textposition='outside',showlegend=False),row=1,col=1)
mp = pg[pg['n']>1].copy()
mp['Groups'] = mp['CANCER_GROUP'].apply(lambda x:', '.join(sorted(str(v) for v in x)))
mp = mp.sort_values('n',ascending=False).head(20)
if len(mp)>0:
    t = mp[['PATIENT_ID','n','Groups']].head(15)
    fig.add_trace(go.Table(header=dict(values=['Patient','N','Groups'],fill_color='lightblue',align='left',font=dict(size=10)),
        cells=dict(values=[t['PATIENT_ID'],t['n'],t['Groups']],align='left',height=22,font=dict(size=9))),row=1,col=2)
fig.update_layout(height=500,title='Multi-Cancer-Group Analysis')
fig.show()
