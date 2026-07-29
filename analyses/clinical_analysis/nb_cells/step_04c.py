top10 = pc.head(10)
dp = top10['Pred'].iloc[0]
wp = dex[dex['pred']==dp]
wp2 = df[df['PATIENT_ID'].isin(wp['PATIENT_ID'].unique())]
fig = make_subplots(rows=2,cols=2,subplot_titles=('Age: With vs Without','Sex','Prevalence %','Summary'),
    specs=[[{'type':'box'},{'type':'bar'}],[{'type':'bar'},{'type':'table'}]])
fig.add_trace(go.Box(y=wp['AGE'].dropna(),name='With',marker_color=px.colors.qualitative.Plotly[0]),row=1,col=1)
fig.add_trace(go.Box(y=df[~df['PATIENT_ID'].isin(wp['PATIENT_ID'])]['AGE'].dropna(),name='Without',marker_color='lightgray'),row=1,col=1)
sx = wp2['SEX'].fillna('NaN').value_counts()
fig.add_trace(go.Bar(x=sx.index,y=sx.values,marker_color='lightcoral',text=sx.values,textposition='outside',showlegend=False),row=1,col=2)
fig.add_trace(go.Bar(x=top10['Pred'],y=top10['%'],marker_color='steelblue',text=top10['%'],textposition='outside',showlegend=False),row=2,col=1)
td = top10.rename(columns={'%':'Pct'})
fig.add_trace(go.Table(header=dict(values=list(td.columns),fill_color='lightblue',align='left'),
    cells=dict(values=[td[c] for c in td.columns],align='left',height=22)),row=2,col=2)
fig.update_layout(height=650,title='Predisposition Explorer')
fig.show()
