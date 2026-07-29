m2 = pg[pg['n']==2].copy()
if len(m2)>0:
    m2['gl'] = m2['CANCER_GROUP'].apply(lambda x:sorted(str(v) for v in x))
    m2['g1'] = m2['gl'].apply(lambda x:x[0])
    m2['g2'] = m2['gl'].apply(lambda x:x[1])
    co = pd.crosstab(m2['g1'],m2['g2'])
    fig = go.Figure(data=go.Heatmap(z=co.values,x=co.columns,y=co.index,text=co.values,
        texttemplate='%{text}',colorscale='Blues'))
    fig.update_layout(title='Cancer Group Co-occurrence (2 groups)',xaxis_tickangle=-45,height=450)
    fig.show()
else: print('No patients with exactly 2 groups.')
