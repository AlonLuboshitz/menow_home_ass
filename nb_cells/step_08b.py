t15 = df['MOLECULAR_SUBTYPE'].value_counts().head(15).index
t10 = df['CANCER_GROUP'].value_counts().head(10).index
ct = pd.crosstab(df['MOLECULAR_SUBTYPE'],df['CANCER_GROUP'])
cf = ct.loc[ct.index.intersection(t15),ct.columns.intersection(t10)]
cn = cf.div(cf.sum(axis=1),axis=0).fillna(0)
fig = go.Figure(data=go.Heatmap(z=cn.values,x=cn.columns,y=cn.index,text=cf.values,texttemplate='%{text}',
    textfont=dict(size=9),colorscale='Blues',colorbar=dict(title='Proportion')))
fig.update_layout(title='Subtype x Cancer Group (Row %)',xaxis_tickangle=-45,height=550,yaxis=dict(autorange='reversed'))
fig.show()
