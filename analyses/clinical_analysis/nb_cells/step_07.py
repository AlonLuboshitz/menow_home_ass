df = read_samples()
df = clean_tf_tp(df)
fig = make_subplots(rows=1,cols=3,subplot_titles=('Tumor Fraction','Tumor Ploidy','Fraction vs Ploidy'),
    specs=[[{'type':'histogram'},{'type':'histogram'},{'type':'scatter'}]])
tf = df[df['TF_group']=='Measured']['TUMOR_FRACTION'].dropna()
fig.add_trace(go.Histogram(x=tf,nbinsx=40,marker_color='steelblue',opacity=0.75),row=1,col=1)
nuk = (df['TF_group']=='Unknown').sum()
fig.add_annotation(text=f'Unknown: {nuk} ({nuk/len(df):.1%})',xref='paper',yref='paper',x=0.5,y=-0.3,showarrow=False,font=dict(color='gray'),row=1,col=1)
tp = df[df['TP_group']=='Measured']['TUMOR_PLOIDY'].dropna()
fig.add_trace(go.Histogram(x=tp,nbinsx=30,marker_color='darkgreen',opacity=0.75),row=1,col=2)
nuk2 = (df['TP_group']=='Unknown').sum()
fig.add_annotation(text=f'Unknown: {nuk2} ({nuk2/len(df):.1%})',xref='paper',yref='paper',x=0.5,y=-0.3,showarrow=False,font=dict(color='gray'),row=1,col=2)
sc = df.dropna(subset=['TUMOR_FRACTION','TUMOR_PLOIDY','CANCER_GROUP'])
top = sc['CANCER_GROUP'].value_counts().head(8).index
for cg in top:
    s = sc[sc['CANCER_GROUP']==cg]
    fig.add_trace(go.Scatter(x=s['TUMOR_FRACTION'],y=s['TUMOR_PLOIDY'],mode='markers',name=cg,marker=dict(size=5,opacity=0.6)),row=1,col=3)
fig.update_layout(height=450,title='Tumor Purity & Ploidy',xaxis3_title='Fraction',yaxis3_title='Ploidy')
fig.show()
