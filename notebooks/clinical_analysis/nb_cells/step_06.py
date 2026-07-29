df = read_samples()
fig = make_subplots(rows=2,cols=2,subplot_titles=('Broad Histology','Cancer Group','CNS Region','Tumor Type'),
    specs=[[{'type':'bar'},{'type':'bar'}],[{'type':'bar'},{'type':'bar'}]])
hc = df['BROAD_HISTOLOGY'].value_counts().head(12)
fig.add_trace(go.Bar(x=hc.index,y=hc.values,marker_color='lightblue',text=hc.values,textposition='outside',showlegend=False),row=1,col=1)
cg = df['CANCER_GROUP'].value_counts().head(12)
fig.add_trace(go.Bar(x=cg.index,y=cg.values,marker_color='lightgreen',text=cg.values,textposition='outside',showlegend=False),row=1,col=2)
cr = df['CNS_REGION'].value_counts()
fig.add_trace(go.Bar(x=cr.index,y=cr.values,marker_color='lightsalmon',text=cr.values,textposition='outside',showlegend=False),row=2,col=1)
tt = df['TUMOR_TYPE'].value_counts()
comm = tt[tt>=30]; rare = tt[tt<30]
if len(rare)>0:
    oth = pd.Series({'Other':rare.sum()})
    tt_plot = pd.concat([comm,oth])
else: tt_plot = comm
fig.add_trace(go.Bar(x=tt_plot.index,y=tt_plot.values,marker_color='plum',text=tt_plot.values,textposition='outside',showlegend=False),row=2,col=2)
fig.update_layout(title='Sample Cancer Type Distributions',height=650,xaxis_tickangle=-45,xaxis2_tickangle=-45,xaxis3_tickangle=-45,xaxis4_tickangle=-45)
fig.show()
