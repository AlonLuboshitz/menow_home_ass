df = read_patients()
df = clean_race_eth(df)
age = df.dropna(subset=['AGE'])
n_sex = age['SEX'].nunique(); n_race = age['RACE'].nunique(); n_eth = age['ETHNICITY'].nunique()
fig = go.Figure()
fig.add_trace(go.Histogram(x=age['AGE'],nbinsx=40,name='All',marker_color='steelblue',opacity=0.75))
colors = px.colors.qualitative.Plotly
for i,cat in enumerate(age['SEX'].value_counts().index):
    d = age[age['SEX']==cat]
    fig.add_trace(go.Histogram(x=d['AGE'],nbinsx=40,name=f'Sex:{cat}',marker_color=colors[i],opacity=0.6,visible=False))
for i,cat in enumerate(age['RACE'].value_counts().index):
    d = age[age['RACE']==cat]
    fig.add_trace(go.Histogram(x=d['AGE'],nbinsx=40,name=f'Race:{cat}',marker_color=colors[i%len(colors)],opacity=0.6,visible=False))
for i,cat in enumerate(age['ETHNICITY'].value_counts().index):
    d = age[age['ETHNICITY']==cat]
    fig.add_trace(go.Histogram(x=d['AGE'],nbinsx=40,name=f'Eth:{cat}',marker_color=colors[i%len(colors)],opacity=0.6,visible=False))
tr = 1 + n_sex + n_race + n_eth
fig.update_layout(updatemenus=[dict(buttons=[
    dict(label='Overall',method='update',args=[{'visible':[True]+[False]*(tr-1)},{'title':'Age -- Overall'}]),
    dict(label='By SEX',method='update',args=[{'visible':[True]+[True]*n_sex+[False]*(n_race+n_eth)},{'title':'Age by SEX'}]),
    dict(label='By RACE',method='update',args=[{'visible':[True]+[False]*n_sex+[True]*n_race+[False]*n_eth},{'title':'Age by RACE'}]),
    dict(label='By ETHNICITY',method='update',args=[{'visible':[True]+[False]*n_sex+[False]*n_race+[True]*n_eth},{'title':'Age by ETHNICITY'}])
],direction='down',showactive=True,x=1.0,y=1.15)],title='Age Distribution',xaxis_title='Age (years)',yaxis_title='Count',height=500,bargap=0.05)
fig.show()
