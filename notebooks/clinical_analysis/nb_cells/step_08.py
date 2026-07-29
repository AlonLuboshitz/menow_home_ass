df = read_samples()
df = clean_subtype(df)
sc = df['MOLECULAR_SUBTYPE'].value_counts().head(20).reset_index()
sc.columns = ['Subtype','Count']
sc['%'] = (sc['Count']/len(df)*100).round(1)
fig = px.bar(sc,y='Subtype',x='Count',orientation='h',title='Top 20 Subtypes',
    text=[f'{c} ({p:.1f}%)' for c,p in zip(sc['Count'],sc['%'])],color='Count',color_continuous_scale='Viridis')
fig.update_traces(textposition='outside')
fig.update_layout(height=600,yaxis={'categoryorder':'total ascending'})
fig.show()
