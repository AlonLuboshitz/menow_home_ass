pc = dex[~dex['pred'].isin(['No predisposition','Unknown'])]
pc = pc['pred'].value_counts().head(15).reset_index()
pc.columns = ['Pred','Count']
pc['%'] = (pc['Count']/total_p*100).round(1)
fig = px.bar(pc,y='Pred',x='Count',orientation='h',title=f'Top 15 Predispositions (N={total_p})',
             text=[f'{c} ({p:.1f}%)' for c,p in zip(pc['Count'],pc['%'])],color='Count',color_continuous_scale='Blues')
fig.update_traces(textposition='outside')
fig.update_layout(height=500,yaxis={'categoryorder':'total ascending'})
fig.show()
