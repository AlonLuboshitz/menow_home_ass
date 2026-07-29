miss = df.isna().mean().sort_values(ascending=False)*100
fig = px.bar(x=miss.index,y=miss.values,title='Patient Data -- Missingness',labels={'x':'Column','y':'% Missing'},text=[f'{v:.1f}%' for v in miss.values])
fig.update_traces(marker_color='crimson',textposition='outside')
fig.update_layout(xaxis_tickangle=-45,height=450)
fig.show()
