miss = df.isna().mean().sort_values(ascending=False)*100
fig = px.bar(x=miss.index,y=miss.values,title='Sample Data -- Missingness',labels={'x':'Column','y':'% Missing'},
    text=[f'{v:.1f}%' for v in miss.values],color=miss.values,color_continuous_scale='Reds')
fig.update_traces(textposition='outside')
fig.update_layout(xaxis_tickangle=-45,height=500)
fig.show()
print('\nMissing > 0:'); print(miss[miss>0].round(1).to_string())
