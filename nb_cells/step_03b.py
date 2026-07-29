df['efs_bin'] = df['efs_event'].map({1:'Event',0:'No Event'}).fillna('Unknown')
fig = make_subplots(rows=1,cols=2,subplot_titles=('EFS Binary','EFS Detailed'))
eb = df['efs_bin'].value_counts().reset_index(); eb.columns=['Status','Count']
fig.add_trace(go.Bar(x=eb['Status'],y=eb['Count'],marker_color=['lightcoral','lightgreen','lightgray'],text=eb['Count'],textposition='outside',showlegend=False),row=1,col=1)
ed = df['efs_detail'].value_counts(dropna=False).reset_index(); ed.columns=['Status','Count']
ed['Status'] = ed['Status'].fillna('Unknown')
fig.add_trace(go.Bar(x=ed['Status'],y=ed['Count'],marker_color='lightcoral',text=ed['Count'],textposition='outside',showlegend=False),row=1,col=2)
fig.update_layout(title='Event-Free Survival Status',height=450,xaxis2_tickangle=-45)
fig.show()
