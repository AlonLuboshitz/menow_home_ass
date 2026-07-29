df = read_patients()
df = clean_os(df)
df = clean_efs(df)
os_c = df['os_label'].value_counts(dropna=False).reset_index()
os_c.columns = ['Status','Count']
os_c['Status'] = os_c['Status'].fillna('Unknown')
fig = px.pie(os_c,values='Count',names='Status',title='Overall Survival Status',hole=0.3)
fig.update_traces(textinfo='label+percent')
fig.show()
