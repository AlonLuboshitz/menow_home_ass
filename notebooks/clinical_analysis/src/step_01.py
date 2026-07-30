df = read_patients()
print(f'Shape: {df.shape[0]} rows x {df.shape[1]} columns\n')
rows = []
for col in df.columns:
    nn = df[col].notna().sum(); nu = df[col].isna().sum(); pct = nu/len(df)*100
    if df[col].dtype == 'object': ex = f'unique={df[col].nunique()}'
    else: ex = f'min={df[col].min()}, max={df[col].max()}'
    rows.append({'Column':col,'Dtype':str(df[col].dtype),'Non-null':nn,'Null':nu,'% Null':f'{pct:.1f}%','Extra':ex})
print(pd.DataFrame(rows).to_string(index=False))
