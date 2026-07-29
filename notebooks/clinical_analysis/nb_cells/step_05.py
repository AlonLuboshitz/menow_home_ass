df = read_samples()
print(f'Shape: {df.shape[0]} rows x {df.shape[1]} columns\n')
rows = []
for col in df.columns:
    nn = df[col].notna().sum(); nu = df[col].isna().sum(); pct = nu/len(df)*100
    if df[col].dtype=='object': ex = f'unique={df[col].nunique()}'
    elif df[col].dtype in ('float64','int64'): ex = f'min={df[col].min():.2f}, max={df[col].max():.2f}'
    else: ex = ''
    rows.append({'Column':col,'Dtype':str(df[col].dtype),'Non-null':nn,'Null':nu,'% Null':f'{pct:.1f}%','Extra':ex})
print(pd.DataFrame(rows).to_string(index=False))
