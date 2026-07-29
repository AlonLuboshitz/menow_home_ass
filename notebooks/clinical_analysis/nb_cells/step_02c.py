print(f'RACE after: {sorted(df["RACE"].unique())}')
print(f'ETHNICITY after: {sorted(df["ETHNICITY"].unique())}')
print(f'Age data: {df["AGE"].notna().sum()}/{len(df)}')
