v = df['AGE'].isna().sum()
print(f'AGE has {v} missing ({v/len(df):.1%}) -- excluded from age plots')
print(f'OS_STATUS null: {df["OS_STATUS"].isna().sum()}')
print(f'EFS_STATUS null: {df["EFS_STATUS"].isna().sum()}')
