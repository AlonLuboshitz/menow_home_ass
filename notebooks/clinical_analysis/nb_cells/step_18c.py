print(f'CNS_REGION missing: {merged["CNS_REGION"].isna().sum()}')
print(merged['CNS_REGION'].value_counts().head(8))
