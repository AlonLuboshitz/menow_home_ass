print(f'Most common CG: {df["CANCER_GROUP"].value_counts().index[0]} ({df["CANCER_GROUP"].value_counts().iloc[0]})')
print(f'Total CGs: {df["CANCER_GROUP"].nunique()}')
