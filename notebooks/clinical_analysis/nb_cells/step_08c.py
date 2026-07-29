unc = (df['MOLECULAR_SUBTYPE']=='Unclassified').sum()
print(f'Unclassified: {unc} ({unc/len(df):.1%})')
print(f'Distinct subtypes: {df["MOLECULAR_SUBTYPE"].nunique()}')
