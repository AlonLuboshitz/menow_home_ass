print('Top 10 subtypes:')
print(merged['MOLECULAR_SUBTYPE'].value_counts().head(10))
print('\nDiversity:',div)
