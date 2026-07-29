print(f'Strategy unique: {df["EXPERIMENT_STRATEGY"].nunique()}')
print(f'Library unique: {df["RNA_LIBRARY_SELECTION"].nunique()}')
print('Strategies:',list(df['EXPERIMENT_STRATEGY'].unique()))
