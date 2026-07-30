print(f'TF missing: {df["TUMOR_FRACTION"].isna().sum()}/{len(df)}')
print(f'TP missing: {df["TUMOR_PLOIDY"].isna().sum()}/{len(df)}')
print(f'TF range: {tf.min():.3f} - {tf.max():.3f}')
print(f'TP range: {tp.min():.1f} - {tp.max():.1f}')
