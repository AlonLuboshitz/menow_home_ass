multi = (sc['N']>1).sum()
print(f'Patients with >1 sample: {multi} ({multi/len(sc):.1%})')
print(f'Mean: {sc["N"].mean():.2f}')
