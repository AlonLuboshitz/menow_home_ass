print(f'Multi-Group patients: {(pg["n"]>1).sum()} / {pg["PATIENT_ID"].nunique()}')
print(f'Max groups: {pg["n"].max()}')
if len(m2)>0: print(f'Exactly 2 groups: {len(m2)}')
