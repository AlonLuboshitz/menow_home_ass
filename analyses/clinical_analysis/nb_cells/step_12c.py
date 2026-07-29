mg = merged.groupby('PATIENT_ID')['CANCER_GROUP'].nunique()
m = mg[mg>1]
print(f'Patients in multiple curves: {len(m)} (entries: {int(m.sum())})')
