# Step 12: Validation -- multi-group patients and per-group test coverage.
mg = merged.groupby('PATIENT_ID')['CANCER_GROUP'].nunique()
m = mg[mg>1]
print(f'Patients in multiple curves: {len(m)} (entries: {int(m.sum())})')
print(f'Groups: {len(os_n)} with OS data, {len(efs_n)} with EFS data')
print(f'Per-group tests (N>={MIN_N}): {len(test_groups_os)} OS, {len(test_groups_efs)} EFS | plotted: {len(plot_groups)}')
