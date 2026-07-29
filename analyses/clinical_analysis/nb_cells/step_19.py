patients = read_patients(); samples = read_samples()
patients = clean_os(patients); patients = clean_efs(patients); patients = clean_pred(patients)
samples = clean_subtype(samples)
n_p = patients['PATIENT_ID'].nunique(); n_s = samples['SAMPLE_ID'].nunique()
n_cg = samples['CANCER_GROUP'].nunique()
a_m = patients['AGE'].median(); a_r = (patients['AGE'].min(),patients['AGE'].max())
p_m = (patients['SEX']=='Male').mean()*100; p_f = (patients['SEX']=='Female').mean()*100
n_pr = (patients['CANCER_PREDISPOSITIONS']!='No predisposition').sum()
t_cg = samples['CANCER_GROUP'].value_counts().index[0]; t_cg_c = samples['CANCER_GROUP'].value_counts().iloc[0]
n_uc = (samples['MOLECULAR_SUBTYPE']=='Unclassified').sum(); n_mt = samples['TUMOR_FRACTION'].isna().sum(); n_mtp = samples['TUMOR_PLOIDY'].isna().sum()
os_c = patients[['OS_MONTHS','os_event']].dropna(); os_m = os_c['OS_MONTHS'].median()
report = f'''# PBTA_RNA Basic Clinical Summary

## Dataset Overview
- Patients: {n_p}
- Samples: {n_s}
- Cancer Groups: {n_cg}
- Most common: {t_cg} ({t_cg_c})

## Demographics
- Age: median {a_m:.1f}y, range {a_r[0]:.0f}-{a_r[1]:.0f}
- Sex: {p_m:.1f}% M, {p_f:.1f}% F

## Survival
- Median OS: {os_m:.1f}mo
- Patients with OS data: {len(os_c)}/{n_p}

## Predispositions
- Known predisposition: {n_pr}/{n_p} ({n_pr/n_p*100:.1f}%)

## Subtypes
- Unclassified: {n_uc}/{n_s} ({n_uc/n_s*100:.1f}%)
- Distinct subtypes: {samples['MOLECULAR_SUBTYPE'].nunique()}

## Tumor Purity
- Missing TF: {n_mt}/{n_s} ({n_mt/n_s*100:.1f}%)
- Missing TP: {n_mtp}/{n_s} ({n_mtp/n_s*100:.1f}%)
'''
with open('/home/alon/menow_home_ass/basic_clinical_summary.md','w') as f: f.write(report)
print('Saved: basic_clinical_summary.md')
print(report)
