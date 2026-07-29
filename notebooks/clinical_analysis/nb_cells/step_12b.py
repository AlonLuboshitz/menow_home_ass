rows = []
for cg in top6:
    s = merged[merged['CANCER_GROUP']==cg]
    os_s = s[['OS_MONTHS','os_event']].dropna()
    ef_s = s[['EFS_MONTHS','efs_event']].dropna()
    rows.append({'CG':cg,'N':len(s),'N_Pat':s['PATIENT_ID'].nunique(),
        'MedOS':f'{os_s["OS_MONTHS"].median():.1f}' if len(os_s)>0 else 'N/A',
        'MedEFS':f'{ef_s["EFS_MONTHS"].median():.1f}' if len(ef_s)>0 else 'N/A'})
print(pd.DataFrame(rows).to_string(index=False))
