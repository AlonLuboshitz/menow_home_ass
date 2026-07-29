t = patients[patients['hp']==True]; f = patients[patients['hp']==False]
os_t = t[['OS_MONTHS','os_event']].dropna(); os_f = f[['OS_MONTHS','os_event']].dropna()
if len(os_t)>3 and len(os_f)>3: print(f'Log-rank OS: p={logrank2(os_t["OS_MONTHS"],os_t["os_event"],os_f["OS_MONTHS"],os_f["os_event"]):.4f}')
ef_t = t[['EFS_MONTHS','efs_event']].dropna(); ef_f = f[['EFS_MONTHS','efs_event']].dropna()
if len(ef_t)>3 and len(ef_f)>3: print(f'Log-rank EFS: p={logrank2(ef_t["EFS_MONTHS"],ef_t["efs_event"],ef_f["EFS_MONTHS"],ef_f["efs_event"]):.4f}')
at = t['AGE'].dropna(); af = f['AGE'].dropna()
if len(at)>3 and len(af)>3:
    s,p = mannwhitneyu(at,af,alternative='two-sided')
    print(f'MW Age: U={s:.0f}, p={p:.4f}, medians: {at.median():.1f} vs {af.median():.1f}')
