ct = pd.crosstab(merged['CANCER_GROUP'],merged['SEX'])
chi2,p,dof,_ = chi2_contingency(ct.fillna(0))
print(f'Chi2: x2={chi2:.2f}, p={p:.4f}, df={dof}')
